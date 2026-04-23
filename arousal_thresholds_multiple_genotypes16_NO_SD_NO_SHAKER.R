# =============================================================================
# Arousal Threshold Analysis — RESTED vs SLEEPING ONLY
# MAXIMUM STATISTICAL RIGOUR VERSION
#
# MULTIPLE COMPARISONS PHILOSOPHY — ERR ON SIDE OF CAUTION
# ─────────────────────────────────────────────────────────────────────────────
# This script deliberately overcorrects rather than undercorrects.
# Every family of tests is defined as broadly as defensible.
#
# CORRECTION STRUCTURE:
#
#   Analysis 1 (rested vs sleeping within each genotype):
#     → 4 pre-specified tests (one per genotype), corrected as ONE family
#     → Holm PRIMARY | BH SENSITIVITY
#     → No omnibus gate needed (one comparison per genotype = no inflation risk
#       within genotype; Holm handles the cross-genotype family)
#
#   Analysis 2 (sleep latency across genotypes):
#     → Kruskal-Wallis omnibus FIRST; pairwise only if KW significant
#     → Holm on all pairwise Mann-Whitney U tests
#
#   Analysis 3 (threshold during sleep, across genotypes):
#     → Omnibus (KW or ANOVA) FIRST; pairwise only if omnibus significant
#     → Holm on all pairwise tests
#
#   Analyses 5 + 10 (latency: HK/canton_s and inc/w1118):
#     → Treated as ONE family of 2 tests (both ask: does mutant latency differ?)
#     → Holm correction across both raw p-values
#
#   Analyses 7 + 11 (threshold during sleep: HK/canton_s and inc/w1118):
#     → Treated as ONE family of 2 tests (both ask: does mutant threshold differ?)
#     → Holm correction across both raw p-values
#
#   Analysis 12 (genotype x sleep condition factorial):
#     → Omnibus interaction test FIRST (ART-ANOVA or two-way ANOVA)
#     → Pairwise genotype comparisons within each condition ONLY if
#       interaction is significant (strict gatekeeping)
#     → If interaction not significant: report main effect of genotype only
#     → Holm correction across conditions (rested + sleeping) as a family
#
#   Analysis 6 (chi-square waking proportion):
#     → Single pre-specified test; Fisher's exact as sensitivity
#     → No correction needed (one test)
#
#   Analysis 8 (censored latency HK vs canton_s):
#     → Included in Analyses 5+10 family if testing same biological question,
#       otherwise standalone (censored data makes it a separate methodological
#       approach to the same question as Analysis 5)
#     → Treated as sensitivity check for Analysis 5; no additional correction
#
# ASSUMPTION CHECKS (printed for every test):
#   Shapiro-Wilk per group + Levene's across groups
#   → Both violated: Mann-Whitney U (non-parametric)
#   → Normality violated, variances equal, large n: note CLT but use MWU
#   → All met: Welch t-test (parametric)
#
# TEST LABELS ON PLOTS:
#   Black brackets = Holm-corrected PRIMARY
#   Blue  brackets = BH-corrected SENSITIVITY (where BH sig but Holm not)
#   Panel subtitle = test used (Mann-Whitney U or Welch t)
# =============================================================================

library(car)
library(emmeans)

if (!requireNamespace("ARTool", quietly = TRUE)) {
  install.packages("ARTool", repos = "https://cloud.r-project.org")
}
library(ARTool)

# =============================================================================
# GENOTYPE ORDER (shaker removed)
# =============================================================================
GENOTYPE_ORDER <- c("canton_s", "w1118", "hk", "inc")

order_genos <- function(f) {
  present <- intersect(GENOTYPE_ORDER, levels(droplevels(factor(f))))
  extra   <- setdiff(levels(droplevels(factor(f))), GENOTYPE_ORDER)
  factor(f, levels = c(present, extra))
}

# =============================================================================
# DATA LOADING
# =============================================================================
AWAKE_FILE    <- "/home/joeh/mutantscreen_round_two/rested_vs_sd_all_genotypes.csv"
SLEEPING_FILE <- "/home/joeh/mutantscreen_round_two/sleeping_all_genotypes.csv"
awake        <- read.csv(AWAKE_FILE,    stringsAsFactors = FALSE)
sleeping_raw <- read.csv(SLEEPING_FILE, stringsAsFactors = FALSE)

# =============================================================================
# ZT FROM FILE CTIME
# =============================================================================
ZT_ZERO_HOUR <- 8L
ctime_to_ZT <- function(filepath) {
  filepath <- trimws(as.character(filepath))
  info     <- file.info(filepath)
  ct       <- info$ctime
  if (is.na(ct)) { warning(paste("No ctime for:", filepath)); return(NA_integer_) }
  hr <- as.integer(format(ct, "%H")) +
    as.numeric(format(ct, "%M")) / 60 +
    as.numeric(format(ct, "%S")) / 3600
  as.integer(hr - ZT_ZERO_HOUR) %% 24L
}
get_ZT_col <- function(file_vec) vapply(file_vec, ctime_to_ZT, integer(1L))
cat("Computing ZT from file ctimes...\n")
awake$ZT_fr        <- get_ZT_col(awake$file)
sleeping_raw$ZT_fr <- get_ZT_col(sleeping_raw$file)
cat("Done.\n\n")

# =============================================================================
# BUILD COMBINED DATA FRAME — SD AND SHAKER EXCLUDED
# =============================================================================
sleeping_export <- sleeping_raw[
  !(is.na(sleeping_raw$genotype) | sleeping_raw$genotype == "") &
    !is.na(sleeping_raw$volume_percent_threshold) &
    trimws(tolower(sleeping_raw$genotype)) != "shaker", ]
sleeping_export$sleep_cond <- "sleeping"
if (!"volume_percent_2sd_200ms" %in% names(sleeping_export)) sleeping_export$volume_percent_2sd_200ms <- NA_real_
if (!"volume_percent_2sd_500ms" %in% names(sleeping_export)) sleeping_export$volume_percent_2sd_500ms <- NA_real_
if (!"volume_percent_2sd_1s"    %in% names(sleeping_export)) sleeping_export$volume_percent_2sd_1s    <- NA_real_
sleeping_export <- sleeping_export[, c("genotype","sleep_cond","volume_percent_threshold",
                                       "sleep_latency","day","volume_percent_2sd_200ms",
                                       "volume_percent_2sd_500ms","volume_percent_2sd_1s","ZT_fr")]

# Rested only from awake file, shaker excluded
awake_export <- awake[!is.na(awake$volume_percent_threshold) &
                        trimws(tolower(awake$sleep_cond)) == "rested" &
                        trimws(tolower(awake$genotype)) != "shaker",
                      c("sleep_cond","genotype","volume_percent_threshold",
                        "volume_percent_2sd_200ms","volume_percent_2sd_500ms",
                        "volume_percent_2sd_1s","ZT_fr")]
awake_export$sleep_latency <- NA_real_
awake_export$day           <- NA_real_
awake_export <- awake_export[, c("genotype","sleep_cond","volume_percent_threshold",
                                 "sleep_latency","day","volume_percent_2sd_200ms",
                                 "volume_percent_2sd_500ms","volume_percent_2sd_1s","ZT_fr")]

combined            <- rbind(sleeping_export, awake_export)
combined$genotype   <- order_genos(factor(trimws(combined$genotype)))
combined$sleep_cond <- factor(trimws(combined$sleep_cond),
                              levels = c("rested","sleeping"))

stopifnot(!any(combined$sleep_cond == "sd", na.rm = TRUE))
stopifnot(!any(trimws(tolower(as.character(combined$genotype))) == "shaker", na.rm = TRUE))
cat("SD excluded confirmed.\n")
cat("Shaker excluded confirmed.\n")

write.csv(combined, "arousal_latency_ZT_rested_sleeping_only.csv", row.names = FALSE)
cat("CSV written: arousal_latency_ZT_rested_sleeping_only.csv\n")
cat("  Rows:", nrow(combined), "\n\n")

# =============================================================================
# COLOUR PALETTE
# =============================================================================
PALETTE_BORDER <- c("#C0392B","#4A7FC1","#27AE60","#8E44AD",
                    "#E67E22","#16A085","#2C3E50","#F39C12")
PALETTE_FILL <- adjustcolor(PALETTE_BORDER, alpha.f = 0.55)
make_geno_colours <- function(geno_levels) {
  n <- length(geno_levels)
  if (n > length(PALETTE_BORDER))
    stop(paste("Add more colours to PALETTE_BORDER — need at least", n))
  border <- setNames(PALETTE_BORDER[seq_len(n)], geno_levels)
  fill   <- setNames(PALETTE_FILL[seq_len(n)],   geno_levels)
  list(border = border, fill = fill)
}

# =============================================================================
# BUILD ANALYSIS DATA FRAMES
# =============================================================================
df <- combined[!is.na(combined$volume_percent_threshold),
               c("genotype","sleep_cond","volume_percent_threshold","day","ZT_fr")]
stopifnot(!any(df$sleep_cond == "sd", na.rm = TRUE))
stopifnot(!any(trimws(tolower(as.character(df$genotype))) == "shaker", na.rm = TRUE))

remove_iqr_outliers <- function(x) {
  Q1 <- quantile(x,0.25,na.rm=TRUE); Q3 <- quantile(x,0.75,na.rm=TRUE)
  x[x >= (Q1-1.5*(Q3-Q1)) & x <= (Q3+1.5*(Q3-Q1))]
}
df_clean <- do.call(rbind, lapply(
  split(df, list(df$genotype,df$sleep_cond), drop=TRUE), function(sub) {
    if (nrow(sub)==0) return(NULL)
    sub[sub$volume_percent_threshold %in%
          remove_iqr_outliers(sub$volume_percent_threshold), ]
  }
))
df_clean$sleep_cond <- factor(df_clean$sleep_cond, levels=c("rested","sleeping"))
df_clean$genotype   <- order_genos(df_clean$genotype)

latency <- combined[
  !(is.na(combined$genotype)|combined$genotype=="") & !is.na(combined$sleep_latency),
  c("genotype","sleep_latency","volume_percent_threshold","day",
    "volume_percent_2sd_200ms","volume_percent_2sd_500ms","volume_percent_2sd_1s","ZT_fr")]
latency$genotype <- order_genos(factor(trimws(as.character(latency$genotype))))

remove_iqr <- function(data, col, group_col) {
  keep <- rep(TRUE, nrow(data))
  for (g in unique(data[[group_col]])) {
    idx  <- which(data[[group_col]]==g)
    vals <- data[[col]][idx]
    q1 <- quantile(vals,0.25,na.rm=TRUE); q3 <- quantile(vals,0.75,na.rm=TRUE)
    iqr <- q3-q1
    keep[idx[vals < (q1-1.5*iqr)|vals > (q3+1.5*iqr)]] <- FALSE
  }
  data[keep,]
}
latency_clean <- remove_iqr(latency,"sleep_latency","genotype")
latency_clean$genotype <- order_genos(latency_clean$genotype)

cat("Data loaded.\n")
cat("df:", nrow(df), "rows | Genotypes:", paste(levels(df$genotype),collapse=", "), "\n")
cat("Sleep conditions:", paste(levels(df$sleep_cond),collapse=", "), "\n")
cat("latency:", nrow(latency), "rows\n\n")

# =============================================================================
# COLOURS FOR CONDITIONS
# =============================================================================
cond_cols      <- c(rested="#1A5CB0", sleeping="#E00000")
cond_cols_fill <- c(rested=adjustcolor("#1A5CB0",alpha.f=0.55),
                    sleeping=adjustcolor("#E00000",alpha.f=0.55))
cond_labels    <- c(rested="Awake", sleeping="Asleep")

par_defaults <- list(cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
                     font=2,font.main=2,font.lab=2,font.axis=2)
do.call(par, par_defaults)

# =============================================================================
# HELPERS
# =============================================================================
sig_label <- function(p) {
  if (is.na(p)||p>=0.05) return(NULL)
  if (p<0.001) return("***")
  if (p<0.01)  return("**")
  return("*")
}

draw_manual_box <- function(xi, vals, fill_col, dot_col, box_width=0.35) {
  vals <- vals[!is.na(vals)]
  if (length(vals)==0) return()
  q1 <- quantile(vals,0.25); med <- quantile(vals,0.50); q3 <- quantile(vals,0.75)
  iqr <- q3-q1
  lo <- max(min(vals),q1-1.5*iqr); hi <- min(max(vals),q3+1.5*iqr)
  segments(xi,lo,xi,q1,lwd=2,col="black"); segments(xi,q3,xi,hi,lwd=2,col="black")
  segments(xi-box_width*0.5,lo,xi+box_width*0.5,lo,lwd=2,col="black")
  segments(xi-box_width*0.5,hi,xi+box_width*0.5,hi,lwd=2,col="black")
  rect(xi-box_width,q1,xi+box_width,q3,col=fill_col,border=NA)
  rect(xi-box_width,q1,xi+box_width,q3,col=NA,border="black",lwd=2)
  segments(xi-box_width,med,xi+box_width,med,lwd=3,col="black")
  points(jitter(rep(xi,length(vals)),amount=box_width*0.35), vals,
         pch=16, col=adjustcolor(dot_col,alpha.f=0.65), cex=1.0)
}

draw_bracket <- function(x1,x2,y,label,tip=1.5,col="black") {
  segments(x1,y,x2,y,lwd=1.5,col=col)
  segments(x1,y-tip,x1,y,lwd=1.5,col=col)
  segments(x2,y-tip,x2,y,lwd=1.5,col=col)
  text((x1+x2)/2,y+tip*0.8,label,cex=1.3,font=2,col=col)
}

two_group_bracket <- function(ylim,x1,x2,p,col="black") {
  if (is.na(p)||p>=0.05) return(invisible(NULL))
  lbl <- if(p<0.001)"***" else if(p<0.01)"**" else "*"
  brk_y <- ylim[2]*0.88; tip <- (ylim[2]-ylim[1])*0.015
  segments(x1,brk_y,x2,brk_y,lwd=1.5,col=col)
  segments(x1,brk_y-tip,x1,brk_y,lwd=1.5,col=col)
  segments(x2,brk_y-tip,x2,brk_y,lwd=1.5,col=col)
  text((x1+x2)/2,brk_y+tip*2,lbl,cex=1.4,font=2,col=col)
}

p_footer <- function(p1,p2,label1="Test 1",label2="Test 2") {
  s1 <- formatC(p1,digits=3,format="g"); s2 <- formatC(p2,digits=3,format="g")
  mtext(paste0(label1," p=",s1,"  |  ",label2," p=",s2),
        side=1,line=3.5,cex=0.95,col="grey40",font=2)
}

cohens_d <- function(x1,x2) {
  n1 <- length(x1); n2 <- length(x2)
  pooled_sd <- sqrt(((n1-1)*var(x1)+(n2-1)*var(x2))/(n1+n2-2))
  (mean(x2)-mean(x1))/pooled_sd
}

cohens_d_ci <- function(x1,x2,conf=0.95) {
  n1 <- length(x1); n2 <- length(x2)
  d  <- cohens_d(x1,x2)
  se <- sqrt((n1+n2)/(n1*n2)+d^2/(2*(n1+n2-2)))
  z  <- qnorm(1-(1-conf)/2)
  c(d=d, lower=round(d-z*se,3), upper=round(d+z*se,3))
}

# =============================================================================
# ASSUMPTION CHECK HELPER
# Prints SW per group and Levene across groups.
# Returns parametric=TRUE only if ALL groups normal AND variances homogeneous.
# =============================================================================
check_assumptions <- function(vals, groups, context_label="") {
  groups  <- droplevels(factor(groups))
  glevels <- levels(groups)
  cat(paste0("\n  --- Assumption checks",
             if(nchar(context_label)>0) paste0(" [",context_label,"]"), " ---\n"))
  
  sw_res <- do.call(rbind, lapply(glevels, function(g) {
    v <- vals[groups==g & !is.na(vals)]
    if (length(v)<3) {
      cat(sprintf("    SW [%s]: n=%d < 3 — untestable, assume non-normal\n",g,length(v)))
      return(data.frame(group=g,n=length(v),SW_W=NA,SW_p=NA,normal=FALSE,
                        stringsAsFactors=FALSE))
    }
    sw <- shapiro.test(v); ok <- sw$p.value>=0.05
    cat(sprintf("    SW [%s]: n=%d, W=%.4f, p=%.4f  → %s\n",
                g,length(v),sw$statistic,sw$p.value,
                if(ok)"normal OK" else "*** NON-NORMAL ***"))
    data.frame(group=g,n=length(v),SW_W=sw$statistic,SW_p=sw$p.value,normal=ok,
               stringsAsFactors=FALSE)
  }))
  
  all_normal <- all(sw_res$normal, na.rm=TRUE)
  
  valid_groups <- sw_res$group[sw_res$n>=2]
  if (length(valid_groups)>=2) {
    lev_df   <- data.frame(v=vals[!is.na(vals)],g=groups[!is.na(vals)])
    lev_df   <- lev_df[lev_df$g %in% valid_groups,]
    lev_df$g <- droplevels(factor(lev_df$g))
    lev      <- car::leveneTest(v~g, data=lev_df, center=median)
    lev_p    <- lev[["Pr(>F)"]][1]; lev_ok <- lev_p>=0.05
    cat(sprintf("    Levene's: F(%d,%d)=%.4f, p=%.4f  → %s\n",
                lev[["Df"]][1],lev[["Df"]][2],lev[["F value"]][1],lev_p,
                if(lev_ok)"homogeneous OK" else "*** HETEROGENEOUS ***"))
  } else { lev_p <- NA; lev_ok <- TRUE
  cat("    Levene's: insufficient groups\n") }
  
  use_param <- all_normal && lev_ok
  cat(sprintf("    → Decision: %s\n",
              if(use_param)"PARAMETRIC (Welch t / ANOVA)"
              else         "NON-PARAMETRIC (Mann-Whitney U / Kruskal-Wallis)"))
  list(sw=sw_res, levene_p=lev_p, normal_ok=all_normal,
       homog_ok=lev_ok, parametric=use_param)
}

# =============================================================================
# RAW TEST HELPER — runs both tests, selects primary by assumption check,
# returns raw p-values WITHOUT applying correction (correction applied
# externally across the appropriate family).
# =============================================================================
raw_two_group_test <- function(vals1, vals2, label1, label2) {
  chk <- check_assumptions(c(vals1,vals2),
                           factor(c(rep(label1,length(vals1)),
                                    rep(label2,length(vals2)))),
                           paste0(label1," vs ",label2))
  tt <- t.test(vals1, vals2)
  wc <- wilcox.test(vals1, vals2, exact=FALSE)
  
  use_wc <- !chk$parametric
  cat(sprintf("    → Primary test: %s\n",
              if(use_wc)"Mann-Whitney U" else "Welch t-test"))
  cat(sprintf("    Welch t: t=%.3f, df=%.1f, raw p=%.4f\n",
              tt$statistic, tt$parameter, tt$p.value))
  cat(sprintf("    Mann-Whitney U: W=%.1f, raw p=%.4f\n",
              wc$statistic, wc$p.value))
  
  list(tt=tt, wc=wc, use_wilcoxon=use_wc,
       primary_p_raw=if(use_wc) wc$p.value else tt$p.value,
       sens_p_raw=if(use_wc) tt$p.value else wc$p.value,
       primary_label=if(use_wc)"Mann-Whitney U [PRIMARY]" else "Welch t [PRIMARY]",
       sens_label=if(use_wc)"Welch t" else "Mann-Whitney U")
}

# =============================================================================
# PLOT HELPER for a single two-group comparison with pre-computed
# Holm-corrected and BH-corrected p-values
# =============================================================================
plot_two_group <- function(vals1, vals2, label1, label2,
                           p_holm, p_bh, use_wilcoxon,
                           ylab="Value", main_title="",
                           ylim_override=NULL) {
  ylim_use <- if(!is.null(ylim_override)) ylim_override else
    c(0, max(c(vals1,vals2),na.rm=TRUE)*1.4)
  gc2 <- make_geno_colours(c(label1,label2))
  
  par(mfrow=c(1,1), mar=c(6,6,5,2),
      cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
      font=2,font.main=2,font.lab=2,font.axis=2)
  plot(0,0,type="n", xlim=c(0.5,2.5), ylim=ylim_use,
       xlab="Genotype", ylab=ylab, main=main_title, xaxt="n", las=1, bty="n")
  axis(1,at=1:2,labels=c(label1,label2),tick=FALSE,cex.axis=1.3,font=2)
  mtext(if(use_wilcoxon)"Mann-Whitney U primary" else "Welch t primary",
        side=3,line=0,cex=0.9,
        col=if(use_wilcoxon)"darkred" else "grey30",font=3)
  draw_manual_box(1,vals1,gc2$fill[label1],gc2$border[label1])
  draw_manual_box(2,vals2,gc2$fill[label2],gc2$border[label2])
  
  # Holm bracket (black) — primary
  lbl_h <- sig_label(p_holm)
  if (!is.null(lbl_h)) {
    par(xpd=NA)
    two_group_bracket(ylim_use,1,2,p_holm,col="black")
    par(xpd=FALSE)
  }
  # BH bracket (blue) — only if BH sig but Holm not
  lbl_b <- sig_label(p_bh)
  if (!is.null(lbl_b) && is.null(lbl_h)) {
    par(xpd=NA)
    brk_y <- ylim_use[2]*0.94; tip <- (ylim_use[2]-ylim_use[1])*0.015
    segments(1,brk_y,2,brk_y,lwd=1.5,col="#2980B9")
    segments(1,brk_y-tip,1,brk_y,lwd=1.5,col="#2980B9")
    segments(2,brk_y-tip,2,brk_y,lwd=1.5,col="#2980B9")
    text(1.5,brk_y+tip*2,paste0(lbl_b,"(BH)"),cex=1.3,font=2,col="#2980B9")
    par(xpd=FALSE)
  }
  legend("topright",legend=c(label1,label2),
         fill=gc2$fill[c(label1,label2)],border="black",
         bty="n",cex=1.2,text.font=2)
}

# =============================================================================
# ANALYSIS 1: Rested vs sleeping within each genotype
# ONE family of 4 tests (one per genotype, shaker excluded) — Holm PRIMARY
# No omnibus gate needed: one comparison per genotype
# =============================================================================
run_analysis1 <- function(d, label) {
  cat(paste0("\n", strrep("=",70), "\n"))
  cat(paste0("ANALYSIS 1 — Awake vs Sleeping per Genotype  [", label, "]\n"))
  cat("CORRECTION: Holm across all 4 genotypes as ONE family (PRIMARY)\n")
  cat("            BH across all 4 genotypes (SENSITIVITY)\n")
  cat(paste0(strrep("=",70), "\n"))
  
  d$genotype   <- order_genos(droplevels(d$genotype))
  d$sleep_cond <- droplevels(factor(d$sleep_cond, levels=c("rested","sleeping")))
  genos <- levels(d$genotype)
  
  # Collect raw p-values and test metadata for all genotypes
  raw_ps     <- setNames(rep(NA_real_, length(genos)), genos)
  use_wc_vec <- setNames(rep(NA, length(genos)), genos)
  r_list     <- list(); s_list <- list()
  
  for (g in genos) {
    sub <- d[d$genotype==g & !is.na(d$volume_percent_threshold), ]
    r_vals <- sub$volume_percent_threshold[sub$sleep_cond=="rested"]
    s_vals <- sub$volume_percent_threshold[sub$sleep_cond=="sleeping"]
    r_vals <- r_vals[!is.na(r_vals)]; s_vals <- s_vals[!is.na(s_vals)]
    
    cat(paste0("\n", strrep("-",60), "\n"))
    cat(sprintf("Genotype: %s  (awake n=%d, sleeping n=%d)\n",
                g, length(r_vals), length(s_vals)))
    
    if (length(r_vals)<2 || length(s_vals)<2) {
      cat("  SKIPPED — insufficient data.\n"); next
    }
    
    r_list[[g]] <- r_vals; s_list[[g]] <- s_vals
    res <- raw_two_group_test(r_vals, s_vals, "rested", "sleeping")
    raw_ps[g]     <- res$primary_p_raw
    use_wc_vec[g] <- res$use_wilcoxon
  }
  
  # Apply corrections across the full family of 4 genotypes
  valid    <- !is.na(raw_ps)
  p_holm   <- setNames(rep(NA_real_,length(genos)), genos)
  p_bh     <- setNames(rep(NA_real_,length(genos)), genos)
  p_holm[valid] <- p.adjust(raw_ps[valid], method="holm")
  p_bh[valid]   <- p.adjust(raw_ps[valid], method="BH")
  
  cat(paste0("\n", strrep("=",70), "\n"))
  cat("ANALYSIS 1 — CORRECTION SUMMARY (family of 4)\n")
  cat(sprintf("  %-12s  %-8s  %-10s  %-10s  %-10s  %s\n",
              "Genotype","Test","Raw p","Holm p","BH p","Sig (Holm)"))
  for (g in genos) {
    if (is.na(raw_ps[g])) next
    cat(sprintf("  %-12s  %-8s  %-10.4f  %-10.4f  %-10.4f  %s\n",
                g,
                if(isTRUE(use_wc_vec[g]))"MWU" else "Welch",
                raw_ps[g], p_holm[g], p_bh[g],
                if(!is.na(p_holm[g])&&p_holm[g]<0.05)"*" else "ns"))
  }
  
  # Cohen's d for each genotype
  cat("\nEffect sizes (Cohen's d, awake vs sleeping):\n")
  for (g in genos) {
    if (is.null(r_list[[g]])) next
    res_ci <- cohens_d_ci(r_list[[g]], s_list[[g]])
    d_mag  <- ifelse(abs(res_ci["d"])>=0.8,"large",
                     ifelse(abs(res_ci["d"])>=0.5,"medium","small"))
    cat(sprintf("  %s: d=%.3f [%.3f, %.3f] (%s)\n",
                g, res_ci["d"], res_ci["lower"], res_ci["upper"], d_mag))
  }
  
  # PLOT
  n_genos <- length(genos)
  par(mfrow=c(1,n_genos), mar=c(9,5.5,6,1), oma=c(0,0,7,0),
      cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
      font=2,font.main=2,font.lab=2,font.axis=2)
  
  for (g in genos) {
    r_vals <- r_list[[g]]; s_vals <- s_list[[g]]
    is_wc  <- isTRUE(use_wc_vec[g])
    
    plot(0,0,type="n", xlim=c(0.5,2.5), ylim=c(0,100),
         xlab="", ylab=if(g==genos[1])"Volume % Threshold" else "",
         main=g, xaxt="n", las=1, bty="n")
    axis(1,at=1:2,labels=FALSE,tick=FALSE)
    par(xpd=NA)
    text(1:2,y=-12,labels=c("Awake","Sleeping"),
         srt=45,adj=1,cex=1.1,font=2)
    mtext(if(is_wc)"Mann-Whitney U" else "Welch t",
          side=3,line=0,cex=0.85,
          col=if(is_wc)"darkred" else "grey30",font=3)
    par(xpd=FALSE)
    
    if (!is.null(r_vals) && !is.null(s_vals)) {
      draw_manual_box(1,r_vals,cond_cols_fill["rested"],  cond_cols["rested"])
      draw_manual_box(2,s_vals,cond_cols_fill["sleeping"],cond_cols["sleeping"])
    }
    
    # Holm bracket (black)
    lbl_h <- sig_label(p_holm[g])
    if (!is.null(lbl_h)) {
      par(xpd=NA)
      draw_bracket(1,2,y=92,label=lbl_h,tip=1.5,col="black")
      par(xpd=FALSE)
    }
    # BH bracket (blue) only where BH sig but Holm not
    lbl_b <- sig_label(p_bh[g])
    if (!is.null(lbl_b) && is.null(lbl_h)) {
      par(xpd=NA)
      draw_bracket(1,2,y=98,label=paste0(lbl_b,"(BH)"),tip=1.5,col="#2980B9")
      par(xpd=FALSE)
    }
    
    # Cohen's d footnote
    if (!is.null(r_vals) && !is.null(s_vals)) {
      d_val  <- round(cohens_d(r_vals,s_vals),3)
      d_size <- ifelse(abs(d_val)>=0.8,"large",ifelse(abs(d_val)>=0.5,"medium","small"))
      mtext(paste0("d=",d_val," (",d_size,")"),
            side=1,line=6,cex=0.9,col="grey25",font=3)
    }
  }
  
  mtext(paste0("Awake vs Sleeping: Arousal Threshold  [",label,"]\n",
               "Holm correction across all 4 genotypes (PRIMARY) | Blue=BH only"),
        outer=TRUE,cex=1.1,font=2,line=4)
  par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  legend("topright",
         legend=c("Awake","Asleep","Holm sig. (black)","BH only (blue)"),
         fill=c(cond_cols_fill["rested"],cond_cols_fill["sleeping"],NA,NA),
         border=c("black","black",NA,NA),
         lty=c(NA,NA,1,1),col=c(NA,NA,"black","#2980B9"),lwd=c(NA,NA,2,2),
         bty="n",cex=1.1,text.font=2)
  
  invisible(list(raw_ps=raw_ps, p_holm=p_holm, p_bh=p_bh,
                 use_wc=use_wc_vec, r_list=r_list, s_list=s_list))
}

datasets_main <- list("All data"=df, "IQR outliers removed"=df_clean)
for (label in names(datasets_main)) run_analysis1(datasets_main[[label]], label)

# =============================================================================
# ANALYSIS 2: Sleep Latency ~ genotype
# Kruskal-Wallis omnibus FIRST; pairwise Mann-Whitney U only if KW significant
# Holm PRIMARY across all pairwise tests | BH SENSITIVITY
# =============================================================================
run_latency_plot <- function(d, label) {
  d$genotype  <- order_genos(droplevels(d$genotype))
  geno_levels <- levels(d$genotype)
  n_genos     <- length(geno_levels)
  gc          <- make_geno_colours(geno_levels)
  
  cat(paste0("\n", strrep("=",70), "\n"))
  cat(paste0("ANALYSIS 2 — Sleep Latency  [", label, "]\n"))
  cat("Kruskal-Wallis omnibus gates pairwise Mann-Whitney U tests\n")
  cat(paste0(strrep("=",70), "\n"))
  cat("N:\n"); print(table(d$genotype))
  
  cat("\nSW normality checks per genotype (latency always non-parametric):\n")
  for (g in geno_levels) {
    v <- d$sleep_latency[d$genotype==g & !is.na(d$sleep_latency)]
    if (length(v)>=3) {
      sw <- shapiro.test(v)
      cat(sprintf("  SW [%s]: n=%d, W=%.4f, p=%.4f  → %s\n",
                  g,length(v),sw$statistic,sw$p.value,
                  if(sw$p.value>=0.05)"normal OK" else "*** NON-NORMAL ***"))
    } else cat(sprintf("  SW [%s]: n=%d < 3\n",g,length(v)))
  }
  cat("  → Non-parametric tests used regardless (latency is right-skewed)\n")
  
  if (n_genos==2) {
    g1 <- d$sleep_latency[d$genotype==geno_levels[1]]
    g2 <- d$sleep_latency[d$genotype==geno_levels[2]]
    wc <- wilcox.test(g1,g2,exact=FALSE); tt <- t.test(g1,g2)
    cat("\nMann-Whitney U [PRIMARY]:\n"); print(wc)
    cat("\nWelch t (sensitivity):\n"); print(tt)
    ylim_l <- c(0,max(d$sleep_latency,na.rm=TRUE)*1.3)
    par(mar=c(6,6,5,2),mfrow=c(1,1),
        cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
        font=2,font.main=2,font.lab=2,font.axis=2)
    plot(0,0,type="n",xlim=c(0.5,2.5),ylim=ylim_l,
         xlab="Genotype",ylab="Sleep Latency (s)",
         main=paste0("Sleep Latency  [",label,"]"),xaxt="n",las=1,bty="n")
    axis(1,at=1:2,labels=geno_levels,tick=FALSE,cex.axis=1.3,font=2)
    for (i in 1:2)
      draw_manual_box(i,d$sleep_latency[d$genotype==geno_levels[i]],
                      gc$fill[geno_levels[i]],gc$border[geno_levels[i]])
    two_group_bracket(ylim_l,1,2,wc$p.value)
    p_footer(wc$p.value,tt$p.value,"Mann-Whitney U [PRIMARY]","Welch t")
    
  } else {
    kw <- kruskal.test(sleep_latency~genotype, data=d)
    cat("\nKruskal-Wallis omnibus:\n"); print(kw)
    
    if (kw$p.value >= 0.05) {
      cat("\n  KW not significant (p=",round(kw$p.value,4),
          ") — pairwise tests NOT run (omnibus gate).\n")
      cat("  Conclusion: no evidence of latency differences across genotypes.\n")
    } else {
      cat("\n  KW significant — proceeding to pairwise Mann-Whitney U tests\n")
      cat("  Holm correction across ALL pairwise tests (PRIMARY)\n")
      pw_holm <- pairwise.wilcox.test(d$sleep_latency,d$genotype,
                                      p.adjust.method="holm",exact=FALSE)
      pw_bh   <- pairwise.wilcox.test(d$sleep_latency,d$genotype,
                                      p.adjust.method="BH",  exact=FALSE)
      cat("\nPairwise Mann-Whitney U — Holm [PRIMARY]:\n");   print(pw_holm)
      cat("\nPairwise Mann-Whitney U — BH [SENSITIVITY]:\n"); print(pw_bh)
      
      ylim_l <- c(0,max(d$sleep_latency,na.rm=TRUE)*1.5)
      par(mar=c(7,6,5,2),mfrow=c(1,1),
          cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
          font=2,font.main=2,font.lab=2,font.axis=2)
      xs <- seq_len(n_genos)
      plot(0,0,type="n",xlim=c(0.5,n_genos+0.5),ylim=ylim_l,
           xlab="",ylab="Sleep Latency (s)",
           main=paste0("Sleep Latency  [",label,"]"),xaxt="n",las=1,bty="n")
      axis(1,at=xs,labels=geno_levels,tick=FALSE,cex.axis=1.2,las=2,font=2)
      for (i in xs)
        draw_manual_box(i,d$sleep_latency[d$genotype==geno_levels[i]],
                        gc$fill[geno_levels[i]],gc$border[geno_levels[i]])
      pw_mat <- pw_holm$p.value; pw_bh_mat <- pw_bh$p.value
      bracket_step <- (ylim_l[2]-ylim_l[1])*0.07; drawn <- 0
      for (r in seq_len(nrow(pw_mat))) for (cc in seq_len(ncol(pw_mat))) {
        p_val <- pw_mat[r,cc]; if(is.na(p_val)||p_val>=0.05) next
        lbl <- sig_label(p_val)
        xi <- which(geno_levels==rownames(pw_mat)[r])
        xj <- which(geno_levels==colnames(pw_mat)[cc])
        if(!length(xi)||!length(xj)) next
        brk_y <- ylim_l[2]*0.82+drawn*bracket_step
        par(xpd=NA); draw_bracket(xj,xi,brk_y,lbl,tip=(ylim_l[2]-ylim_l[1])*0.012)
        par(xpd=FALSE); drawn <- drawn+1
      }
      for (r in seq_len(nrow(pw_bh_mat))) for (cc in seq_len(ncol(pw_bh_mat))) {
        bh_p <- pw_bh_mat[r,cc]; hl_p <- pw_mat[r,cc]
        if(is.na(bh_p)||bh_p>=0.05||(!is.na(hl_p)&&hl_p<0.05)) next
        lbl <- sig_label(bh_p)
        xi <- which(geno_levels==rownames(pw_bh_mat)[r])
        xj <- which(geno_levels==colnames(pw_bh_mat)[cc])
        if(!length(xi)||!length(xj)) next
        brk_y <- ylim_l[2]*0.82+drawn*bracket_step
        par(xpd=NA)
        draw_bracket(xj,xi,brk_y,lbl,tip=(ylim_l[2]-ylim_l[1])*0.012,col="#2980B9")
        par(xpd=FALSE); drawn <- drawn+1
      }
      mtext(paste0("KW p=",formatC(kw$p.value,digits=3,format="g"),
                   "  |  Black=Holm, Blue=BH-only (Mann-Whitney U pairwise)"),
            side=1,line=5,cex=0.95,col="grey40",font=2)
      legend("topright",legend=geno_levels,
             fill=gc$fill[geno_levels],border="black",bty="n",cex=1.2,text.font=2)
    }
  }
}
datasets_lat <- list("All data"=latency, "IQR outliers removed"=latency_clean)
for (label in names(datasets_lat)) run_latency_plot(datasets_lat[[label]], label)

# =============================================================================
# ANALYSIS 3: Threshold during sleep, across genotypes
# Omnibus (KW or ANOVA) FIRST; pairwise only if omnibus significant
# Holm PRIMARY | BH SENSITIVITY
# =============================================================================
run_sleeping_plot <- function(d, label) {
  d$genotype  <- order_genos(droplevels(d$genotype))
  geno_levels <- levels(d$genotype)
  n_genos     <- length(geno_levels)
  gc          <- make_geno_colours(geno_levels)
  
  cat(paste0("\n", strrep("=",70), "\n"))
  cat(paste0("ANALYSIS 3 — Threshold (Sleeping Only)  [", label, "]\n"))
  cat("Omnibus test gates pairwise comparisons\n")
  cat(paste0(strrep("=",70), "\n"))
  cat("N:\n"); print(table(d$genotype))
  for (g in geno_levels) {
    v <- d$volume_percent_threshold[d$genotype==g]
    cat(sprintf("  %s: mean=%.2f, SD=%.2f\n",g,mean(v,na.rm=TRUE),sd(v,na.rm=TRUE)))
  }
  
  chk3 <- check_assumptions(d$volume_percent_threshold, d$genotype,
                            paste0("Analysis 3 ",label))
  
  ylim_s <- c(0,max(d$volume_percent_threshold,na.rm=TRUE)*1.6)
  par(mar=c(7,6,5,2),mfrow=c(1,1),
      cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
      font=2,font.main=2,font.lab=2,font.axis=2)
  xs <- seq_len(n_genos)
  plot(0,0,type="n",xlim=c(0.5,n_genos+0.5),ylim=ylim_s,
       xlab="",ylab="Volume % Threshold",
       main=paste0("Threshold — Sleeping Flies  [",label,"]"),
       xaxt="n",las=1,bty="n")
  axis(1,at=xs,labels=geno_levels,tick=FALSE,cex.axis=1.2,las=2,font=2)
  for (i in xs)
    draw_manual_box(i,d$volume_percent_threshold[d$genotype==geno_levels[i]],
                    gc$fill[geno_levels[i]],gc$border[geno_levels[i]])
  
  if (chk3$parametric) {
    aov_m <- aov(volume_percent_threshold~genotype, data=d)
    omnibus_p <- summary(aov_m)[[1]][["Pr(>F)"]][1]
    cat("\nOne-way ANOVA:\n"); print(summary(aov_m))
    omnibus_str_base <- paste0("ANOVA p=",formatC(omnibus_p,digits=3,format="g"))
  } else {
    kw3 <- kruskal.test(volume_percent_threshold~genotype, data=d)
    omnibus_p <- kw3$p.value
    cat("\nKruskal-Wallis:\n"); print(kw3)
    omnibus_str_base <- paste0("KW p=",formatC(omnibus_p,digits=3,format="g"))
  }
  
  if (omnibus_p >= 0.05) {
    cat("\n  Omnibus not significant — pairwise tests NOT run (gate).\n")
    mtext(paste0(omnibus_str_base," (ns) — no pairwise tests"),
          side=1,line=5,cex=0.95,col="grey40",font=2)
  } else {
    cat("\n  Omnibus significant — running pairwise tests\n")
    if (chk3$parametric) {
      pw_holm <- pairwise.t.test(d$volume_percent_threshold,d$genotype,
                                 p.adjust.method="holm",pool.sd=FALSE)
      pw_bh   <- pairwise.t.test(d$volume_percent_threshold,d$genotype,
                                 p.adjust.method="BH",  pool.sd=FALSE)
      cat("\nPairwise Welch — Holm [PRIMARY]:\n"); print(pw_holm)
      cat("\nPairwise Welch — BH [SENSITIVITY]:\n"); print(pw_bh)
    } else {
      pw_holm <- pairwise.wilcox.test(d$volume_percent_threshold,d$genotype,
                                      p.adjust.method="holm",exact=FALSE)
      pw_bh   <- pairwise.wilcox.test(d$volume_percent_threshold,d$genotype,
                                      p.adjust.method="BH",  exact=FALSE)
      cat("\nPairwise Mann-Whitney U — Holm [PRIMARY]:\n"); print(pw_holm)
      cat("\nPairwise Mann-Whitney U — BH [SENSITIVITY]:\n"); print(pw_bh)
    }
    pw_mat_h <- pw_holm$p.value; pw_mat_b <- pw_bh$p.value
    bracket_step <- (ylim_s[2]-ylim_s[1])*0.07; drawn <- 0
    for (r in seq_len(nrow(pw_mat_h))) for (cc in seq_len(ncol(pw_mat_h))) {
      p_val <- pw_mat_h[r,cc]; if(is.na(p_val)||p_val>=0.05) next
      lbl <- sig_label(p_val)
      xi <- which(geno_levels==rownames(pw_mat_h)[r])
      xj <- which(geno_levels==colnames(pw_mat_h)[cc])
      if(!length(xi)||!length(xj)) next
      brk_y <- ylim_s[2]*0.80+drawn*bracket_step
      par(xpd=NA); draw_bracket(xj,xi,brk_y,lbl,tip=(ylim_s[2]-ylim_s[1])*0.012)
      par(xpd=FALSE); drawn <- drawn+1
    }
    for (r in seq_len(nrow(pw_mat_b))) for (cc in seq_len(ncol(pw_mat_b))) {
      bh_p <- pw_mat_b[r,cc]; hl_p <- pw_mat_h[r,cc]
      if(is.na(bh_p)||bh_p>=0.05||(!is.na(hl_p)&&hl_p<0.05)) next
      lbl <- sig_label(bh_p)
      xi <- which(geno_levels==rownames(pw_mat_b)[r])
      xj <- which(geno_levels==colnames(pw_mat_b)[cc])
      if(!length(xi)||!length(xj)) next
      brk_y <- ylim_s[2]*0.80+drawn*bracket_step
      par(xpd=NA)
      draw_bracket(xj,xi,brk_y,lbl,tip=(ylim_s[2]-ylim_s[1])*0.012,col="#2980B9")
      par(xpd=FALSE); drawn <- drawn+1
    }
    mtext(paste0(omnibus_str_base,"  |  Black=Holm, Blue=BH-only"),
          side=1,line=5,cex=0.95,col="grey40",font=2)
  }
  legend("topright",legend=geno_levels,
         fill=gc$fill[geno_levels],border="black",bty="n",cex=1.2,text.font=2)
}
sleeping_only       <- df[df$sleep_cond=="sleeping",]
sleeping_only_clean <- remove_iqr(sleeping_only,"volume_percent_threshold","genotype")
for (d_name in c("All data","IQR outliers removed")) {
  d_use <- if(d_name=="All data") sleeping_only else sleeping_only_clean
  run_sleeping_plot(d_use, d_name)
}

# =============================================================================
# ANALYSIS 4: Correlation — Spearman (unchanged, always appropriate)
# =============================================================================
cor_df       <- latency[!is.na(latency$volume_percent_threshold),]
cor_df_clean <- remove_iqr(remove_iqr(cor_df,"sleep_latency","genotype"),
                           "volume_percent_threshold","genotype")
for (label in c("All data","IQR outliers removed")) {
  d           <- if(label=="All data") cor_df else cor_df_clean
  d$genotype  <- order_genos(droplevels(d$genotype))
  geno_levels <- levels(d$genotype); n_genos <- length(geno_levels)
  gc          <- make_geno_colours(geno_levels)
  cat(paste0("\n", strrep("=",70), "\n"))
  cat(paste0("ANALYSIS 4 — Spearman Correlation  [", label, "]\n"))
  cat(paste0(strrep("=",70), "\n"))
  ncols <- min(n_genos,3); nrows <- ceiling(n_genos/ncols)
  par(mfrow=c(nrows,ncols),mar=c(5.5,5.5,4,2),oma=c(0,0,4,0),
      cex=1.3,cex.main=1.5,cex.lab=1.3,cex.axis=1.2,
      font=2,font.main=2,font.lab=2,font.axis=2)
  for (g in geno_levels) {
    sub <- d[d$genotype==g,]
    x <- sub$sleep_latency; y <- sub$volume_percent_threshold
    sct <- cor.test(x,y,method="spearman")
    pct <- cor.test(x,y,method="pearson")
    cat(paste0("\n--- ",g," ---\n"))
    cat("Spearman [PRIMARY]:\n"); print(sct)
    cat("Pearson (sensitivity):\n"); print(pct)
    sp_p <- sct$p.value; r_val <- round(sct$estimate,3)
    p_str <- if(sp_p<0.001)"p<0.001" else paste0("p=",round(sp_p,3))
    sp_lbl<- if(sp_p<0.001)"***" else if(sp_p<0.01)"**" else if(sp_p<0.05)"*" else "ns"
    plot(x,y,pch=16,col=adjustcolor(gc$border[g],alpha.f=0.7),cex=1.1,
         xlab="Sleep Latency (s)",ylab="Volume % Threshold",main=g,las=1,bty="n")
    if(length(x)>2) abline(lm(y~x),col="black",lwd=2,lty=2)
    legend("topright",legend=c(paste0("\u03c1=",r_val),paste0(p_str," (",sp_lbl,")")),
           bty="n",cex=1.1,text.font=2)
  }
  mtext(paste0("Latency vs Threshold  [",label,"] — Spearman"),
        outer=TRUE,cex=1.2,font=2,line=1)
}

# =============================================================================
# ANALYSES 5 + 10: LATENCY — HK/canton_s and inc/w1118
# Treated as ONE FAMILY of 2 tests — Holm across both
# =============================================================================
cat(paste0("\n", strrep("=",70), "\n"))
cat("ANALYSES 5 + 10 — Sleep Latency: HK/Canton-S and Inc/w1118\n")
cat("CORRECTION: Holm across BOTH tests as one family (PRIMARY)\n")
cat("            BH across both (SENSITIVITY)\n")
cat(paste0(strrep("=",70), "\n"))

lat_hk_cs <- latency[trimws(tolower(as.character(latency$genotype))) %in%
                       c("hk","canton_s"),]
lat_hk_cs$genotype <- order_genos(droplevels(factor(trimws(
  as.character(lat_hk_cs$genotype)))))

lat_inc_w <- latency[trimws(tolower(as.character(latency$genotype))) %in%
                       c("inc","w1118"),]
lat_inc_w$genotype <- order_genos(droplevels(factor(trimws(
  as.character(lat_inc_w$genotype)))))

# Collect raw p-values for both comparisons
lat_raw_ps   <- c(hk_cs=NA_real_, inc_w118=NA_real_)
lat_use_wc   <- c(hk_cs=NA, inc_w118=NA)
lat_res_list <- list()

if (nlevels(lat_hk_cs$genotype)>=2) {
  gv <- levels(lat_hk_cs$genotype)
  cat(paste0("\n--- Analysis 5: ",gv[1]," vs ",gv[2]," ---\n"))
  res5 <- raw_two_group_test(
    lat_hk_cs$sleep_latency[lat_hk_cs$genotype==gv[1]],
    lat_hk_cs$sleep_latency[lat_hk_cs$genotype==gv[2]],
    gv[1], gv[2])
  lat_raw_ps["hk_cs"]  <- res5$primary_p_raw
  lat_use_wc["hk_cs"]  <- res5$use_wilcoxon
  lat_res_list[["hk_cs"]] <- list(res=res5, gv=gv,
                                  v1=lat_hk_cs$sleep_latency[lat_hk_cs$genotype==gv[1]],
                                  v2=lat_hk_cs$sleep_latency[lat_hk_cs$genotype==gv[2]])
}

if (nlevels(lat_inc_w$genotype)>=2) {
  gv <- levels(lat_inc_w$genotype)
  cat(paste0("\n--- Analysis 10: ",gv[1]," vs ",gv[2]," ---\n"))
  res10 <- raw_two_group_test(
    lat_inc_w$sleep_latency[lat_inc_w$genotype==gv[1]],
    lat_inc_w$sleep_latency[lat_inc_w$genotype==gv[2]],
    gv[1], gv[2])
  lat_raw_ps["inc_w118"] <- res10$primary_p_raw
  lat_use_wc["inc_w118"] <- res10$use_wilcoxon
  lat_res_list[["inc_w118"]] <- list(res=res10, gv=gv,
                                     v1=lat_inc_w$sleep_latency[lat_inc_w$genotype==gv[1]],
                                     v2=lat_inc_w$sleep_latency[lat_inc_w$genotype==gv[2]])
}

# Apply Holm and BH across the family of 2
valid_lat <- !is.na(lat_raw_ps)
lat_p_holm <- lat_raw_ps; lat_p_bh <- lat_raw_ps
lat_p_holm[valid_lat] <- p.adjust(lat_raw_ps[valid_lat], method="holm")
lat_p_bh[valid_lat]   <- p.adjust(lat_raw_ps[valid_lat], method="BH")

cat(paste0("\n", strrep("-",60), "\n"))
cat("LATENCY FAMILY CORRECTION SUMMARY\n")
cat(sprintf("  %-12s  %-8s  %-10s  %-10s  %-10s\n",
            "Comparison","Test","Raw p","Holm p","BH p"))
for (nm in names(lat_raw_ps)) {
  if (is.na(lat_raw_ps[nm])) next
  cat(sprintf("  %-12s  %-8s  %-10.4f  %-10.4f  %-10.4f\n",
              nm,
              if(isTRUE(lat_use_wc[nm]))"MWU" else "Welch",
              lat_raw_ps[nm], lat_p_holm[nm], lat_p_bh[nm]))
}

# Plot Analysis 5
if (!is.null(lat_res_list[["hk_cs"]])) {
  r  <- lat_res_list[["hk_cs"]]
  plot_two_group(r$v1, r$v2, r$gv[1], r$gv[2],
                 lat_p_holm["hk_cs"], lat_p_bh["hk_cs"],
                 isTRUE(lat_use_wc["hk_cs"]),
                 ylab="Sleep Latency (s)",
                 main_title=paste0("Sleep Latency: ",r$gv[1]," vs ",r$gv[2],
                                   "\n[Holm-corrected within latency family]"))
}

# Plot Analysis 10
if (!is.null(lat_res_list[["inc_w118"]])) {
  r  <- lat_res_list[["inc_w118"]]
  plot_two_group(r$v1, r$v2, r$gv[1], r$gv[2],
                 lat_p_holm["inc_w118"], lat_p_bh["inc_w118"],
                 isTRUE(lat_use_wc["inc_w118"]),
                 ylab="Sleep Latency (s)",
                 main_title=paste0("Sleep Latency: ",r$gv[1]," vs ",r$gv[2],
                                   "\n[Holm-corrected within latency family]"))
}

# =============================================================================
# ANALYSIS 6: Waking proportion HK vs Canton-S — chi-square + Fisher's exact
# Single pre-specified test; no correction needed
# =============================================================================
cat(paste0("\n", strrep("=",70), "\n"))
cat("ANALYSIS 6 — Proportion waking: HK vs Canton-S\n")
cat("Single pre-specified test — no multiple comparison correction needed\n")
cat(paste0(strrep("=",70), "\n"))
react_df <- sleeping_raw[
  trimws(tolower(as.character(sleeping_raw$genotype))) %in% c("hk","canton_s") &
    !is.na(sleeping_raw$reaction), ]
react_df$genotype <- order_genos(droplevels(factor(trimws(
  as.character(react_df$genotype)))))
if (nlevels(react_df$genotype)>=2) {
  genos_r <- levels(react_df$genotype)
  ct <- table(Genotype=react_df$genotype,
              Woke=factor(react_df$reaction,levels=c(0,1),
                          labels=c("No wake","Woke")))
  cat("\nContingency table:\n"); print(ct)
  cat("\nExpected counts:\n"); print(chisq.test(ct)$expected)
  cat("\nRow proportions:\n"); print(round(prop.table(ct,margin=1)*100,1))
  cs <- chisq.test(ct,correct=FALSE); fe <- fisher.test(ct)
  cat("\nChi-square [PRIMARY]:\n"); print(cs)
  cat("\nFisher's Exact [SENSITIVITY — preferred when expected counts < 5]:\n")
  print(fe)
  props_woke <- prop.table(ct,margin=1)[,"Woke"]*100
  par(mfrow=c(1,1),mar=c(6,6,5,4),
      cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
      font=2,font.main=2,font.lab=2,font.axis=2)
  bp6 <- barplot(props_woke,col=c("#C0392B","#4A7FC1")[seq_along(genos_r)],
                 border="black",ylim=c(0,100),names.arg=genos_r,
                 ylab="Flies that woke up (%)",
                 main="Proportion Waking: HK vs Canton-S",las=1,bty="n",lwd=2)
  text(bp6,props_woke-4,paste0(round(props_woke,1),"%"),col="white",font=2,cex=1.2)
  cs_str <- if(cs$p.value<0.001)"p<0.001" else paste0("p=",round(cs$p.value,3))
  fe_str <- if(fe$p.value<0.001)"p<0.001" else paste0("p=",round(fe$p.value,3))
  cs_sig  <- if(cs$p.value<0.001)"***" else if(cs$p.value<0.01)"**" else
    if(cs$p.value<0.05)"*" else "ns"
  mtext(paste0("Chi-sq ",cs_str," | Fisher ",fe_str," (",cs_sig,")"),
        side=1,line=4,cex=1.0,col="grey40",font=2)
  if(cs$p.value<0.05) {
    par(xpd=NA); draw_bracket(bp6[1],bp6[2],y=92,label=cs_sig,tip=1.5)
    par(xpd=FALSE)
  }
}

# =============================================================================
# ANALYSES 7 + 11: THRESHOLD DURING SLEEP — HK/canton_s and inc/w1118
# Treated as ONE FAMILY of 2 tests — Holm across both
# =============================================================================
cat(paste0("\n", strrep("=",70), "\n"))
cat("ANALYSES 7 + 11 — Arousal Threshold (Sleeping): HK/Canton-S and Inc/w1118\n")
cat("CORRECTION: Holm across BOTH tests as one family (PRIMARY)\n")
cat("            BH across both (SENSITIVITY)\n")
cat(paste0(strrep("=",70), "\n"))

sleeping_only <- df[df$sleep_cond=="sleeping",]

thresh_hk_cs <- sleeping_only[
  trimws(tolower(as.character(sleeping_only$genotype))) %in% c("hk","canton_s"),]
thresh_hk_cs$genotype <- order_genos(droplevels(factor(trimws(
  as.character(thresh_hk_cs$genotype)))))

thresh_inc_w <- sleeping_only[
  trimws(tolower(as.character(sleeping_only$genotype))) %in% c("inc","w1118"),]
thresh_inc_w$genotype <- order_genos(droplevels(factor(trimws(
  as.character(thresh_inc_w$genotype)))))

thr_raw_ps   <- c(hk_cs=NA_real_, inc_w118=NA_real_)
thr_use_wc   <- c(hk_cs=NA, inc_w118=NA)
thr_res_list <- list()

if (nlevels(thresh_hk_cs$genotype)>=2) {
  gv7 <- levels(thresh_hk_cs$genotype)
  cat(paste0("\n--- Analysis 7: ",gv7[1]," vs ",gv7[2]," (sleeping) ---\n"))
  res7 <- raw_two_group_test(
    thresh_hk_cs$volume_percent_threshold[thresh_hk_cs$genotype==gv7[1]],
    thresh_hk_cs$volume_percent_threshold[thresh_hk_cs$genotype==gv7[2]],
    gv7[1], gv7[2])
  thr_raw_ps["hk_cs"]  <- res7$primary_p_raw
  thr_use_wc["hk_cs"]  <- res7$use_wilcoxon
  thr_res_list[["hk_cs"]] <- list(res=res7, gv=gv7,
                                  v1=thresh_hk_cs$volume_percent_threshold[thresh_hk_cs$genotype==gv7[1]],
                                  v2=thresh_hk_cs$volume_percent_threshold[thresh_hk_cs$genotype==gv7[2]])
}

if (nlevels(thresh_inc_w$genotype)>=2) {
  gv11 <- levels(thresh_inc_w$genotype)
  cat(paste0("\n--- Analysis 11: ",gv11[1]," vs ",gv11[2]," (sleeping) ---\n"))
  res11 <- raw_two_group_test(
    thresh_inc_w$volume_percent_threshold[thresh_inc_w$genotype==gv11[1]],
    thresh_inc_w$volume_percent_threshold[thresh_inc_w$genotype==gv11[2]],
    gv11[1], gv11[2])
  thr_raw_ps["inc_w118"] <- res11$primary_p_raw
  thr_use_wc["inc_w118"] <- res11$use_wilcoxon
  thr_res_list[["inc_w118"]] <- list(res=res11, gv=gv11,
                                     v1=thresh_inc_w$volume_percent_threshold[thresh_inc_w$genotype==gv11[1]],
                                     v2=thresh_inc_w$volume_percent_threshold[thresh_inc_w$genotype==gv11[2]])
}

# Holm and BH across the family of 2
valid_thr <- !is.na(thr_raw_ps)
thr_p_holm <- thr_raw_ps; thr_p_bh <- thr_raw_ps
thr_p_holm[valid_thr] <- p.adjust(thr_raw_ps[valid_thr], method="holm")
thr_p_bh[valid_thr]   <- p.adjust(thr_raw_ps[valid_thr], method="BH")

cat(paste0("\n", strrep("-",60), "\n"))
cat("THRESHOLD (SLEEPING) FAMILY CORRECTION SUMMARY\n")
cat(sprintf("  %-12s  %-8s  %-10s  %-10s  %-10s\n",
            "Comparison","Test","Raw p","Holm p","BH p"))
for (nm in names(thr_raw_ps)) {
  if (is.na(thr_raw_ps[nm])) next
  cat(sprintf("  %-12s  %-8s  %-10.4f  %-10.4f  %-10.4f\n",
              nm,
              if(isTRUE(thr_use_wc[nm]))"MWU" else "Welch",
              thr_raw_ps[nm], thr_p_holm[nm], thr_p_bh[nm]))
}

# Cohen's d
for (nm in names(thr_res_list)) {
  r   <- thr_res_list[[nm]]
  cd  <- cohens_d_ci(r$v1, r$v2)
  d_mag <- ifelse(abs(cd["d"])>=0.8,"large",ifelse(abs(cd["d"])>=0.5,"medium","small"))
  cat(sprintf("  Cohen's d [%s]: %.3f [%.3f, %.3f] (%s)\n",
              nm, cd["d"], cd["lower"], cd["upper"], d_mag))
}

# Plot Analysis 7
if (!is.null(thr_res_list[["hk_cs"]])) {
  r <- thr_res_list[["hk_cs"]]
  plot_two_group(r$v1, r$v2, r$gv[1], r$gv[2],
                 thr_p_holm["hk_cs"], thr_p_bh["hk_cs"],
                 isTRUE(thr_use_wc["hk_cs"]),
                 ylab="Volume % Threshold",
                 main_title=paste0("Arousal Threshold (Sleeping): ",
                                   r$gv[1]," vs ",r$gv[2],
                                   "\n[Holm-corrected within threshold family]"),
                 ylim_override=c(0,100))
}

# Plot Analysis 11
if (!is.null(thr_res_list[["inc_w118"]])) {
  r <- thr_res_list[["inc_w118"]]
  plot_two_group(r$v1, r$v2, r$gv[1], r$gv[2],
                 thr_p_holm["inc_w118"], thr_p_bh["inc_w118"],
                 isTRUE(thr_use_wc["inc_w118"]),
                 ylab="Volume % Threshold",
                 main_title=paste0("Arousal Threshold (Sleeping): ",
                                   r$gv[1]," vs ",r$gv[2],
                                   "\n[Holm-corrected within threshold family]"),
                 ylim_override=c(0,100))
}

# =============================================================================
# ANALYSIS 8: HK vs Canton-S — censored latency
# Sensitivity check for Analysis 5 — reported separately, Mann-Whitney U always
# =============================================================================
cat(paste0("\n", strrep("=",70), "\n"))
cat("ANALYSIS 8 — Sleep Latency incl. censored: HK vs Canton-S\n")
cat("Sensitivity check for Analysis 5 (censored data — Mann-Whitney U mandatory)\n")
cat(paste0(strrep("=",70), "\n"))
lat_hk_cs_cap <- sleeping_raw[
  trimws(tolower(as.character(sleeping_raw$genotype))) %in% c("hk","canton_s") &
    !is.na(sleeping_raw$sleep_latency) & sleeping_raw$sleep_latency<5000, ]
lat_hk_cs_cap$genotype <- order_genos(droplevels(factor(trimws(
  as.character(lat_hk_cs_cap$genotype)))))
if (nlevels(lat_hk_cs_cap$genotype)>=2) {
  gv8 <- levels(lat_hk_cs_cap$genotype)
  g1l8 <- lat_hk_cs_cap$sleep_latency[lat_hk_cs_cap$genotype==gv8[1]]
  g2l8 <- lat_hk_cs_cap$sleep_latency[lat_hk_cs_cap$genotype==gv8[2]]
  n_c1 <- sum(g1l8==4000,na.rm=TRUE); n_c2 <- sum(g2l8==4000,na.rm=TRUE)
  cat(sprintf("N: %s=%d (%d censored) | %s=%d (%d censored)\n",
              gv8[1],length(g1l8),n_c1,gv8[2],length(g2l8),n_c2))
  wc8 <- wilcox.test(g1l8,g2l8,exact=FALSE); tt8 <- t.test(g1l8,g2l8)
  cat("\nMann-Whitney U [PRIMARY — censored data]:\n"); print(wc8)
  cat("\nWelch t (sensitivity only):\n"); print(tt8)
  ylim_8 <- c(0,5000*1.30); gc8 <- make_geno_colours(gv8)
  par(mfrow=c(1,1),mar=c(6,6,5,2),
      cex=1.4,cex.main=1.6,cex.lab=1.4,cex.axis=1.3,
      font=2,font.main=2,font.lab=2,font.axis=2)
  plot(0,0,type="n",xlim=c(0.5,2.5),ylim=ylim_8,
       xlab="Genotype",ylab="Sleep Latency (s)",
       main="Sleep Latency (incl. censored): HK vs Canton-S\n[Mann-Whitney U — sensitivity check for Analysis 5]",
       xaxt="n",las=1,bty="n")
  axis(1,at=1:2,labels=gv8,tick=FALSE,cex.axis=1.3,font=2)
  draw_manual_box(1,g1l8,gc8$fill[gv8[1]],gc8$border[gv8[1]])
  draw_manual_box(2,g2l8,gc8$fill[gv8[2]],gc8$border[gv8[2]])
  two_group_bracket(ylim_8,1,2,wc8$p.value)
  p_footer(wc8$p.value,tt8$p.value,"Mann-Whitney U [PRIMARY]","Welch t")
  mtext(paste0("n=",length(g1l8)," (",n_c1," censored)"),
        side=1,at=1,line=2,cex=0.95,col="grey35",font=2)
  mtext(paste0("n=",length(g2l8)," (",n_c2," censored)"),
        side=1,at=2,line=2,cex=0.95,col="grey35",font=2)
  abline(h=4000,lty=2,col="grey60",lwd=1.5)
  legend("topright",legend=gv8,fill=gc8$fill[gv8],border="black",
         bty="n",cex=1.2,text.font=2)
}

# =============================================================================
# ANALYSIS 9: Threshold >1%, rested vs sleeping
# Same as Analysis 1 but filtered to >1% — same correction structure (family of 4)
# =============================================================================
cat(paste0("\n", strrep("=",70), "\n"))
cat("ANALYSIS 9 — Threshold >1%: Awake vs Sleeping\n")
cat(paste0(strrep("=",70), "\n"))
df_above1 <- df[!is.na(df$volume_percent_threshold) & df$volume_percent_threshold>1,]
df_above1_clean <- df_clean[!is.na(df_clean$volume_percent_threshold) &
                              df_clean$volume_percent_threshold>1,]
for (label in c("All data","IQR outliers removed")) {
  d_use <- if(label=="All data") df_above1 else df_above1_clean
  run_analysis1(d_use, paste0(label," [>1% threshold]"))
}

# =============================================================================
# DESCRIPTIVE STATISTICS
# =============================================================================
cat(paste0("\n", strrep("=",70), "\n"))
cat("DESCRIPTIVE STATISTICS\n")
cat(paste0(strrep("=",70), "\n"))
desc_stats <- do.call(rbind, lapply(levels(df$genotype), function(g) {
  do.call(rbind, lapply(levels(df$sleep_cond), function(cond) {
    vals <- df$volume_percent_threshold[df$genotype==g & df$sleep_cond==cond]
    vals <- vals[!is.na(vals)]
    if(length(vals)==0) return(NULL)
    data.frame(genotype=g,sleep_cond=cond,n=length(vals),
               mean=round(mean(vals),2),sd=round(sd(vals),2),
               stringsAsFactors=FALSE)
  }))
}))
print(desc_stats, row.names=FALSE)
# =============================================================================
# ANALYSIS 12: Two-way factorial — Genotype x Sleep Condition
#
# Omnibus interaction test (ART-ANOVA or two-way ANOVA) reported for context.
#
# CORRECTION STRUCTURE — separate family of 4 for each mutant/control pair:
#
#   Analysis 12A (HK vs canton_s):
#     (1) HK vs canton_s        | awake
#     (2) HK vs canton_s        | sleeping
#     (3) awake vs sleeping     | HK
#     (4) awake vs sleeping     | canton_s
#     → Holm across all 4 (PRIMARY) | BH across all 4 (SENSITIVITY)
#
#   Analysis 12B (inc vs w1118) — identical structure, corrected independently:
#     (1) inc vs w1118          | awake
#     (2) inc vs w1118          | sleeping
#     (3) awake vs sleeping     | inc
#     (4) awake vs sleeping     | w1118
#     → Holm across all 4 (PRIMARY) | BH across all 4 (SENSITIVITY)
#
# Black brackets = Holm (PRIMARY) | Blue = BH only (SENSITIVITY)
# =============================================================================
cat(paste0("\n", strrep("=",70), "\n"))
cat("ANALYSIS 12 — Two-way Factorial: Genotype x Sleep Condition\n")
cat("Family of 4 per pair (genotype contrast x2 + condition contrast x2)\n")
cat("12A (HK/canton_s) and 12B (inc/w1118) corrected as SEPARATE families\n")
cat(paste0(strrep("=",70), "\n"))

# ---------------------------------------------------------------------------
# PHASE 1: Collect all 4 raw p-values for one mutant/control pair.
# Returns everything needed for plotting once corrected p-values are known.
# ---------------------------------------------------------------------------
collect_analysis12 <- function(geno_mutant, geno_control, d,
                               label = paste0(geno_mutant, " vs ", geno_control)) {
  suffix <- if (geno_mutant == "hk") "A" else "B"
  cat(paste0("\n", strrep("-",60), "\n"))
  cat(paste0("ANALYSIS 12", suffix, " — ", label, " [raw p collection]\n"))
  cat(paste0(strrep("-",60), "\n"))
  
  sub <- d[trimws(tolower(as.character(d$genotype))) %in%
             c(tolower(geno_mutant), tolower(geno_control)) &
             !is.na(d$volume_percent_threshold) & d$volume_percent_threshold > 1, ]
  if (nrow(sub) == 0) { cat("  No data.\n"); return(NULL) }
  
  sub$genotype   <- factor(trimws(as.character(sub$genotype)),
                           levels = c(geno_control, geno_mutant))
  sub$sleep_cond <- factor(sub$sleep_cond, levels = c("rested", "sleeping"))
  sub <- sub[!is.na(sub$sleep_cond), ]
  
  cat("\nSample sizes:\n"); print(table(sub$genotype, sub$sleep_cond))
  cat("\nCell means:\n")
  print(round(tapply(sub$volume_percent_threshold,
                     list(sub$genotype, sub$sleep_cond), mean, na.rm = TRUE), 2))
  
  # Assumption checks for omnibus
  cat(paste0("\n", strrep("-",60), "\n"))
  cat("ASSUMPTION CHECKS (omnibus)\n")
  cat(paste0(strrep("-",60), "\n"))
  model_lm <- lm(volume_percent_threshold ~ genotype * sleep_cond, data = sub)
  resids   <- residuals(model_lm)
  sw_resid <- shapiro.test(resids); sw_ok <- sw_resid$p.value >= 0.05
  cat(sprintf("\n  (1) SW on residuals: W=%.4f, p=%.4f  -> %s\n",
              sw_resid$statistic, sw_resid$p.value,
              if (sw_ok) "NORMAL" else "*** NON-NORMAL ***"))
  sub$cell <- interaction(sub$genotype, sub$sleep_cond)
  lev      <- car::leveneTest(volume_percent_threshold ~ cell, data = sub, center = median)
  lev_p    <- lev[["Pr(>F)"]][1]; lev_ok <- lev_p >= 0.05
  cat(sprintf("  (2) Levene's: F(%d,%d)=%.4f, p=%.4f  -> %s\n",
              lev[["Df"]][1], lev[["Df"]][2], lev[["F value"]][1], lev_p,
              if (lev_ok) "HOMOGENEOUS" else "*** HETEROGENEOUS ***"))
  
  use_art <- !(sw_ok && lev_ok)
  cat(sprintf("  -> Omnibus: %s\n\n",
              if (use_art) "ART-ANOVA (assumptions violated)"
              else         "Standard two-way ANOVA (assumptions met)"))
  
  int_p <- NA
  if (!use_art) {
    aov_t2  <- car::Anova(model_lm, type = "II")
    cat("Type II ANOVA:\n"); print(aov_t2)
    int_row <- rownames(aov_t2)[grepl("genotype:sleep_cond", rownames(aov_t2), fixed = TRUE)]
    int_p   <- if (length(int_row) == 1) aov_t2[int_row, "Pr(>F)"] else NA
  } else {
    cat("ART-ANOVA:\n")
    art_model   <- ARTool::art(volume_percent_threshold ~ genotype * sleep_cond, data = sub)
    art_aov     <- anova(art_model); print(art_aov)
    int_row_art <- rownames(art_aov)[grepl("genotype:sleep_cond",
                                           rownames(art_aov), fixed = TRUE)]
    if (length(int_row_art) == 1) int_p <- art_aov[int_row_art, "Pr(>F)"]
  }
  
  int_sig_str <- if (is.na(int_p))    "N/A"
  else if (int_p < 0.001)             "< 0.001 ***"
  else if (int_p < 0.01)             sprintf("%.4f **",  int_p)
  else if (int_p < 0.05)             sprintf("%.4f *",   int_p)
  else                               sprintf("%.4f (ns)", int_p)
  cat(sprintf("\n  -> Interaction p = %s\n", int_sig_str))
  
  # -------------------------------------------------------------------
  # Collect raw p-values for all 4 tests:
  #   geno_rested   — mutant vs control within awake
  #   geno_sleeping — mutant vs control within sleeping
  #   cond_mutant   — awake vs sleeping within mutant
  #   cond_control  — awake vs sleeping within control
  # -------------------------------------------------------------------
  raw_ps  <- setNames(rep(NA_real_, 4),
                      c("geno_rested", "geno_sleeping",
                        "cond_mutant", "cond_control"))
  use_wc  <- setNames(rep(NA, 4), names(raw_ps))
  v_store <- list()
  
  cat(paste0("\n", strrep("-",60), "\n"))
  cat("RAW PAIRWISE TESTS — genotype contrasts within condition\n")
  cat(paste0(strrep("-",60), "\n"))
  
  for (cond in c("rested", "sleeping")) {
    key    <- paste0("geno_", cond)
    v_ctrl <- sub$volume_percent_threshold[sub$genotype == geno_control &
                                             sub$sleep_cond == cond]
    v_mut  <- sub$volume_percent_threshold[sub$genotype == geno_mutant  &
                                             sub$sleep_cond == cond]
    v_ctrl <- v_ctrl[!is.na(v_ctrl)]; v_mut <- v_mut[!is.na(v_mut)]
    v_store[[key]] <- list(ctrl = v_ctrl, mut = v_mut)
    if (length(v_ctrl) < 2 || length(v_mut) < 2) {
      cat(sprintf("  %s vs %s [%s]: SKIPPED\n", geno_control, geno_mutant, cond)); next
    }
    cat(sprintf("\n  --- %s vs %s | %s ---\n", geno_control, geno_mutant, cond))
    res          <- raw_two_group_test(v_ctrl, v_mut, geno_control, geno_mutant)
    raw_ps[key]  <- res$primary_p_raw
    use_wc[key]  <- res$use_wilcoxon
  }
  
  cat(paste0("\n", strrep("-",60), "\n"))
  cat("RAW PAIRWISE TESTS — condition contrasts within genotype\n")
  cat(paste0(strrep("-",60), "\n"))
  
  for (geno in c(geno_mutant, geno_control)) {
    key <- paste0("cond_", if (geno == geno_mutant) "mutant" else "control")
    v_r <- sub$volume_percent_threshold[sub$genotype == geno &
                                          sub$sleep_cond == "rested"]
    v_s <- sub$volume_percent_threshold[sub$genotype == geno &
                                          sub$sleep_cond == "sleeping"]
    v_r <- v_r[!is.na(v_r)]; v_s <- v_s[!is.na(v_s)]
    v_store[[key]] <- list(rested = v_r, sleeping = v_s)
    if (length(v_r) < 2 || length(v_s) < 2) {
      cat(sprintf("  awake vs sleeping [%s]: SKIPPED\n", geno)); next
    }
    cat(sprintf("\n  --- awake vs sleeping | %s ---\n", geno))
    res         <- raw_two_group_test(v_r, v_s, "rested", "sleeping")
    raw_ps[key] <- res$primary_p_raw
    use_wc[key] <- res$use_wilcoxon
  }
  
  list(
    label        = label,
    suffix       = suffix,
    sub          = sub,
    geno_mutant  = geno_mutant,
    geno_control = geno_control,
    raw_ps       = raw_ps,
    use_wc       = use_wc,
    v_store      = v_store,
    int_p        = int_p,
    int_sig_str  = int_sig_str,
    use_art      = use_art
  )
}

# ---------------------------------------------------------------------------
# PHASE 2: Apply correction across the family of 4, print summary, plot.
# ---------------------------------------------------------------------------
run_analysis12 <- function(info) {
  if (is.null(info)) return(invisible(NULL))
  
  geno_mutant  <- info$geno_mutant
  geno_control <- info$geno_control
  label        <- info$label
  raw_ps       <- info$raw_ps
  use_wc       <- info$use_wc
  int_sig_str  <- info$int_sig_str
  omnibus_type <- if (info$use_art) "ART-ANOVA" else "Two-way ANOVA"
  
  p_holm <- raw_ps; p_bh <- raw_ps
  
  geno_keys <- c("geno_rested", "geno_sleeping")
  cond_keys_corr <- c("cond_mutant", "cond_control")
  
  valid_geno <- !is.na(raw_ps[geno_keys])
  valid_cond <- !is.na(raw_ps[cond_keys_corr])
  
  p_holm[geno_keys[valid_geno]] <- p.adjust(raw_ps[geno_keys[valid_geno]], method = "holm")
  p_bh[geno_keys[valid_geno]]   <- p.adjust(raw_ps[geno_keys[valid_geno]], method = "BH")
  
  p_holm[cond_keys_corr[valid_cond]] <- p.adjust(raw_ps[cond_keys_corr[valid_cond]], method = "holm")
  p_bh[cond_keys_corr[valid_cond]]   <- p.adjust(raw_ps[cond_keys_corr[valid_cond]], method = "BH")
  
  test_labels <- c(
    geno_rested   = paste0(geno_mutant, " vs ", geno_control, " | awake"),
    geno_sleeping = paste0(geno_mutant, " vs ", geno_control, " | sleeping"),
    cond_mutant   = paste0("awake vs sleeping | ", geno_mutant),
    cond_control  = paste0("awake vs sleeping | ", geno_control)
  )
  
  cat(paste0("\n", strrep("=",70), "\n"))
  cat(paste0("ANALYSIS 12", info$suffix, " — ", label,
             " FAMILY-OF-4 CORRECTION SUMMARY\n"))
  cat(sprintf("  %-42s  %-8s  %-10s  %-10s  %-10s  %s\n",
              "Test", "Method", "Raw p", "Holm p", "BH p", "Sig(Holm)"))
  for (nm in names(raw_ps)) {
    if (is.na(raw_ps[nm])) next
    cat(sprintf("  %-42s  %-8s  %-10.4f  %-10.4f  %-10.4f  %s\n",
                test_labels[nm],
                if (isTRUE(use_wc[nm])) "MWU" else "Welch",
                raw_ps[nm], p_holm[nm], p_bh[nm],
                if (!is.na(p_holm[nm]) && p_holm[nm] < 0.05) "*" else "ns"))
  }
  cat(paste0(strrep("=",70), "\n"))
  
  conds_12      <- c("rested", "sleeping")
  box_offset    <- 0.22; bw <- 0.20
  group_centres <- c(1, 3)
  CTRL_BORDER   <- "#1A5CB0"; MUT_BORDER  <- "#D10000"
  REST_BORDER   <- "#1A5CB0"; SLEEP_BORDER <- "#E00000"
  CTRL_FILL     <- adjustcolor(CTRL_BORDER,  alpha.f = 0.55)
  MUT_FILL      <- adjustcolor(MUT_BORDER,   alpha.f = 0.55)
  REST_FILL     <- adjustcolor(REST_BORDER,  alpha.f = 0.55)
  SLEEP_FILL    <- adjustcolor(SLEEP_BORDER, alpha.f = 0.55)
  
  # Helper to draw bracket + significance label on a grouped box plot
  add_bracket <- function(x1, x2, p_h, p_b, y_holm = 95, y_bh = 102, tip = 2) {
    lbl_h <- sig_label(p_h); lbl_b <- sig_label(p_b)
    if (!is.null(lbl_h)) {
      par(xpd = NA)
      segments(x1, y_holm, x2, y_holm, lwd = 1.8)
      segments(x1, y_holm - tip, x1, y_holm, lwd = 1.8)
      segments(x2, y_holm - tip, x2, y_holm, lwd = 1.8)
      text((x1 + x2)/2, y_holm + tip, lbl_h, cex = 1.4, font = 2)
      par(xpd = FALSE)
    } else if (!is.null(lbl_b)) {
      par(xpd = NA)
      segments(x1, y_bh, x2, y_bh, lwd = 1.8, col = "#2980B9")
      segments(x1, y_bh - tip, x1, y_bh, lwd = 1.8, col = "#2980B9")
      segments(x2, y_bh - tip, x2, y_bh, lwd = 1.8, col = "#2980B9")
      text((x1 + x2)/2, y_bh + tip, paste0(lbl_b, "(BH)"),
           cex = 1.3, font = 2, col = "#2980B9")
      par(xpd = FALSE)
    }
  }
  
  footer <- paste0(omnibus_type, " interaction p=", int_sig_str,
                   "  |  Holm/BH family of 4")
  
  # -------------------------------------------------------------------
  # PLOT 1: Genotype contrast within each condition
  # -------------------------------------------------------------------
  x_ctrl <- group_centres - box_offset
  x_mut  <- group_centres + box_offset
  
  par(mfrow = c(1,1), mar = c(8,6,7,2),
      cex = 1.4, cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.3,
      font = 2, font.main = 2, font.lab = 2, font.axis = 2)
  plot(0, 0, type = "n", xlim = c(0.3, 4.2), ylim = c(0, 115),
       xlab = "", ylab = "Volume % Threshold",
       main = paste0(label, " — Genotype Contrast\n",
                     "[Holm family-of-4 | Black=Holm, Blue=BH-only]"),
       xaxt = "n", las = 1, bty = "n")
  par(xpd = NA)
  text(group_centres, y = -10, labels = c("Awake", "Asleep"),
       adj = 0.5, cex = 1.2, font = 2)
  par(xpd = FALSE)
  abline(v = 2, lty = 3, col = "grey80", lwd = 1)
  
  for (i in seq_along(conds_12)) {
    key <- paste0("geno_", conds_12[i])
    cv  <- info$v_store[[key]]$ctrl
    mv  <- info$v_store[[key]]$mut
    if (!is.null(cv)) draw_manual_box(x_ctrl[i], cv, CTRL_FILL, CTRL_BORDER, box_width = bw)
    if (!is.null(mv)) draw_manual_box(x_mut[i],  mv, MUT_FILL,  MUT_BORDER,  box_width = bw)
    add_bracket(x_ctrl[i], x_mut[i], p_holm[key], p_bh[key])
  }
  
  mtext(footer, side = 1, line = 6, cex = 0.95, col = "grey35", font = 2)
  legend("topleft",
         legend = c(paste0(geno_control, " (control)"),
                    paste0(geno_mutant,  " (mutant)"),
                    "Holm sig. (black)", "BH only (blue)"),
         fill   = c(CTRL_FILL, MUT_FILL, NA, NA),
         border = c(CTRL_BORDER, MUT_BORDER, NA, NA),
         lty    = c(NA, NA, 1, 1), col = c(NA, NA, "black", "#2980B9"),
         lwd    = c(NA, NA, 2, 2),
         bty = "n", cex = 1.1, text.font = 2)
  
  # -------------------------------------------------------------------
  # PLOT 2: Condition contrast within each genotype
  # -------------------------------------------------------------------
  genos_12 <- c(geno_control, geno_mutant)
  x_rest   <- group_centres - box_offset
  x_sleep  <- group_centres + box_offset
  cond_keys <- c(
    setNames("cond_control", geno_control),
    setNames("cond_mutant",  geno_mutant)
  )
  
  par(mfrow = c(1,1), mar = c(8,6,7,2),
      cex = 1.4, cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.3,
      font = 2, font.main = 2, font.lab = 2, font.axis = 2)
  plot(0, 0, type = "n", xlim = c(0.3, 4.2), ylim = c(0, 115),
       xlab = "", ylab = "Volume % Threshold",
       main = paste0(label, " — Condition Contrast\n",
                     "[Holm family-of-4 | Black=Holm, Blue=BH-only]"),
       xaxt = "n", las = 1, bty = "n")
  par(xpd = NA)
  text(group_centres, y = -10, labels = genos_12, adj = 0.5, cex = 1.2, font = 2)
  par(xpd = FALSE)
  abline(v = 2, lty = 3, col = "grey80", lwd = 1)
  
  for (i in seq_along(genos_12)) {
    key <- cond_keys[genos_12[i]]
    vr  <- info$v_store[[key]]$rested
    vs  <- info$v_store[[key]]$sleeping
    if (!is.null(vr)) draw_manual_box(x_rest[i],  vr, REST_FILL,  REST_BORDER,  box_width = bw)
    if (!is.null(vs)) draw_manual_box(x_sleep[i], vs, SLEEP_FILL, SLEEP_BORDER, box_width = bw)
    add_bracket(x_rest[i], x_sleep[i], p_holm[key], p_bh[key])
  }
  
  mtext(footer, side = 1, line = 6, cex = 0.95, col = "grey35", font = 2)
  legend("topleft",
         legend = c("Awake", "Asleep",
                    "Holm sig. (black)", "BH only (blue)"),
         fill   = c(REST_FILL, SLEEP_FILL, NA, NA),
         border = c(REST_BORDER, SLEEP_BORDER, NA, NA),
         lty    = c(NA, NA, 1, 1), col = c(NA, NA, "black", "#2980B9"),
         lwd    = c(NA, NA, 2, 2),
         bty = "n", cex = 1.1, text.font = 2)
  
  # -------------------------------------------------------------------
  # PLOT 3: Interaction plot (means +/- SE)
  # -------------------------------------------------------------------
  sub         <- info$sub
  geno_levels <- levels(sub$genotype)
  gc12        <- make_geno_colours(geno_levels)
  summ <- do.call(rbind, lapply(geno_levels, function(g) {
    do.call(rbind, lapply(conds_12, function(cond) {
      vals <- sub$volume_percent_threshold[sub$genotype == g & sub$sleep_cond == cond]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) return(NULL)
      data.frame(genotype = g, sleep_cond = cond,
                 mean = mean(vals), se = sd(vals)/sqrt(length(vals)),
                 stringsAsFactors = FALSE)
    }))
  }))
  
  x_pos <- seq_along(conds_12)
  par(mfrow = c(1,1), mar = c(6,6,5,2),
      cex = 1.4, cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.3,
      font = 2, font.main = 2, font.lab = 2, font.axis = 2)
  plot(0, 0, type = "n", xlim = c(0.7, length(conds_12) + 0.3), ylim = c(0, 100),
       xlab = "", ylab = "Mean Volume % Threshold (+/-SE)",
       main = paste0("Interaction Plot: ", label, " [", omnibus_type, "]"),
       xaxt = "n", las = 1, bty = "n")
  par(xpd = NA)
  text(x_pos, y = -8, labels = c("Awake", "Asleep"),
       adj = 0.5, cex = 1.1, font = 2)
  par(xpd = FALSE)
  for (g in geno_levels) {
    sub_g <- summ[summ$genotype == g, ]
    sub_g <- sub_g[match(conds_12, sub_g$sleep_cond), ]
    arrows(x_pos, sub_g$mean - sub_g$se, x_pos, sub_g$mean + sub_g$se,
           angle = 90, code = 3, length = 0.07, col = gc12$border[g], lwd = 2)
    lines(x_pos, sub_g$mean, col = gc12$border[g], lwd = 2.5,
          lty = if (g == geno_control) 1 else 2)
    points(x_pos, sub_g$mean, pch = if (g == geno_control) 16 else 17,
           col = gc12$border[g], cex = 1.6)
  }
  mtext(paste0(omnibus_type, " interaction p=", int_sig_str),
        side = 1, line = 4.5, cex = 1.0, col = "grey35", font = 2)
  legend("topleft", legend = geno_levels, col = gc12$border[geno_levels],
         lty = c(1,2), pch = c(16,17), lwd = 2, bty = "n", cex = 1.2, text.font = 2)
}

# ---------------------------------------------------------------------------
# Run 12A and 12B — each with its own independent family of 4
# ---------------------------------------------------------------------------
info_hk  <- collect_analysis12("hk",  "canton_s", df)
info_inc <- collect_analysis12("inc", "w1118",    df)

run_analysis12(info_hk)
run_analysis12(info_inc)

cat("\n=== ALL ANALYSES COMPLETE ===\n")
cat("\nMULTIPLE COMPARISONS SUMMARY:\n")
cat("  Analysis 1:    Holm across 4 genotypes (awake vs sleeping family)\n")
cat("  Analysis 2:    KW gate -> Holm pairwise Mann-Whitney U\n")
cat("  Analysis 3:    Omnibus gate -> Holm pairwise\n")
cat("  Analyses 5+10: Holm across latency family of 2\n")
cat("  Analysis 6:    Single test, no correction\n")
cat("  Analyses 7+11: Holm across threshold-sleeping family of 2\n")
cat("  Analysis 8:    Sensitivity check for Analysis 5\n")
cat("  Analysis 12A:  Holm across family of 4 (HK/canton_s: geno x2 + cond x2)\n")
cat("  Analysis 12B:  Holm across family of 4 (inc/w1118:   geno x2 + cond x2)\n")
cat("\nBlack brackets = Holm (PRIMARY) | Blue = BH only (SENSITIVITY)\n")
cat("\nNOTE: shaker genotype excluded from all analyses.\n")