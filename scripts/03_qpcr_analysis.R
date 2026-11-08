# =============================================================================
# 03_qpcr_analysis.R — symbiont abundance by quantitative PCR
#
# READS   data/raw/qpcr/qPCR_Cq_values.csv
# WRITES  results/tables/Table1_qPCR_foldchange.csv       (Table 1)
#         results/tables/qPCR_per_sample_results.csv
#         results/tables/qPCR_between_species_pairwise.csv
#         results/reports/qPCR_statistics_report.txt
#         results/figures/Figure5_qPCR_foldchange.*       (Figure 5)
#
# REPRODUCES  Table 1 and all qPCR statistics in Results; Figure 5.
#
# METHOD  Livak 2^-ddCq, with the coral host actin gene Cni_ACTB as reference
#         and the species-matched control group as calibrator. All inferential
#         tests are run on ddCq, never on fold change: fold change is
#         log-normally distributed and testing it directly inflates false
#         positives.
#
# DESIGN CONSTRAINT
#   One tank per treatment. Control-versus-heat p-values test whether two
#   tanks differed, not whether temperature caused the difference, and are
#   reported as descriptive effect sizes only. Comparisons AMONG SPECIES
#   WITHIN THE HEATED TANK are valid inference.
#
# Run from the scripts/ directory:  Rscript 03_qpcr_analysis.R
# =============================================================================

source("00_common.R")

REF_GENE      <- "Cni_ACTB"
TARGETS       <- c("Endozoicomonas","Cladocopium","Durusdinium")
CONTROL_LEVEL <- "Control"
MAX_TECH_DIFF <- 0.5   # flag technical replicates differing by more than this
CQ_WARN       <- 32    # at or beyond the limit of reliable quantification
CQ_EXCLUDE    <- 36    # treated as non-detects

cq_file <- file.path(DIR_RAW, "qpcr", "qPCR_Cq_values.csv")
if (!file.exists(cq_file)) stop("Missing: ", cq_file)

raw <- read_csv(cq_file, show_col_types = FALSE) |>
  mutate(across(c(Sample_ID, Species, Treatment, Target), \(x) trimws(as.character(x))),
         Cq = suppressWarnings(as.numeric(Cq)))

# ── 1. TECHNICAL REPLICATES ──────────────────────────────────────────────────
cat("=====================================================================\n")
cat(" 1. TECHNICAL REPLICATES\n")
cat("=====================================================================\n")
tech <- raw |> filter(!is.na(Cq)) |>
  group_by(Sample_ID, Species, Treatment, Target) |>
  summarise(n_tech = n(), mean_Cq = mean(Cq), range_Cq = diff(range(Cq)),
            .groups = "drop") |>
  mutate(tech_flag = range_Cq > MAX_TECH_DIFF)
cat(sprintf("Reactions: %d | samples x targets: %d | replicates differing by >%.1f Cq: %d\n",
            nrow(raw), nrow(tech), MAX_TECH_DIFF, sum(tech$tech_flag)))
if (any(tech$tech_flag))
  print(as.data.frame(filter(tech, tech_flag) |>
          mutate(across(where(is.numeric), \(x) round(x, 3)))))

cat(sprintf("\nLimit of quantification (warn > %d, exclude > %d):\n", CQ_WARN, CQ_EXCLUDE))
loq <- tech |> mutate(LOQ_flag = case_when(mean_Cq > CQ_EXCLUDE ~ "EXCLUDE",
                                           mean_Cq > CQ_WARN    ~ "near LOQ",
                                           TRUE                 ~ "ok"))
print(as.data.frame(count(loq, Target, LOQ_flag)))

# ── 2. LIVAK 2^-ddCq ─────────────────────────────────────────────────────────
sample_cq <- tech |> select(Sample_ID, Species, Treatment, Target,
                            Cq = mean_Cq, n_tech, range_Cq)
ref <- sample_cq |> filter(Target == REF_GENE) |> select(Sample_ID, Cq_ref = Cq)

dcq <- sample_cq |> filter(Target %in% TARGETS) |>
  left_join(ref, by = "Sample_ID") |>
  mutate(dCq = Cq - Cq_ref,
         LOQ_flag = case_when(Cq > CQ_EXCLUDE ~ "EXCLUDE",
                              Cq > CQ_WARN    ~ "near LOQ",
                              TRUE            ~ "ok")) |>
  filter(!is.na(dCq))

# Calibrator is the control group of the SAME species: baseline symbiont load
# differs greatly between species, so cross-species fold change is meaningless.
calib <- dcq |> filter(Treatment == CONTROL_LEVEL) |>
  group_by(Species, Target) |>
  summarise(dCq_cal = mean(dCq), n_cal = n(), .groups = "drop")

res <- dcq |> left_join(calib, by = c("Species","Target")) |>
  mutate(ddCq = dCq - dCq_cal, FoldChange = 2^(-ddCq),
         Species = factor(Species, levels = sp_levels),
         Target  = factor(Target,  levels = TARGETS),
         Treatment = factor(Treatment, levels = c(CONTROL_LEVEL, "Heat")))
write_csv(res |> mutate(across(where(is.numeric), \(x) round(x, 4))),
          file.path(DIR_TAB, "qPCR_per_sample_results.csv"))

# Fold-change bounds come from transforming the ddCq +/- SD limits. Taking the
# SD of fold-change values directly is wrong: the distribution is log-normal,
# so the error bars are inherently asymmetric.
summ <- res |> group_by(Species, Treatment, Target) |>
  summarise(n = n(), mean_ddCq = mean(ddCq),
            sd_ddCq = ifelse(n() > 1, sd(ddCq), NA_real_),
            n_flagged = sum(LOQ_flag != "ok"), .groups = "drop") |>
  mutate(FoldChange = 2^(-mean_ddCq),
         FC_lo_SD   = 2^(-(mean_ddCq + sd_ddCq)),
         FC_hi_SD   = 2^(-(mean_ddCq - sd_ddCq)),
         Reliability = case_when(n_flagged == n ~ "UNRELIABLE - all near/over LOQ",
                                 n_flagged > 0  ~ "caution - some near LOQ",
                                 TRUE           ~ "ok"))

# ── 3. STATISTICS ────────────────────────────────────────────────────────────
heat <- filter(res, Treatment != CONTROL_LEVEL)

omnibus <- map_dfr(TARGETS, function(tg){
  sub <- filter(heat, Target == tg)
  if (n_distinct(sub$Species) < 2) return(NULL)
  a  <- summary(aov(ddCq ~ Species, data = sub))[[1]]
  kwt <- kruskal.test(ddCq ~ Species, data = sub)
  tibble(Target = tg, df1 = a[["Df"]][1], df2 = a[["Df"]][2],
         F = a[["F value"]][1], p = a[["Pr(>F)"]][1],
         KW_H = unname(kwt$statistic), KW_p = kwt$p.value)
})

pairwise_species <- map_dfr(TARGETS, function(tg){
  sub <- filter(heat, Target == tg)
  sps <- sort(unique(as.character(sub$Species)))
  if (length(sps) < 2) return(NULL)
  out <- map_dfr(combn(sps, 2, simplify = FALSE), function(pr){
    a <- sub$ddCq[sub$Species == pr[1]]; b <- sub$ddCq[sub$Species == pr[2]]
    if (length(a) < 2 || length(b) < 2) return(NULL)
    tt <- t.test(a, b)
    tibble(Target = tg, Species_A = pr[1], Species_B = pr[2],
           n_A = length(a), n_B = length(b),
           mean_ddCq_A = mean(a), mean_ddCq_B = mean(b), p_raw = tt$p.value)
  })
  if (nrow(out) == 0) return(NULL)
  mutate(out, p_BH = p.adjust(p_raw, method = "BH"))
}) |>
  mutate(signif = case_when(p_BH < 0.001 ~ "***", p_BH < 0.01 ~ "**",
                            p_BH < 0.05 ~ "*", TRUE ~ "ns"))
write_csv(pairwise_species |> mutate(across(where(is.numeric), \(x) round(x, 4))),
          file.path(DIR_TAB, "qPCR_between_species_pairwise.csv"))

tab1 <- summ |> filter(Treatment != CONTROL_LEVEL) |>
  select(Target, Species, n, ddCq = mean_ddCq, FoldChange, Reliability) |>
  left_join(omnibus |> select(Target, F, p), by = "Target") |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  arrange(Target, Species)
write_csv(tab1, file.path(DIR_TAB, "Table1_qPCR_foldchange.csv"))

sink(file.path(DIR_REP, "qPCR_statistics_report.txt"), split = TRUE)
cat("=====================================================================\n")
cat(" qPCR STATISTICAL REPORT  ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("=====================================================================\n\n")
cat("SECTION 1. CONTROL vs HEAT -- DESCRIPTIVE ONLY\n")
cat("---------------------------------------------------------------------\n")
cat("One tank per treatment, so fragments within a tank are pseudoreplicates\n")
cat("with respect to temperature. A control-vs-heat p-value here tests whether\n")
cat("two tanks differed, not whether temperature caused the difference. Report\n")
cat("the fold changes; do NOT report such p-values as evidence of a\n")
cat("temperature effect.\n\n")
print(as.data.frame(summ |> filter(Treatment != CONTROL_LEVEL) |>
        select(Species, Target, n, mean_ddCq, sd_ddCq, FoldChange,
               FC_lo_SD, FC_hi_SD, Reliability) |>
        mutate(across(where(is.numeric), \(x) round(x, 3)))))

cat("\n\nSECTION 2. BETWEEN SPECIES, HEATED TANK -- VALID INFERENCE\n")
cat("---------------------------------------------------------------------\n")
cat("All three species were held in the same heated tank, so tank identity is\n")
cat("held constant and species differences are not confounded by tank. These\n")
cat("tests are the basis for the statistical claims in the manuscript.\n")
cat("Tests are run on ddCq, not on fold change.\n\n")
print(as.data.frame(omnibus |> mutate(across(where(is.numeric), \(x) round(x, 4)))))

cat("\n\nSECTION 3. PAIRWISE, Welch t-test with Benjamini-Hochberg correction\n")
cat("---------------------------------------------------------------------\n")
print(as.data.frame(pairwise_species |>
        mutate(across(where(is.numeric), \(x) round(x, 4)))))

cat("\n\nSECTION 4. LIMIT OF QUANTIFICATION\n")
cat("---------------------------------------------------------------------\n")
cat("Endozoicomonas and Durusdinium in A. cf. manni returned Cq 32-35, at or\n")
cat("beyond reliable quantification. Those fold changes are relative patterns,\n")
cat("not precise estimates. Sclerophytum sp. is n = 2 per group throughout,\n")
cat("following mortality of two of the four colonies.\n")
sink()

# ── 4. VERIFICATION ──────────────────────────────────────────────────────────
cat("\n--- VERIFICATION against the published manuscript ---\n")
gf <- function(tg) omnibus$F[omnibus$Target == tg]
check("Endozoicomonas omnibus F(2,7)", gf("Endozoicomonas"), PUB$qpcr_endo_F, 0.05)
check("Durusdinium omnibus F(2,7)",    gf("Durusdinium"),    PUB$qpcr_duru_F, 0.05)
check("Cladocopium omnibus F(2,7)",    gf("Cladocopium"),    PUB$qpcr_clad_F, 0.05)
fc <- function(sp, tg) summ$FoldChange[summ$Species == sp & summ$Target == tg &
                                       summ$Treatment == "Heat"]
check("Endozoicomonas fold change, Sclerophytum sp.",
      fc("Sclerophytum","Endozoicomonas"), 6.41, 0.02)
check("Durusdinium fold change, A. cf. manni",
      fc("Acropora","Durusdinium"), 4.21, 0.02)
cat("  all qPCR checks passed\n\n")

# ── 5. FIGURE 5 ──────────────────────────────────────────────────────────────
# Heated fragments only, against the control reference line at 1. Controls are
# the calibrator and sit at unity by definition, so plotting them alongside
# would invite exactly the control-vs-heat comparison the design cannot support.
fig_dat  <- heat |> mutate(Species = factor(as.character(Species), levels = sp_levels))
fig_summ <- summ |> filter(Treatment != CONTROL_LEVEL) |>
  mutate(Species = factor(as.character(Species), levels = sp_levels))

# vertical room for the brackets above and the n labels below, on a log2 axis
y_hi <- max(c(fig_dat$FoldChange, fig_summ$FC_hi_SD), na.rm = TRUE)
y_lo <- min(c(fig_dat$FoldChange, fig_summ$FC_lo_SD), na.rm = TRUE)

# significant between-species contrasts, drawn as brackets inside each panel
brk5 <- pairwise_species |> filter(p_BH < 0.05) |>
  mutate(x = match(Species_A, sp_levels), xend = match(Species_B, sp_levels)) |>
  group_by(Target) |> arrange(abs(xend - x), .by_group = TRUE) |>
  mutate(y = y_hi * 2^(0.55 * row_number())) |> ungroup()

# omnibus ANOVA per target, as plotmath so that p is italic
anno5 <- omnibus |>
  mutate(lab = sprintf('paste("ANOVA  ", italic(F), " = %.2f,  ", %s)', F,
                       ifelse(p < 0.001, 'italic(p) < 0.001',
                              sprintf('italic(p) == %.3f', p))),
         y = y_hi * 2^2.0)

# group sizes, which differ because two Sclerophytum colonies were lost
nlab5 <- fig_summ |> transmute(Target, Species, lab = paste0("n=", n),
                               y = y_lo * 2^-1.0)

fig5 <- ggplot(fig_dat, aes(Species, FoldChange, colour = Species)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  ggbeeswarm::geom_quasirandom(width = 0.16, size = 1.6, alpha = 0.6, shape = 16) +
  geom_errorbar(data = fig_summ, aes(x = Species, y = FoldChange,
                                     ymin = FC_lo_SD, ymax = FC_hi_SD),
                width = 0.14, linewidth = 0.5, colour = "grey25", inherit.aes = FALSE) +
  geom_point(data = fig_summ, aes(Species, FoldChange, fill = Species),
             size = 4.0, shape = 21, colour = "grey15", stroke = 0.9,
             inherit.aes = FALSE) +
  { if (nrow(brk5))
      list(geom_segment(data = brk5, aes(x = x, xend = xend, y = y, yend = y),
                        inherit.aes = FALSE, colour = "grey25", linewidth = 0.4),
           geom_text(data = brk5, aes(x = (x + xend)/2, y = y, label = signif),
                     inherit.aes = FALSE, vjust = -0.25, size = 4,
                     colour = "grey15"))
    else NULL } +
  geom_text(data = anno5, aes(x = 2, y = y, label = lab), inherit.aes = FALSE,
            parse = TRUE, size = 3.1, colour = "grey25") +
  geom_text(data = nlab5, aes(x = Species, y = y, label = lab), inherit.aes = FALSE,
            size = 3.0, colour = "grey45") +
  scale_colour_manual(values = sp_cols) + scale_fill_manual(values = sp_cols) +
  scale_x_discrete(labels = sp_axis_labels_2line()) +
  # explicit labels: the default formatter pads every break to the width of
  # the most precise one, giving "8.000" instead of "8"
  scale_y_continuous(trans = "log2", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8),
                     labels = c("0.125","0.25","0.5","1","2","4","8")) +
  coord_cartesian(ylim = c(y_lo * 2^-1.6, y_hi * 2^2.5)) +
  facet_wrap(~Target, nrow = 1) +
  labs(x = NULL,
       y = expression("Fold change relative to control (2"^-Delta*Delta*"Cq)")) +
  theme_pub() +
  # species names are plotmath, so no face = "italic" here
  theme(legend.position = "none",
        strip.text = element_text(face = "bold.italic"))
save_fig(fig5, "Figure5_qPCR_foldchange", w = 10, h = 4.8)

cat("03_qpcr_analysis.R complete.\n")
