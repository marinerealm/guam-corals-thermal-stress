# =============================================================================
# 02_pam_analysis.R — photosynthetic efficiency (Fv/Fm)
#
# READS   data/processed/PAM_BASELINE_10Mar.csv, PAM_BASELINE.csv,
#         PAM_HEAT.csv, PAM_RECOVERY.csv
#         data/raw/logger/*            (for the cumulative heat dose)
# WRITES  results/tables/Table_PAM_descriptive.csv     (Supplementary Table S2)
#         results/tables/Table_PAM_effect_sizes.csv    (Supplementary Table S5)
#         results/reports/PAM_statistics_report.txt    (Tables S3, S4, S6)
#         results/figures/Figure2_*, Figure3_*, Figure4_*
#         results/figures/FigureS1_raw_fvfm.*         (Supplementary Fig. S1)
#
# REPRODUCES  all Fv/Fm statistics in Results; Supplementary Tables S2-S5;
#             Figures 2, 3 and 4; Supplementary Figure S1.
#
# DESIGN CONSTRAINT
#   One tank per treatment, so fragments within a tank are pseudoreplicates
#   with respect to temperature. Control-versus-heated differences are
#   reported as effect sizes only. Comparisons AMONG SPECIES WITHIN THE
#   HEATED TANK are valid inference, because all three shared the same water.
#
# Run from the scripts/ directory:  Rscript 02_pam_analysis.R
# =============================================================================

source("00_common.R")

# ── 1. LOAD ──────────────────────────────────────────────────────────────────
read_pam <- function(f){
  read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
    mutate(FvFm = suppressWarnings(as.numeric(FvFm)),
           Fo   = suppressWarnings(as.numeric(Fo)),
           Fm   = suppressWarnings(as.numeric(Fm)),
           Date = as.Date(Date))
}
files <- file.path(DIR_PROC, c("PAM_BASELINE_10Mar.csv","PAM_BASELINE.csv",
                               "PAM_HEAT.csv","PAM_RECOVERY.csv"))
miss <- files[!file.exists(files)]
if (length(miss)) stop("Missing PAM file(s):\n  ", paste(miss, collapse = "\n  "))

pam <- bind_rows(lapply(files, read_pam)) |> filter(!is.na(FvFm)) |>
  mutate(
    Elev_meas = unname(ELEV_MAP[as.character(Date)]),
    ExpDay    = as.integer(Date - as.Date("2025-03-10")) + 1L,
    Species   = factor(Species, levels = sp_levels),
    Tank      = factor(Tank, levels = c("Control","Heat")),
    Phase     = case_when(Date <= as.Date("2025-03-11") ~ "Baseline",
                          Date <= as.Date("2025-03-16") ~ "Warming",
                          TRUE                          ~ "Recovery"),
    Phase     = factor(Phase, levels = c("Baseline","Warming","Recovery")))

cat("=====================================================================\n")
cat(" 1. DATA LOADED  (phases defined by realised temperature)\n")
cat("=====================================================================\n")
print(as.data.frame(count(pam, Date, Elev_meas, Phase, name = "n")))

# cumulative heat dose, from the logger record
dose <- dose_by_night(cumulative_dose(load_logger()))

# ── 2. QUALITY CONTROL ───────────────────────────────────────────────────────
cat("\n=====================================================================\n")
cat(" 2. QUALITY CONTROL\n")
cat("=====================================================================\n")
bad <- filter(pam, FvFm < 0.15 | FvFm > 0.80)
cat(sprintf("Implausible readings removed (<0.15 or >0.80): %d\n", nrow(bad)))
if (nrow(bad)) print(as.data.frame(select(bad, Date, Species, Tank, FvFm)))

dat <- filter(pam, FvFm >= 0.15, FvFm <= 0.80) |>
  group_by(Date, Species, Tank) |>
  mutate(Outlier = n() >= 5 &
           (FvFm < quantile(FvFm,.25) - 1.5*IQR(FvFm) |
            FvFm > quantile(FvFm,.75) + 1.5*IQR(FvFm))) |>
  ungroup()
cat(sprintf("IQR flagged (groups of 5+ only): %d\n", sum(dat$Outlier)))
dat <- filter(dat, !Outlier) |> select(-Outlier)
cat(sprintf("Retained: %d readings\n", nrow(dat)))

# ── 3. BASELINE AND NORMALISATION ────────────────────────────────────────────
# Species differ substantially in absolute Fv/Fm for reasons unrelated to
# condition (morphology, tissue optics, fragment geometry), so every
# comparison is made on values normalised to each species-by-tank group's
# own pre-warming baseline.
base <- dat |> filter(Date %in% BASE_DATES) |>
  group_by(Species, Tank) |>
  summarise(Baseline = mean(FvFm), n_base = n(), .groups = "drop")
cat("\nBaseline used for normalisation (", BASELINE_NIGHTS, "):\n")
print(as.data.frame(base |> mutate(Baseline = round(Baseline, 3))))

dat <- left_join(dat, base, by = c("Species","Tank")) |>
  mutate(Pct = 100 * FvFm / Baseline)

# ── 4. DESCRIPTIVE TABLES ────────────────────────────────────────────────────
desc <- dat |>
  group_by(Date, ExpDay, Elev_meas, Phase, Species, Tank) |>
  summarise(n = n(), Mean = mean(FvFm), SD = sd(FvFm), SE = sd(FvFm)/sqrt(n()),
            Pct_mean = mean(Pct), Pct_SE = sd(Pct)/sqrt(n()), .groups = "drop")
write_csv(desc |> mutate(across(where(is.numeric), \(x) round(x, 3))),
          file.path(DIR_TAB, "Table_PAM_descriptive.csv"))

eff <- desc |>
  select(Date, Elev_meas, Phase, Species, Tank, Mean, Pct_mean, n) |>
  pivot_wider(names_from = Tank, values_from = c(Mean, Pct_mean, n)) |>
  mutate(Diff_FvFm = Mean_Heat - Mean_Control,
         Diff_Pct  = Pct_mean_Heat - Pct_mean_Control) |>
  left_join(dose, by = "Date")
write_csv(eff |> mutate(across(where(is.numeric), \(x) round(x, 3))),
          file.path(DIR_TAB, "Table_PAM_effect_sizes.csv"))

# ── 5. STATISTICS ────────────────────────────────────────────────────────────
sink(file.path(DIR_REP, "PAM_statistics_report.txt"), split = TRUE)
cat("=====================================================================\n")
cat(" PAM STATISTICAL REPORT  ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat(" Baseline nights: ", BASELINE_NIGHTS, "\n")
cat("=====================================================================\n")
cat("\nDESIGN CONSTRAINT\n")
cat("One tank per treatment, so fragments within a tank are not independent\n")
cat("with respect to temperature. Control-vs-heat differences are reported as\n")
cat("effect sizes only. Between-species comparisons within the heated tank are\n")
cat("valid, because all species shared the same water.\n")

cat("\n--- 5a. Baseline: species differences and tank offset [Table S3] ---\n")
bl <- filter(dat, Date %in% BASE_DATES)
print(summary(aov(FvFm ~ Species, data = bl)))
cat("\nTukey HSD:\n"); print(round(TukeyHSD(aov(FvFm ~ Species, data = bl))$Species, 4))
cat("\nTank offset within each species at baseline (should be near zero):\n")
for (s in sp_levels){
  x <- filter(bl, Species == s)
  if (n_distinct(x$Tank) < 2 || min(table(x$Tank)) < 2) next
  tt <- t.test(FvFm ~ Tank, data = x)
  cat(sprintf("  %-14s control %.3f vs heat %.3f  diff %+.3f  p = %.3f\n",
              s, tt$estimate[1], tt$estimate[2],
              tt$estimate[2] - tt$estimate[1], tt$p.value))
}

cat("\n--- 5b. Control stability across the experiment ---\n")
for (s in sp_levels){
  x <- filter(dat, Species == s, Tank == "Control")
  if (nrow(x) < 6) next
  ct <- summary(lm(Pct ~ ExpDay, data = x))$coefficients
  cat(sprintf("  %-14s slope %+6.2f %%/day  p = %.3f  (%s)\n", s,
              ct["ExpDay","Estimate"], ct["ExpDay","Pr(>|t|)"],
              ifelse(ct["ExpDay","Pr(>|t|)"] > 0.05, "stable", "DRIFTING")))
}

cat("\n--- 5c. BETWEEN-SPECIES WITHIN HEATED TANK  [VALID INFERENCE] ---\n")
hw <- filter(dat, Tank == "Heat", Phase == "Warming")
cat(sprintf("Pooled across warming nights (n = %d):\n", nrow(hw)))
kw <- kruskal.test(Pct ~ Species, data = hw)
cat(sprintf("  Kruskal-Wallis H = %.3f, df = %d, p = %.4f\n",
            kw$statistic, kw$parameter, kw$p.value))
a1 <- aov(Pct ~ Species, data = hw); print(summary(a1))
cat("\nTukey HSD:\n"); print(round(TukeyHSD(a1)$Species, 4))

cat("\nPer warming night:\n")
print(as.data.frame(map_dfr(sort(unique(hw$Date)), function(d){
  s <- filter(hw, Date == d); if (n_distinct(s$Species) < 2) return(NULL)
  a <- summary(aov(Pct ~ Species, data = s))[[1]]
  tibble(Date = d, n = nrow(s), F = round(a[["F value"]][1], 3),
         p = round(a[["Pr(>F)"]][1], 4))
})))

rec <- filter(dat, Tank == "Heat", Phase == "Recovery")
if (nrow(rec) > 6){
  cat("\nRecovery night, between species:\n")
  a2 <- aov(Pct ~ Species, data = rec); print(summary(a2))
  cat("\nTukey HSD:\n"); print(round(TukeyHSD(a2)$Species, 4))
}

cat("\n--- 5d. Trajectory within each species and tank [Table S4] ---\n")
print(as.data.frame(map_dfr(sp_levels, function(s){
  map_dfr(c("Control","Heat"), function(tk){
    x <- filter(dat, Species == s, Tank == tk, Phase != "Baseline")
    if (nrow(x) < 6) return(NULL)
    ct <- summary(lm(Pct ~ ExpDay, data = x))$coefficients
    tibble(Species = s, Tank = tk, n = nrow(x),
           slope = round(ct["ExpDay","Estimate"], 3),
           p = round(ct["ExpDay","Pr(>|t|)"], 4))
  })})))

cat("\n--- 5e. Effect size vs cumulative excess heat [Table S5] ---\n")
cat("Heated minus control, in percentage points of each group's own baseline,\n")
cat("against cumulative excess heat. Reported as effect sizes: with one tank\n")
cat("per treatment this difference cannot be tested inferentially.\n\n")
print(as.data.frame(eff |>
  select(Date, Elev_meas, Dose, Species, Diff_FvFm, Diff_Pct) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))))

# The regression spans 11-16 March. The pre-warming night of 11 March has zero
# dose and is retained because it anchors the relationship at zero exposure.
# The recovery night (17 March) is excluded: cumulative exposure keeps rising
# after the heater is switched off, while the fragments are no longer heated,
# so dose cannot predict the recovery state.
DOSE_SET <- eff |> filter(Date >= as.Date("2025-03-11"), Date <= as.Date("2025-03-16"))
cat("\nDoes the effect scale with cumulative excess heat? (11-16 March, n = 6)\n")
dose_fit <- map_dfr(sp_levels, function(s){
  x <- filter(DOSE_SET, Species == s)
  if (nrow(x) < 3) return(NULL)
  m  <- lm(Diff_Pct ~ Dose, data = x); sm <- summary(m)
  cat(sprintf("  %-14s %+6.2f pp per C-day  p = %.3f  R2 = %.2f  intercept %+5.1f\n",
              s, coef(m)[2], coef(sm)[2,4], sm$r.squared, coef(m)[1]))
  tibble(Species = s, slope = coef(m)[2], p = coef(sm)[2,4], r2 = sm$r.squared)
})
r_ed <- with(DOSE_SET |> filter(Species == "Acropora"), cor(Elev_meas, Dose))
cat(sprintf("\nCorrelation between nightly elevation and cumulative dose: r = %.2f\n", r_ed))
cat("The tank warmed monotonically, so magnitude and duration of warming are\n")
cat("confounded and these data cannot separate their effects.\n")
sink()

# ── 6. VERIFICATION ──────────────────────────────────────────────────────────
cat("\n--- VERIFICATION against the published manuscript ---\n")
check("baseline species ANOVA F(2,83)",
      summary(aov(FvFm ~ Species, data = bl))[[1]][["F value"]][1], PUB$baseline_F, 0.05)
check("between-species warming F(2,137)",
      summary(a1)[[1]][["F value"]][1], PUB$warming_F, 0.05)
check("between-species warming Kruskal-Wallis H",
      unname(kw$statistic), PUB$warming_KW, 0.05)
# Tolerances are wide enough to absorb rounding in the published values, but
# far tighter than any structural error (wrong QC, wrong dates, wrong time
# zone) would produce — those shift these numbers by whole units.
check("dose-response slope, P. decussata",
      dose_fit$slope[dose_fit$Species == "Pavona"], PUB$dose_slope_pav, 0.15)
check("dose-response p, P. decussata",
      dose_fit$p[dose_fit$Species == "Pavona"], PUB$dose_p_pav, 0.004)
cat("  all PAM checks passed\n\n")

# ── 7. FIGURES ───────────────────────────────────────────────────────────────
nightly_pam <- dat |> group_by(Date, Species, Tank) |>
  summarise(n = n(), Mean = mean(Pct), SE = sd(Pct)/sqrt(n()), .groups = "drop")

dfrac <- function(s, off = 0) structure(as.numeric(as.Date(s)) + off, class = "Date")
axis_lab <- nightly_pam |> distinct(Date) |> arrange(Date) |>
  left_join(dose, by = "Date") |>
  mutate(Elev = unname(ELEV_MAP[as.character(Date)]),
         Label = paste0(format(Date, "%d %b"), "\n",
                        ifelse(is.na(Elev), "(no logger)", sprintf("+%.1f°C", Elev)), "\n",
                        ifelse(is.na(Dose), "—", sprintf("%.1f °C-d", Dose))))
phase_bands <- tibble(
  Phase = factor(c("Baseline","Thermal Stress","Recovery"),
                 levels = c("Baseline","Thermal Stress","Recovery")),
  xmin = c(dfrac("2025-03-10", -0.5), dfrac("2025-03-11", 0.5), dfrac("2025-03-16", 0.5)),
  xmax = c(dfrac("2025-03-11",  0.5), dfrac("2025-03-16", 0.5), dfrac("2025-03-17", 0.5)))
phase_bands$xmid <- structure((as.numeric(phase_bands$xmin) +
                               as.numeric(phase_bands$xmax))/2, class = "Date")

# --- Figure 2: time course, percentage of own baseline ---
fig2 <- ggplot() +
  geom_rect(data = phase_bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Phase),
            inherit.aes = FALSE, alpha = 0.55, show.legend = FALSE) +
  scale_fill_manual(values = phase_cols) +
  geom_vline(xintercept = c(dfrac("2025-03-11", 0.5), dfrac("2025-03-16", 0.5)),
             linetype = "dotted", colour = "grey45", linewidth = 0.45) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey40") +
  geom_text(data = phase_bands, aes(x = xmid, y = Inf, label = Phase),
            inherit.aes = FALSE, vjust = 1.8, size = 2.9,
            fontface = "italic", colour = "grey35") +
  geom_point(data = dat, aes(Date, Pct, colour = Tank), inherit.aes = FALSE,
             position = position_jitter(width = 0.13, height = 0),
             size = 0.75, alpha = 0.33, show.legend = FALSE) +
  geom_line(data = nightly_pam, aes(Date, Mean, colour = Tank, group = Tank),
            linewidth = 0.95, inherit.aes = FALSE) +
  geom_errorbar(data = nightly_pam, aes(Date, Mean, ymin = Mean-SE, ymax = Mean+SE),
                width = 0.1, linewidth = 0.5, colour = "grey25", inherit.aes = FALSE) +
  geom_point(data = nightly_pam, aes(Date, Mean, colour = Tank),
             size = 3.0, shape = 21, fill = "grey45", stroke = 1.0, inherit.aes = FALSE) +
  scale_colour_manual(values = tank_cols, labels = c("Control","Heat")) +
  scale_x_date(breaks = axis_lab$Date, labels = axis_lab$Label,
               limits = c(dfrac("2025-03-10", -0.5), dfrac("2025-03-17", 0.5)),
               expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  facet_wrap(~Species, nrow = 1, labeller = labeller(Species = sp_labeller)) +
  labs(x = NULL, y = LAB_FVFM) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2, size = 3))) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6.6, lineheight = 1.0))
save_fig(fig2, "Figure2_pct_baseline", w = 13, h = 5.6)

# --- Supplementary Figure S1: the same series without baseline normalisation ---
# Species differ several-fold in absolute Fv/Fm, which is why the main figures
# are normalised; S1 shows the raw values so that the normalisation can be
# checked. Date labels carry the temperature elevation but not the dose.
nightly_raw <- dat |> group_by(Date, Species, Tank) |>
  summarise(n = n(), Mean = mean(FvFm), SE = sd(FvFm)/sqrt(n()), .groups = "drop")
axis_lab_raw <- axis_lab |>
  mutate(Label = paste0(format(Date, "%d %b"), "\n",
                        ifelse(is.na(Elev), "(no logger)", sprintf("+%.1f°C", Elev))))

figS1 <- ggplot() +
  geom_rect(data = phase_bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Phase),
            inherit.aes = FALSE, alpha = 0.55, show.legend = FALSE) +
  scale_fill_manual(values = phase_cols) +
  geom_vline(xintercept = c(dfrac("2025-03-11", 0.5), dfrac("2025-03-16", 0.5)),
             linetype = "dotted", colour = "grey45", linewidth = 0.45) +
  geom_text(data = phase_bands, aes(x = xmid, y = Inf, label = Phase),
            inherit.aes = FALSE, vjust = 1.8, size = 2.9,
            fontface = "italic", colour = "grey35") +
  geom_point(data = dat, aes(Date, FvFm, colour = Tank), inherit.aes = FALSE,
             position = position_jitter(width = 0.13, height = 0),
             size = 0.75, alpha = 0.33, show.legend = FALSE) +
  geom_line(data = nightly_raw, aes(Date, Mean, colour = Tank, group = Tank),
            linewidth = 0.95, inherit.aes = FALSE) +
  geom_errorbar(data = nightly_raw, aes(Date, Mean, ymin = Mean-SE, ymax = Mean+SE),
                width = 0.1, linewidth = 0.5, colour = "grey25", inherit.aes = FALSE) +
  geom_point(data = nightly_raw, aes(Date, Mean, colour = Tank),
             size = 3.0, shape = 21, fill = "grey45", stroke = 1.0, inherit.aes = FALSE) +
  scale_colour_manual(values = tank_cols, labels = c("Control","Heat")) +
  scale_x_date(breaks = axis_lab_raw$Date, labels = axis_lab_raw$Label,
               limits = c(dfrac("2025-03-10", -0.5), dfrac("2025-03-17", 0.5)),
               expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  facet_wrap(~Species, nrow = 1, labeller = labeller(Species = sp_labeller)) +
  labs(x = NULL, y = LAB_FVFM_RAW) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2, size = 3))) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6.6, lineheight = 1.0))
save_fig(figS1, "FigureS1_raw_fvfm", w = 13, h = 5.6)

# --- Figure 3: between-species comparison, heated tank, warming period ---
tuk <- as.data.frame(TukeyHSD(a1)$Species)
cat("\nFigure 3 Tukey contrasts:\n"); print(round(tuk, 4))
star <- function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else
                    if (p < 0.05) "*" else NA_character_
# Bracket label as plotmath, so that p is italic and "< 0.001" is not rounded
# to "= 0.000": paste() joins the quoted text to the italic symbol.
plab <- function(st, p) if (p < 0.001)
    sprintf('paste("%s  ", italic(p) < 0.001)', st) else
    sprintf('paste("%s  ", italic(p) == %.3f)', st, p)
ymax <- max(hw$Pct, na.rm = TRUE)
brk <- tibble(contrast = rownames(tuk), p = tuk$`p adj`) |>
  filter(!is.na(sapply(p, star))) |>
  mutate(lab = mapply(plab, sapply(p, star), p),
         g1 = sub("-.*", "", contrast), g2 = sub(".*-", "", contrast),
         x = match(g2, sp_levels), xend = match(g1, sp_levels)) |>
  arrange(abs(xend - x)) |>
  mutate(y = ymax + seq_len(n()) * 0.07 * diff(range(hw$Pct, na.rm = TRUE)))

# Omnibus ANOVA annotation in the top-left corner
a1s  <- summary(a1)[[1]]
a1_F <- a1s[["F value"]][1]; a1_p <- a1s[["Pr(>F)"]][1]; a1_df <- a1s[["Df"]][2]
anova_lab <- sprintf('paste("ANOVA  ", italic(F)[2*","*%d], " = %.2f,  ", %s)',
                     a1_df, a1_F,
                     if (a1_p < 0.001) 'italic(p) < 0.001'
                     else sprintf('italic(p) == %.3f', a1_p))

fig3 <- ggplot(hw, aes(Species, Pct, colour = Species, fill = Species)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey45") +
  geom_boxplot(width = 0.55, alpha = 0.18, outlier.shape = NA, linewidth = 0.6) +
  ggbeeswarm::geom_quasirandom(width = 0.18, size = 1.1, alpha = 0.55, shape = 16) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2,
               fill = "white", stroke = 1.1) +
  { if (nrow(brk))
      list(geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = y),
                        inherit.aes = FALSE, colour = "grey25", linewidth = 0.4),
           geom_text(data = brk, aes(x = (x+xend)/2, y = y, label = lab),
                     inherit.aes = FALSE, vjust = -0.3, size = 3.1,
                     colour = "grey15", parse = TRUE))
    else NULL } +
  annotate("text", x = 0.55, y = ymax + 0.26 * diff(range(hw$Pct, na.rm = TRUE)),
           label = anova_lab, parse = TRUE, hjust = 0, size = 3.1,
           fontface = "italic", colour = "grey25") +
  scale_colour_manual(values = sp_cols) + scale_fill_manual(values = sp_cols) +
  scale_x_discrete(labels = sp_axis_labels()) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = NULL, y = LAB_FVFM) +
  theme_pub() +
  # no face = "italic" here: the labels are plotmath, which sets its own
  # italics and keeps "cf." and "sp." upright
  theme(legend.position = "none")
save_fig(fig3, "Figure3_between_species", w = 7.5, h = 5.5)

# --- Figure 4: effect size against cumulative excess heat ---
# Per-panel fit statistics for Figure 4. A reader cannot judge a regression
# line without them: Sclerophytum's line is visibly steeper than Acropora's
# yet neither is significant, and the slope is a difference between tanks,
# not a decline of the heated fragments (see the report above).
dose_fits <- map_dfr(sp_levels, function(s){
  x  <- filter(DOSE_SET, Species == s)
  m  <- lm(Diff_Pct ~ Dose, data = x); sm <- summary(m)
  pv <- coef(sm)[2, 4]
  tibble(Species = factor(s, levels = sp_levels),
         lab = sprintf('paste(italic(R)^2, " = %.2f,  ", italic(p), " %s")',
                       sm$r.squared,
                       if (pv < 0.001) "< 0.001" else sprintf("= %.3f", pv)))
})

fig4 <- ggplot(DOSE_SET, aes(Dose, Diff_Pct)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_smooth(method = "lm", se = TRUE, colour = NA, fill = "grey80", alpha = 0.55) +
  geom_smooth(aes(colour = Species), method = "lm", se = FALSE, linewidth = 1.0) +
  geom_point(aes(colour = Species), size = 2.8, shape = 21,
             fill = "white", stroke = 1.3) +
  geom_text(data = dose_fits, aes(x = -Inf, y = Inf, label = lab),
            inherit.aes = FALSE, parse = TRUE, hjust = -0.06, vjust = 1.5,
            size = 3.0, colour = "grey25") +
  scale_colour_manual(values = sp_cols) +
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15),
                     expand = expansion(mult = 0.07)) +
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.16))) +
  facet_wrap(~Species, nrow = 1, labeller = labeller(Species = sp_labeller)) +
  labs(x = expression("Cumulative excess heat ("*degree*"C-days)"),
       # plotmath, so the minus comes from the maths font: a literal U+2212
       # is silently dropped by some graphics devices whose text font
       # lacks that glyph, leaving "Heat ... control"
       y = expression("Heat" - "control (percentage points)")) +
  theme_pub() + theme(legend.position = "none")
save_fig(fig4, "Figure4_effect_vs_dose", w = 12, h = 4.6)

cat("02_pam_analysis.R complete.\n")
