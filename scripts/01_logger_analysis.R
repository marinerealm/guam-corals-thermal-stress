# =============================================================================
# 01_logger_analysis.R — thermal history and cumulative heat exposure
#
# READS   data/raw/logger/VK_Control_C.csv, VK_Treatment_T.csv
# WRITES  results/tables/TableS1_nightly_temperature.csv   (Supplementary Table S1)
#         results/tables/dose_by_night.csv
#         results/reports/logger_report.txt
#         results/figures/Figure1_thermal_history_and_dose.*  (Figure 1)
#         results/figures/FigureS3_complete_logger_record.*   (Supplementary Fig. S3)
#
# REPRODUCES  Supplementary Table S1; the thermal values in Results
#             ("Realised thermal conditions"); Figure 1; Supplementary Fig. S3.
#
# Run from the scripts/ directory:  Rscript 01_logger_analysis.R
#
# FIGURE 1 AND FIGURE S3 ARE THE SAME DATA AT TWO SPANS
#   Figure 1 ends on 19 March (PANEL_LAST in 00_common.R): the experiment plus
#   one day of convergence, which is enough to show the tanks were thermally
#   equivalent once the heater was off. Figure S3 plots the complete record to
#   27 March. Both carry the cumulative excess heat curve on a right-hand axis.
#
# THREE THINGS THIS SCRIPT GETS RIGHT THAT A NAIVE VERSION DOES NOT
#   1. The phase bands are anchored to heater switch-on at 19:30 on 11 March,
#      established from the paired 15-minute records printed in the report, not
#      to a round date. The script refuses to save a figure whose shading does
#      not match the thermal record.
#   2. Night-time means are plotted at their own mean timestamp (night_means()
#      returns `At`), not at midday. See the note in 00_common.R.
#   3. Dose zero sits at the panel floor, below both tanks' starting
#      temperatures, because the right-hand axis is anchored to the observed
#      temperature range. That is axis arithmetic, not a measurement.
# =============================================================================

source("00_common.R")

logger  <- load_logger()
paired  <- cumulative_dose(logger)
dose    <- dose_by_night(paired)
nightly <- night_means(logger)

# ---------------------------------------------------------------- Table S1
tabS1 <- nightly |>
  select(Date, Tank, n, Temp) |>
  pivot_wider(names_from = Tank, values_from = c(n, Temp)) |>
  left_join(dose, by = "Date") |>
  transmute(Date,
            n         = n_Control,
            Control   = round(Temp_Control, 2),
            Treatment = round(Temp_Heated, 2),
            Elevation = round(Temp_Heated - Temp_Control, 2),
            Cumulative_Cdays = round(Dose, 2)) |>
  arrange(Date)
write_csv(tabS1, file.path(DIR_TAB, "TableS1_nightly_temperature.csv"))
write_csv(dose |> mutate(Dose = round(Dose, 3)),
          file.path(DIR_TAB, "dose_by_night.csv"))

# ---------------------------------------------------------------- report
sink(file.path(DIR_REP, "logger_report.txt"), split = TRUE)
cat("=====================================================================\n")
cat(" LOGGER REPORT  ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("=====================================================================\n\n")
cat("Onset HOBO UA-002-64 Pendant, 15-minute interval, accurate to +/-0.2 C.\n")
cat("Serial numbers 21086852 (control tank) and 21086851 (heated tank).\n")
cat("Deployed 16:00 on 11 March; readings before 16:45 were taken in air and\n")
cat("are excluded. No record exists for 10 March or the earlier part of 11 March.\n\n")

cat("--- Night-time daily means (illuminance zero) ---\n")
print(as.data.frame(tabS1))

cat("\n--- Heater switch-on, evening of 11 March ---\n")
win <- paired |> filter(DateTime >= gt("2025-03-11 19:00:00"),
                        DateTime <= gt("2025-03-11 21:15:00")) |>
  transmute(Time = format(DateTime, "%H:%M"),
            Control = round(Control, 2), Heated = round(Heated, 2),
            Difference = round(dT, 2))
print(as.data.frame(win))
cat("\nThe two tanks were at the same temperature at 19:15 and diverged\n")
cat("immediately afterwards, placing switch-on at approximately 19:30.\n")
cat("That night's PAM session ran 19:04-19:17 and so precedes the heater.\n")
cat("The shaded thermal-stress band in Figure 1 begins at that timestamp.\n")

cat("\n--- Composition of the 11 March night-time mean ---\n")
n11 <- logger |> filter(!is.na(Light), Light == 0,
                        guam_date(DateTime) == as.Date("2025-03-11"))
print(as.data.frame(
  n11 |> mutate(Period = ifelse(DateTime < HEATER_ON, "before heater",
                                "after heater")) |>
    group_by(Tank, Period) |>
    summarise(n = n(), Temp = round(mean(Temp), 2), .groups = "drop")))
cat(sprintf("Markers for this night are drawn at %s, the mean timestamp of its\n",
            format(nightly$At[nightly$Date == as.Date("2025-03-11")][1], "%H:%M")))
cat("readings, which falls after the heater came on. The night pools a few\n")
cat("pre-heater readings, when both tanks were identical, with the rest after\n")
cat("switch-on, so its mean is not a baseline value and is not treated as one.\n")

cat("\n--- Cumulative excess heat at 19:00 each night (C-days) ---\n")
print(as.data.frame(dose |> mutate(Dose = round(Dose, 2))))
cat(sprintf("\nTotal over the whole record: %.2f C-days\n", max(paired$CumDose)))
cat("Defined relative to the concurrent unheated control, NOT to a regional\n")
cat("climatological maximum monthly mean. This is not Degree Heating Weeks\n")
cat("and must not be compared with published DHW values.\n")

cat("\n--- Daily range and daytime maxima ---\n")
rng <- logger |> mutate(Date = guam_date(DateTime)) |>
  filter(Date >= as.Date("2025-03-12"), Date <= as.Date("2025-03-16")) |>
  group_by(Tank, Date) |>
  summarise(Min = min(Temp), Max = max(Temp), Range = Max - Min, .groups = "drop")
print(as.data.frame(rng |> mutate(across(where(is.numeric), \(x) round(x, 2)))))
sink()

# ---------------------------------------------------------------- verify
cat("\n--- VERIFICATION against the published manuscript ---\n")
w <- tabS1 |> filter(Date >= as.Date("2025-03-12"), Date <= as.Date("2025-03-16"))
check("night-time means, heated tank (Table S1)",  w$Treatment, PUB$night_heat, 0.02)
check("night-time means, control tank (Table S1)", w$Control,   PUB$night_ctrl, 0.02)
check("cumulative excess heat, warming nights",
      (dose |> filter(Date >= as.Date("2025-03-12"),
                      Date <= as.Date("2025-03-16")))$Dose, PUB$dose_warming, 0.05)

# The shading self-test. A figure whose bands sit in the wrong place is worse
# than no figure, so this stops the script rather than warning.
base_dT <- mean(paired$dT[paired$DateTime <  HEATER_ON])
warm_dT <- mean(paired$dT[paired$DateTime >= HEATER_ON & paired$DateTime < WARM_END])
post_dT <- mean(paired$dT[paired$DateTime >= RECOV_END])
cat(sprintf("  mean difference before the band starts : %+.2f C\n", base_dT))
cat(sprintf("  mean difference inside the warming band: %+.2f C\n", warm_dT))
cat(sprintf("  mean difference after the experiment   : %+.2f C\n", post_dT))
if (abs(base_dT) > 0.35)
  stop("Shading error: the tanks already differ before the band starts.")
if (warm_dT < 1.5)
  stop("Shading error: the warming band does not cover the heated period.")
if (abs(post_dT) > 0.35)
  stop("Unexpected difference after the experiment; check the record.")
cat("  [ok  ] shading matches the thermal record\n")
cat("  all logger checks passed\n\n")

# ------------------------------------------------- right-hand axis transform
# Anchored to the observed temperature range, not adjusted by eye. Dose zero
# therefore lands on the panel floor, below both tanks.
t_lo <- floor(min(logger$Temp)); t_hi <- ceiling(max(logger$Temp))
d_max <- max(paired$CumDose)
B <- t_lo; K <- (t_hi - t_lo) / d_max
to_primary <- function(d) d * K + B
to_dose    <- function(y) (y - B) / K
cat(sprintf("Secondary axis: 0 C-days -> %.1f C ; %.2f C-days -> %.1f C\n\n",
            B, d_max, to_primary(d_max)))

bands <- tibble(
  Phase = factor(names(phase_cols), levels = names(phase_cols)),
  xmin  = c(PANEL_START, HEATER_ON, WARM_END),
  xmax  = c(HEATER_ON,   WARM_END,  RECOV_END))

# ------------------------------------------------------------ figure builder
# `last` is the right-hand edge; everything else is identical between Figure 1
# and Figure S3, so neither can drift away from the other.
thermal_figure <- function(last, brk = "2 days"){
  panel_end <- min(max(paired$DateTime), gt(paste(last, "00:00:00")))
  lg  <- filter(logger,  DateTime <= panel_end)
  ng  <- filter(nightly, At       <= panel_end)
  d_h <- filter(paired,  DateTime <= RECOV_END)
  d_p <- filter(paired,  DateTime >= RECOV_END, DateTime <= panel_end)

  ggplot() +
    geom_rect(data = bands,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Phase),
              inherit.aes = FALSE, alpha = 0.8) +
    scale_fill_manual(values = phase_cols, name = NULL) +
    # the interval for which no thermal record exists
    annotate("segment", x = PANEL_START, xend = SUBMERGED,
             y = t_lo + 0.35, yend = t_lo + 0.35, colour = "grey55",
             linewidth = 0.4,
             arrow = arrow(ends = "both", length = unit(0.05, "in"))) +
    annotate("text", x = PANEL_START + as.difftime(20, units = "hours"),
             y = t_lo + 0.35, label = "no logger record", vjust = -0.7,
             size = 2.8, colour = "grey40", fontface = "italic") +
    geom_vline(xintercept = HEATER_ON, linetype = "dashed",
               colour = "grey25", linewidth = 0.45) +
    annotate("text", x = HEATER_ON + as.difftime(3, units = "hours"), y = Inf,
             label = "heater on, 19:30 on 11 Mar", hjust = 0, vjust = 1.6,
             size = 2.9, colour = "grey25") +
    # raw 15-minute records
    geom_line(data = lg, aes(DateTime, Temp, colour = Tank),
              linewidth = 0.27, alpha = 0.45) +
    # night-time daily means, at the mean timestamp of their own readings
    geom_line(data = ng, aes(At, Temp, colour = Tank), linewidth = 0.9) +
    geom_point(data = ng, aes(At, Temp, colour = Tank),
               size = 1.9, shape = 21, fill = "white", stroke = 1) +
    # cumulative excess heat, right-hand axis; grey dashed once the heater is off
    geom_line(data = d_h, aes(DateTime, to_primary(CumDose)),
              linewidth = 1.0, colour = DOSE_COL) +
    geom_line(data = d_p, aes(DateTime, to_primary(CumDose)),
              linewidth = 0.8, colour = "grey60", linetype = "dashed") +
    geom_point(data = dose |> filter(Date <= as.Date("2025-03-17")),
               aes(gt(paste(Date, "19:00:00")), to_primary(Dose)),
               size = 1.9, shape = 21, fill = "white", colour = DOSE_COL,
               stroke = 1) +
    annotate("text", x = HEATER_ON + as.difftime(4, units = "hours"),
             y = to_primary(0), label = "cumulative excess heat",
             hjust = 0, vjust = -0.9, size = 3.0, colour = DOSE_COL,
             fontface = "italic") +
    scale_colour_manual(values = tank_cols[c("Control", "Heated")], name = NULL) +
    scale_x_datetime(limits = c(PANEL_START, panel_end),
                     breaks = seq(PANEL_START, panel_end, by = brk),
                     date_labels = "%d %b", expand = expansion(mult = 0.012)) +
    scale_y_continuous(
      name = expression("Temperature ("*degree*"C)"),
      sec.axis = sec_axis(~ to_dose(.),
                   name = expression("Cumulative excess heat ("*degree*"C-days)"))) +
    labs(x = "Date (2025)") +
    guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1)) +
    theme_pub() +
    theme(axis.title.y.right = element_text(colour = DOSE_COL),
          axis.text.y.right  = element_text(colour = DOSE_COL),
          axis.ticks.y.right = element_line(colour = DOSE_COL),
          legend.box = "horizontal")
}

# ---------------------------------------------------------------- Figure 1
save_fig(thermal_figure(PANEL_LAST),
         "Figure1_thermal_history_and_dose", w = 10.2, h = 5.4)

# ------------------------------------------------- Supplementary Figure S3
# The complete record, to the end of logging on 27 March.
save_fig(thermal_figure(format(max(guam_date(paired$DateTime)) + 1), brk = "3 days"),
         "FigureS3_complete_logger_record", w = 10.2, h = 5.4)

cat("01_logger_analysis.R complete.\n")
