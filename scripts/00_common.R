# =============================================================================
# 00_common.R — shared settings, helpers and verification targets
#
# Sourced by 01, 02 and 03. Do not run on its own.
#
# TIME ZONE WARNING
#   as.Date() applied to a POSIXct silently converts to UTC, not to the
#   object's own time zone. Guam is UTC+10, so that shifts every timestamp by
#   ten hours and corrupts both the night-time means and the dose integration.
#   Use guam_date() and guam_hour() below, never as.Date(), on a POSIXct.
# =============================================================================

TZ <- "Pacific/Guam"

# ENCODING WARNING
#   The figures contain the degree sign and a true minus. In a session whose
#   locale is not UTF-8 the graphics devices cannot map them and silently draw
#   dots instead, so an axis reads "+2.5..C". Nothing errors, so check this.
if (!isTRUE(l10n_info()[["UTF-8"]]))
  warning("This R session is not using a UTF-8 locale (LC_CTYPE = ",
          Sys.getlocale("LC_CTYPE"), "). Characters such as the degree sign ",
          "will be dropped from the figures. Restart R with, for example, ",
          "LANG=en_US.UTF-8 before generating figures.", call. = FALSE)

# Repository root. Resolved automatically when the script is run from
# scripts/, otherwise edit this one line.
if (!exists("REPO")) REPO <- normalizePath("..", mustWork = FALSE)
if (!dir.exists(file.path(REPO, "data")))
  stop("REPO does not point at the repository root: ", REPO,
       "\nSet REPO <- \"/path/to/guam-coral-thermal-stress\" before sourcing.")

DIR_RAW   <- file.path(REPO, "data", "raw")
DIR_PROC  <- file.path(REPO, "data", "processed")
DIR_TAB   <- file.path(REPO, "results", "tables")
DIR_REP   <- file.path(REPO, "results", "reports")
DIR_FIG   <- file.path(REPO, "results", "figures")
for (d in c(DIR_TAB, DIR_REP, DIR_FIG))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("readr","dplyr","tidyr","ggplot2","ggbeeswarm","purrr","scales")
invisible(lapply(pkgs, function(p){
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))
# guard against plyr / Hmisc masking the dplyr verbs
filter <- dplyr::filter; select <- dplyr::select; mutate <- dplyr::mutate
summarise <- dplyr::summarise; arrange <- dplyr::arrange; count <- dplyr::count

# ---- time-zone safe date handling -------------------------------------------
guam_date <- function(x) as.Date(format(x, "%Y-%m-%d", tz = TZ))
guam_hour <- function(x) as.numeric(format(x, "%H", tz = TZ)) +
                         as.numeric(format(x, "%M", tz = TZ)) / 60
gt        <- function(s) as.POSIXct(s, tz = TZ)

# ---- experiment constants ---------------------------------------------------
SUBMERGED <- gt("2025-03-11 16:45:00")  # loggers were in air before this
HEATER_ON <- gt("2025-03-11 19:30:00")  # established from the paired record
EXP_START <- gt("2025-03-11 00:00:00")
EXP_END   <- gt("2025-03-17 23:59:59")

# Phase boundaries used for the shaded bands in Figure 1. The baseline phase
# of the experiment is 10-11 March, two nights of PAM measurements, and begins
# before the loggers were deployed; the interval with no thermal record is
# marked as such on the figure rather than left blank.
PANEL_START <- gt("2025-03-10 00:00:00")  # start of the baseline phase
WARM_END    <- gt("2025-03-17 00:00:00")  # start of the recovery night
RECOV_END   <- gt("2025-03-18 00:00:00")  # tanks converged after this

# Last date shown in Figure 1. The complete record runs to 27 March and is
# plotted in full as Supplementary Figure S3.
PANEL_LAST  <- "2025-03-19"

DOSE_COL <- "#7F0000"                     # cumulative excess heat curve

# Raster resolution for the PNG copies. 300 is fine for drafts; set 900 before
# a submission that asks for it:  FIG_DPI <- 900  (or Sys.setenv(FIG_DPI=900))
if (!exists("FIG_DPI"))
  FIG_DPI <- as.numeric(Sys.getenv("FIG_DPI", unset = "300"))

BASELINE_NIGHTS <- "10and11"   # "10and11" = primary; "10only" = sensitivity
BASE_DATES <- if (BASELINE_NIGHTS == "10only") as.Date("2025-03-10") else
  as.Date(c("2025-03-10", "2025-03-11"))

sp_levels <- c("Acropora","Pavona","Sclerophytum")
sp_labs   <- c(Acropora = "Acropora cf. manni", Pavona = "Pavona decussata",
               Sclerophytum = "Sclerophytum sp.")

# Plotmath versions of the same labels. A species name is italic but the
# connecting terms "cf." and "sp." are not, which plain italic text cannot
# express. Use these wherever a species name is drawn:
#   facet_wrap(..., labeller = labeller(Species = sp_labeller))
#   scale_x_discrete(labels = sp_axis_labels())
# and do NOT also set face = "italic" on that axis, or the whole label,
# "cf." included, is slanted again.
sp_expr <- c(
  Acropora     = 'italic("Acropora")~"cf."~italic("manni")',
  Pavona       = 'italic("Pavona decussata")',
  Sclerophytum = 'italic("Sclerophytum")~"sp."')
# A plain function returning expressions is NOT enough for a facet strip:
# ggplot2's labeller machinery coerces the result back to character, and the
# strip then shows the plotmath source instead of rendering it. as_labeller()
# with default = label_parsed is the form that actually parses.
sp_labeller    <- as_labeller(sp_expr, default = label_parsed)
sp_axis_labels <- function(levels = sp_levels) parse(text = sp_expr[levels])

# Two-line versions, for axes where the names would otherwise collide
sp_expr_2line <- c(
  Acropora     = 'atop(italic("Acropora"), "cf."~italic("manni"))',
  Pavona       = 'atop(italic("Pavona"), italic("decussata"))',
  Sclerophytum = 'atop(italic("Sclerophytum"), "sp.")')
sp_axis_labels_2line <- function(levels = sp_levels)
  parse(text = sp_expr_2line[levels])

# Axis titles that need subscripts, matching the manuscript's typography.
LAB_FVFM     <- expression(F[v]/F[m]~"(% of baseline)")
LAB_FVFM_RAW <- expression(F[v]/F[m])
sp_cols   <- c(Acropora = "#0072B2", Pavona = "#D55E00", Sclerophytum = "#009E73")
# "Heat" is the PAM treatment label, "Heated" the logger tank label; both are
# the same tank and take the same colour.
tank_cols <- c(Control = "#4393C3", Heat = "#D6604D", Heated = "#D6604D")
phase_cols <- c(Baseline = "#EDEDED", `Thermal Stress` = "#FBE3D6",
                Recovery = "#E4EEF7")

# Measured night-time elevation of the heated tank above control (Table S1)
ELEV_MAP <- c("2025-03-10" = NA, "2025-03-11" = 0.0, "2025-03-12" = 2.0,
              "2025-03-13" = 2.1, "2025-03-14" = 2.5, "2025-03-15" = 2.9,
              "2025-03-16" = 3.5, "2025-03-17" = 1.9)

# ---- values published in the manuscript, used to verify each script ---------
# If a script stops here, the output no longer matches the paper. Do not
# publish figures or tables from a run that failed verification.
PUB <- list(
  dose_warming   = c(2.17, 4.45, 6.99, 9.98, 13.87),   # 19:00, 12-16 March
  night_heat     = c(30.77, 30.73, 31.12, 31.62, 32.30),
  night_ctrl     = c(28.74, 28.63, 28.64, 28.69, 28.83),
  baseline_F     = 59.27,                               # F(2,83), baseline ANOVA
  warming_F      = 7.82,                                # F(2,137), between species
  warming_KW     = 13.33,
  dose_slope_pav = -2.45, dose_p_pav = 0.006,
  qpcr_endo_F    = 8.57,  qpcr_duru_F = 10.94, qpcr_clad_F = 0.62
)

check <- function(label, got, want, tol, fatal = TRUE){
  ok <- all(abs(got - want) <= tol)
  cat(sprintf("  [%s] %-44s got %s | published %s\n",
              if (ok) "ok  " else "FAIL", label,
              paste(round(got, 3), collapse = ", "),
              paste(round(want, 3), collapse = ", ")))
  if (!ok && fatal)
    stop("Verification failed: ", label,
         "\nThe output no longer reproduces the manuscript. Do not use it.")
  invisible(ok)
}

# ---- figure saving: PNG + SVG + PDF -----------------------------------------
# SVG DEVICE, AND WHY IT MATTERS
#   svglite writes real <text> elements, so every label stays editable in
#   Inkscape or Illustrator. grDevices::svg() is cairo-backed and converts
#   each character to an outline path: the figure looks identical but the
#   text becomes groups of letter-shaped paths that cannot be retyped.
#   svglite is used when it is installed, and the fallback says what was lost.
#   If svglite fails to load over a systemfonts version clash, reinstall both
#   from source together:
#     install.packages(c("systemfonts","svglite"), type = "source")
save_fig <- function(p, name, w, h, dpi = FIG_DPI){
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white")
  svg_path <- file.path(DIR_FIG, paste0(name, ".svg"))
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(svg_path, p, width = w, height = h,
           device = svglite::svglite, bg = "white")
  } else {
    message("svglite is not available: writing ", name,
            ".svg with grDevices::svg(), which converts text to outlines. ",
            "Install svglite if you need to edit the labels afterwards.")
    svg(svg_path, width = w, height = h, bg = "white"); print(p); dev.off()
  }
  ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  message("saved: ", name, " (.png .svg .pdf)")
}

theme_pub <- function(base = 11){
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.5),
          axis.text = element_text(colour = "grey15"),
          strip.background = element_rect(fill = "grey93", colour = "grey55",
                                          linewidth = 0.5),
          strip.text = element_text(face = "bold.italic", size = base),
          legend.position = "bottom", legend.title = element_blank(),
          legend.key = element_blank())
}

# ---- logger reading and cumulative excess heat ------------------------------
read_hobo <- function(path, label){
  raw  <- read_csv(path, skip = 1, show_col_types = FALSE, name_repair = "minimal")
  names(raw) <- make.unique(names(raw))
  tcol <- grep("^Temp",      names(raw))[1]
  icol <- grep("^Intensity", names(raw))[1]
  dcol <- grep("Date Time",  names(raw))[1]
  tibble(DateTime = as.POSIXct(raw[[dcol]], format = "%m/%d/%y %I:%M:%S %p", tz = TZ),
         Temp     = suppressWarnings(as.numeric(raw[[tcol]])),
         Light    = suppressWarnings(as.numeric(raw[[icol]])),
         Tank     = label) |>
    filter(!is.na(DateTime), !is.na(Temp))
}

load_logger <- function(){
  f_c <- file.path(DIR_RAW, "logger", "VK_Control_C.csv")
  f_h <- file.path(DIR_RAW, "logger", "VK_Treatment_T.csv")
  miss <- c(f_c, f_h)[!file.exists(c(f_c, f_h))]
  if (length(miss)) stop("Missing logger file(s):\n  ", paste(miss, collapse = "\n  "))
  bind_rows(read_hobo(f_c, "Control"), read_hobo(f_h, "Heated")) |>
    filter(DateTime >= SUBMERGED)
}

# Cumulative excess heat, in degree C-days, integrated over the paired
# 15-minute records:  sum( (T_treatment - T_control) * 0.25/24 ).
# Defined against the concurrent control, NOT a climatological maximum
# monthly mean, so it is not Degree Heating Weeks and must not be compared
# with published DHW values.
cumulative_dose <- function(logger){
  logger |> select(DateTime, Tank, Temp) |>
    pivot_wider(names_from = Tank, values_from = Temp) |>
    filter(!is.na(Control), !is.na(Heated)) |> arrange(DateTime) |>
    mutate(dT = Heated - Control, CumDose = cumsum(dT * (0.25 / 24)))
}

# Dose accumulated by 19:00 on each night, the reference point used throughout
# the manuscript (PAM sessions ran 18:52-19:32).
dose_by_night <- function(paired){
  paired |>
    mutate(Date = guam_date(DateTime), hr = guam_hour(DateTime)) |>
    filter(hr <= 19) |> group_by(Date) |>
    summarise(Dose = dplyr::last(CumDose), .groups = "drop")
}

# Night-time daily means: readings taken when logger illuminance was zero.
# Solar heating cannot contribute to these, so they isolate the heater signal.
#
# `At` is the mean timestamp of the readings that produced each mean, and is
# what Figure 1 plots the markers at. Plotting them at midday instead puts the
# 11 March point before the loggers were in the water and to the LEFT of the
# heater switch-on line, making the heated tank look warm before heating began:
# that night's readings run 18:30-23:45 and average 21:07. Full nights land at
# 11:29, which is indistinguishable from midday.
night_means <- function(logger){
  logger |> filter(!is.na(Light), Light == 0) |>
    mutate(Date = guam_date(DateTime)) |>
    group_by(Date, Tank) |>
    summarise(n = n(), Temp = mean(Temp),
              At = as.POSIXct(mean(as.numeric(DateTime)),
                              origin = "1970-01-01", tz = TZ),
              .groups = "drop")
}
