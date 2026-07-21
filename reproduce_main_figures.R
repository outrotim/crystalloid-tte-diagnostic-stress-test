#!/usr/bin/env Rscript

# Reproduce the three main figures from non-disclosive aggregate inputs only.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "reproduce_main_figures.R"
root <- dirname(normalizePath(script_path, mustWork = TRUE))
input_path <- file.path(root, "aggregate_figure_data.csv")
output_dir <- file.path(root, "reproduced_figures")
dir.create(output_dir, showWarnings = FALSE)

x <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("record_type", "exposure_definition", "population_variant",
              "risk_difference", "rd_lower", "rd_upper", "study",
              "estimate", "lower", "upper")
if (!all(required %in% names(x))) stop("Aggregate input schema is incomplete")

theme_publication <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# Figure 1: reusable identification-diagnostic cascade.
cascade <- subset(x, record_type == "diagnostic_cascade")
cascade$diagnostic_order <- as.numeric(cascade$diagnostic_order)
cascade <- cascade[order(cascade$diagnostic_order), ]
cascade$label <- paste0(cascade$diagnostic_order, ". ", cascade$diagnostic_step,
                        "\n", cascade$empirical_signal)
p1 <- ggplot(cascade, aes(x = diagnostic_order, y = 1, label = label)) +
  geom_tile(width = 0.88, height = 0.62, fill = "#F5F5F2", color = "#6B7280") +
  geom_text(size = 2.65, lineheight = 0.92) +
  geom_segment(data = head(cascade, -1), aes(xend = diagnostic_order + 0.54,
               yend = 1), arrow = arrow(length = grid::unit(0.08, "inches"))) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  labs(title = "A reusable diagnostic cascade for EHR target trial emulation") +
  theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold"))

# Figure 2: post-baseline volume eligibility and exposure operationalization.
rd <- subset(x, record_type == "risk_difference")
for (column in c("risk_difference", "rd_lower", "rd_upper")) rd[[column]] <- 100 * as.numeric(rd[[column]])
panel_a <- subset(rd, population_variant == "full" & exposure_definition == "initial")
panel_a$threshold <- factor(panel_a$volume_threshold_ml,
  levels = c(1000, 500, 0), labels = c("≥1,000 mL", "≥500 mL", "No minimum"))
p2a <- ggplot(panel_a, aes(risk_difference, threshold)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#666666") +
  geom_errorbar(aes(xmin = rd_lower, xmax = rd_upper), orientation = "y", width = 0) +
  geom_point(size = 2.7, color = "#0072B2") +
  labs(title = "A  Post-time-zero volume eligibility", y = NULL,
       x = "MAKE-30 risk difference (percentage points)") + theme_publication

panel_b <- subset(rd, as.numeric(volume_threshold_ml) == 1000)
panel_b$definition <- factor(panel_b$exposure_definition,
  levels = c("initial", "dominant_70", "dominant_80", "dominant_90"),
  labels = c("Initial assignment", "48-h dominant ≥70%", "48-h dominant ≥80%", "48-h dominant ≥90%"))
panel_b$population <- factor(panel_b$population_variant,
  levels = c("full", "same_cohort", "landmark"),
  labels = c("Full eligible cohort", "Common classifiable cohort", "48-h landmark cohort"))
p2b <- ggplot(panel_b, aes(risk_difference, definition, color = population)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#666666") +
  geom_errorbar(aes(xmin = rd_lower, xmax = rd_upper), orientation = "y",
                width = 0, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.2) +
  labs(title = "B  Exposure operationalization at ≥1,000 mL", y = NULL,
       x = "MAKE-30 risk difference (percentage points)", color = NULL) + theme_publication
p2 <- p2a / p2b + plot_layout(heights = c(0.75, 1.25), guides = "collect")

# Figure 3: randomized-evidence calibration; effect scales are displayed, not pooled.
cal <- subset(x, record_type == "rct_calibration")
cal$study <- factor(cal$study, levels = rev(unique(cal$study)))
p3 <- ggplot(cal, aes(estimate, study, color = source_type, shape = source_type)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#666666") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.16) +
  geom_point(size = 2.5) + facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_x_log10() +
  labs(title = "Randomized-evidence calibration of observational estimates",
       subtitle = "Outcome-specific intervals are displayed, not pooled",
       x = "Effect ratio (balanced crystalloid vs saline; log scale)", y = NULL,
       color = NULL, shape = NULL) + theme_publication

figures <- list(Figure_1 = p1, Figure_2 = p2, Figure_3 = p3)
sizes <- list(Figure_1 = c(11.0, 3.0), Figure_2 = c(7.2, 6.9), Figure_3 = c(7.2, 4.6))
for (name in names(figures)) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")), figures[[name]],
         width = sizes[[name]][1], height = sizes[[name]][2], device = cairo_pdf)
  ggsave(file.path(output_dir, paste0(name, ".png")), figures[[name]],
         width = sizes[[name]][1], height = sizes[[name]][2], dpi = 300, bg = "white")
}

message("Created six figure files in: ", output_dir)
