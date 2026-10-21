# KRTINY: FUNGAL COMMUNITY CHANGE, 2013-2023
# =========================================
# A worked example using ITS2 operational taxonomic units (OTUs).
# Based on Krtiny_ITS2_temporal_turnover_2013_2023_upd_TM.Rmd
# by Tijana Martinovic. Data source: Krtiny_ITS2_2013_to_2023.xlsx.
#
# Questions
#   1. Do samples from different sampling rounds separate in an NMDS plot?
#   2. Is community composition associated with sampling year?
#   3. Are some years more variable than others?
#
# START HERE
# Open this script in RStudio. Choose Session > Set Working Directory >
# To Source File Location. Run one section at a time with Ctrl+Enter,
# or use Source to run everything. No original files are needed.
#
# The companion workbook has exactly three sheets, with headers in row 1:
#   OTU_table       : one OTU per row, one sample per column; read counts.
#   Sample_metadata : one sample per row; Sample_ID links to OTU_table.
#   Taxonomy        : one OTU per row; OTU_ID links to OTU_table.
#
# The one-time preparation of this workbook is documented at the end.
# Data cleaning has already been done. Statistical choices remain visible here.


# 1. Packages and settings -------------------------------------------------

# Install these packages once, if needed (remove the first # to run):
# install.packages(c("readxl", "vegan", "ggplot2"))
# vegan also installs permute, used below to restrict permutations by site.
library(readxl)
library(vegan)
library(ggplot2)

# Change this to "Saprotroph", "EcM", or "Yeast" to repeat the analysis.
# EcM means ectomycorrhizal fungi. "All" includes every retained OTU.
selected_guild <- "All"
substrates <- c("Soil", "Litter")
nmds_dimensions <- 3       # The original analysis used three dimensions.
n_permutations <- 999
save_outputs <- FALSE     # TRUE creates a results folder beside the workbook.
run_turnover_extension <- FALSE  # Optional descriptive time-decay plot.

input_file <- "Krtiny_ITS2_student_data.xlsx"
# This alternative also lets you run the script from the parent folder.
if (!file.exists(input_file)) {
  input_file <- file.path("example_for_students", input_file)
}
if (!file.exists(input_file)) {
  stop("Set the working directory to the folder containing the student workbook.")
}
stopifnot(selected_guild %in% c("All", "Saprotroph", "EcM", "Yeast"),
          nmds_dimensions %in% c(2, 3))

sampling_rounds <- c("2013", "2016", "2017a", "2017b", "2020", "2023")
year_colours <- c("2013" = "#66513D", "2016" = "#A67C52",
                 "2017a" = "#74805A", "2017b" = "#A4AF8B",
                 "2020" = "#C17D58", "2023" = "#466C68")


# 2. Read the three sheets -------------------------------------------------

otu_table <- as.data.frame(read_excel(input_file, sheet = "OTU_table"))
metadata <- as.data.frame(read_excel(input_file, sheet = "Sample_metadata"))
taxonomy <- as.data.frame(read_excel(input_file, sheet = "Taxonomy"))

# Column 1 contains OTU IDs. The other columns contain numeric counts.
stopifnot(!anyNA(otu_table$OTU_ID), !anyDuplicated(otu_table$OTU_ID),
          !anyNA(metadata$Sample_ID), !anyDuplicated(metadata$Sample_ID),
          !anyNA(taxonomy$OTU_ID), !anyDuplicated(taxonomy$OTU_ID),
          all(vapply(otu_table[-1], is.numeric, logical(1))))
counts <- as.matrix(otu_table[-1])
rownames(counts) <- otu_table$OTU_ID
counts <- t(counts)        # vegan needs SAMPLES in rows and OTUs in columns.

# Match by IDs, even if someone has sorted a sheet differently in Excel.
stopifnot(setequal(rownames(counts), metadata$Sample_ID),
          setequal(colnames(counts), taxonomy$OTU_ID))
metadata <- metadata[match(rownames(counts), metadata$Sample_ID), ]
taxonomy <- taxonomy[match(colnames(counts), taxonomy$OTU_ID), ]
stopifnot(identical(rownames(counts), metadata$Sample_ID),
          identical(colnames(counts), taxonomy$OTU_ID),
          !anyNA(counts), all(is.finite(counts)), all(counts >= 0),
          all(counts == floor(counts)), all(rowSums(counts) > 0),
          all(rowSums(counts) == metadata$Reads_in_OTU_table),
          !anyNA(metadata), all(metadata$Substrate %in% c("Soil", "Litter")),
          all(metadata$Sampling_round %in% sampling_rounds),
          all(metadata$Year %in% c(2013, 2016, 2017, 2020, 2023)),
          all(as.matrix(taxonomy[c("Saprotroph", "EcM", "Yeast")]) %in% 0:1))

# Factors are categorical variables; Year stays numeric for the trend test.
metadata$Site <- factor(metadata$Site)
metadata$Sampling_round <- factor(metadata$Sampling_round,
                                 levels = sampling_rounds)
metadata$Year_group <- factor(metadata$Year)
# Year_group combines the two 2017 rounds for PERMANOVA, as in the source Rmd.
# Both 2017 observations are retained; counts are not pooled or averaged.

cat("Samples:", nrow(counts), "  OTUs:", ncol(counts), "\n")
print(table(metadata$Substrate, metadata$Sampling_round))
print(colSums(taxonomy[c("Saprotroph", "EcM", "Yeast")]))


# 3. Select a guild and analyse each substrate -----------------------------

# Each guild flag in Taxonomy is 1 (assigned) or 0 (not assigned).
# Missing taxonomic names do not remove an OTU from the all-fungi analysis.
keep_otus <- rep(TRUE, ncol(counts))
if (selected_guild != "All") {
  keep_otus <- taxonomy[[selected_guild]] == 1
}

results <- list()          # Keeps soil and litter results separately in R.
output_dir <- file.path(dirname(input_file), "results", selected_guild)
if (save_outputs) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# A for-loop repeats the same steps: first Soil, then Litter.
# In RStudio, run this whole loop together (from 'for' to its closing '}').
for (substrate in substrates) {
  cat("\n---", substrate, ":", selected_guild, "---\n")
  meta <- droplevels(metadata[metadata$Substrate == substrate, ])
  community <- counts[meta$Sample_ID, keep_otus, drop = FALSE]

  # A guild can be absent from a sample. Distances for an empty community
  # are not meaningful, so report and remove such samples for this guild only.
  empty <- rowSums(community) == 0
  cat("Samples with no detected OTUs in this selection:", sum(empty), "\n")
  print(meta$Sample_ID[empty])
  community <- community[!empty, , drop = FALSE]
  meta <- droplevels(meta[!empty, ])
  community <- community[, colSums(community) > 0, drop = FALSE]
  stopifnot(identical(rownames(community), meta$Sample_ID),
            nrow(community) > 2 * nmds_dimensions + 1,
            ncol(community) > 1)
  # Guild-specific removal can make the repeated series incomplete.
  # Keep the remaining observations; do not silently drop an entire site.
  print(table(meta$Year_group))


  # 4. Presence/absence and Sorensen dissimilarity ------------------------

  # We ask WHICH OTUs occur, not how many reads each OTU received.
  presence_absence <- 1L * (community > 0)
  dissimilarity <- vegdist(presence_absence, method = "bray", binary = TRUE)
  # Bray-Curtis on binary data equals Sorensen dissimilarity:
  #   1 - 2 * shared_OTUs / (OTUs_in_sample_1 + OTUs_in_sample_2).
  # 0 = identical detected OTUs; 1 = no shared detected OTUs.
  stopifnot(all(is.finite(dissimilarity)))
  cat("Samples:", nrow(community), "  OTUs:", ncol(community), "\n")
  cat("Fraction of sample pairs sharing no OTUs:",
      round(mean(dissimilarity == 1), 3), "\n")


  # 5. Fit NMDS and inspect its quality -----------------------------------

  # NMDS tries to preserve the RANK ORDER of community dissimilarities.
  # Multiple random starts help avoid a poor local solution.
  # We pass the distance object, so no hidden count transformation is used.
  set.seed(4)
  nmds <- metaMDS(dissimilarity, k = nmds_dimensions,
                  trymax = 200, autotransform = FALSE, trace = FALSE)
  print(nmds)
  if (nmds$converged == 0) {
    warning("The best NMDS solution was not repeated. Inspect convergence and try more starts.")
  }
  if (nmds$stress > 0.2) {
    warning("NMDS stress is above 0.2. Interpret the plot cautiously and inspect the fit.")
  }
  # Lower stress means a better representation. 0.2 is a rough warning,
  # not a universal pass/fail rule. Also inspect the Shepard (stress) plot.
  # Stress refers to ALL fitted dimensions, not just the two plotted axes.
  stressplot(nmds, main = paste(substrate, selected_guild, "- NMDS fit"))

  # Extract all fitted coordinates and attach the matching metadata.
  coordinates <- as.data.frame(scores(nmds, display = "sites",
                                      choices = seq_len(nmds_dimensions)))
  plot_data <- cbind(meta[match(rownames(coordinates), meta$Sample_ID), ],
                     coordinates)
  nmds_plot <- ggplot(plot_data, aes(NMDS1, NMDS2, colour = Sampling_round)) +
    geom_point(size = 2.7, alpha = 0.85) +
    scale_colour_manual(values = year_colours, drop = FALSE) +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "right") +
    labs(title = paste(substrate, "-", selected_guild),
         subtitle = sprintf("%dD NMDS; stress = %.3f; %d samples",
                            nmds_dimensions, nmds$stress, nrow(meta)),
         x = "NMDS1", y = "NMDS2", colour = "Sampling round",
         caption = "Axes 1 and 2 shown. Inspect NMDS3 as well when using a 3D solution.")
  print(nmds_plot)
  # NMDS axes have no units or percentages of explained variance.
  # Soil and litter have independent ordinations: their coordinates cannot
  # be compared directly. Separation in a plot is not a significance test.
  # In a 3D solution, points that overlap here may separate along NMDS3.
  # Do not remove an extreme point just to make the plot look tidier.


  # 6. PERMANOVA: is composition associated with year? --------------------

  # The same sites were sampled repeatedly. Shuffle observations WITHIN
  # sites, never between sites. This follows the original Rmd design.
  # It assumes exchangeability within a site under the null hypothesis;
  # blocking alone does not account for every form of temporal dependence.
  permutation_design <- permute::how(blocks = meta$Site,
                                     nperm = n_permutations)
  set.seed(4)
  year_test <- adonis2(dissimilarity ~ Year_group, data = meta,
                       permutations = permutation_design)
  set.seed(4)
  trend_test <- adonis2(dissimilarity ~ Year, data = meta,
                        permutations = permutation_design)
  cat("\nCategorical year (2017a and 2017b combined):\n")
  print(year_test)
  cat("\nNumeric year (linear temporal association):\n")
  print(trend_test)
  # Read R2 (variation attributed to year), pseudo-F and Pr(>F) (p-value).
  # A categorical effect does not identify which pairs of years differ.
  # A numeric-year effect does not prove a monotonic decline in similarity.
  # These are two related questions, not independent confirmations.
  # Trying many guilds/tests also requires attention to multiple testing.


  # 7. Check within-year dispersion --------------------------------------

  # PERMANOVA can respond to different group spreads as well as locations.
  # betadisper measures distances to each year's spatial median.
  # Bias adjustment helps with unequal group sizes (2017 has two rounds).
  dispersion <- betadisper(dissimilarity, group = meta$Year_group,
                           type = "median", bias.adjust = TRUE)
  set.seed(4)
  dispersion_test <- permutest(dispersion, permutations = permutation_design)
  print(dispersion_test)
  boxplot(dispersion, main = paste(substrate, selected_guild, "- dispersion"),
          xlab = "Year", ylab = "Distance to the year spatial median")
  # If dispersion differs, interpret PERMANOVA with additional care.
  # A non-significant dispersion test does not prove equal dispersion.


  # 8. Optional extension: similarity versus elapsed years ---------------

  turnover <- NULL
  if (run_turnover_extension) {
    # Use each unordered sample pair once, and compare only the SAME site.
    pairs <- which(upper.tri(as.matrix(dissimilarity)), arr.ind = TRUE)
    pairs <- pairs[meta$Site[pairs[, 1]] == meta$Site[pairs[, 2]], , drop = FALSE]
    turnover <- data.frame(
      Site = meta$Site[pairs[, 1]],
      Sample_1 = meta$Sample_ID[pairs[, 1]],
      Sample_2 = meta$Sample_ID[pairs[, 2]],
      Elapsed_years = abs(meta$Year[pairs[, 1]] - meta$Year[pairs[, 2]]),
      Similarity = 1 - as.matrix(dissimilarity)[pairs]
    )
    turnover_plot <- ggplot(turnover, aes(Elapsed_years, Similarity)) +
      geom_point(alpha = 0.25, position = position_jitter(width = 0.08,
                                                        height = 0, seed = 4)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw(base_size = 12) +
      labs(title = paste(substrate, "-", selected_guild, "temporal turnover"),
           x = "Elapsed time (years)", y = "Sorensen similarity")
    print(turnover_plot)
    # 2017a versus 2017b has a zero CALENDAR-year difference, not necessarily
    # zero elapsed time. Dates within 2017 are not available in this workbook.
    # Pairs share samples and sites, so they are NOT independent replicates.
    # This plot is descriptive. The original mixed-effects models and their
    # uncertainty estimates are outside this introductory NMDS example.
  }


  # 9. Keep the results and optionally save them --------------------------

  results[[substrate]] <- list(
    metadata = meta, presence_absence = presence_absence,
    dissimilarity = dissimilarity, nmds = nmds, coordinates = plot_data,
    nmds_plot = nmds_plot, year_test = year_test, trend_test = trend_test,
    dispersion = dispersion, dispersion_test = dispersion_test,
    turnover = turnover
  )
  if (save_outputs) {
    prefix <- file.path(output_dir, substrate)
    ggsave(paste0(prefix, "_NMDS.png"), nmds_plot,
           width = 8, height = 5.5, dpi = 300)
    write.csv(plot_data, paste0(prefix, "_NMDS_coordinates.csv"), row.names = FALSE)
    capture.output(nmds, year_test, trend_test, dispersion_test,
                   file = paste0(prefix, "_model_summary.txt"))
    pdf(paste0(prefix, "_diagnostics.pdf"), width = 7, height = 5)
    stressplot(nmds, main = paste(substrate, "- NMDS fit"))
    boxplot(dispersion, xlab = "Year", ylab = "Distance to spatial median")
    dev.off()
    if (run_turnover_extension) {
      write.csv(turnover, paste0(prefix, "_temporal_pairs.csv"), row.names = FALSE)
      ggsave(paste0(prefix, "_temporal_turnover.png"), turnover_plot,
             width = 7, height = 5, dpi = 300)
    }
  }
}

# For example, inspect results$Soil$year_test or results$Litter$nmds$stress.
if (save_outputs) {
  saveRDS(results, file.path(output_dir, "analysis_objects.rds"))
  capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
}
sessionInfo()             # Record R and package versions for reproducibility.

# Checkpoints for the supplied workbook (seed 4, vegan 2.7.6, k = 3):
# All fungi: stress is approximately 0.177 for Soil and 0.146 for Litter.
# Soil EcM (~0.201) and both yeast ordinations (~0.222-0.223) trigger the
# stress warning. This is a limitation of their 3D representation, not an
# import failure. Inspect the fit before interpreting apparent separation.
# Most year comparisons also show unequal dispersion in this dataset.
# Discuss this alongside PERMANOVA, rather than interpreting its p-value alone.
# Coordinates can rotate/reflect, and fits can vary with software versions.


# 10. Student exercises ---------------------------------------------------
# A. Compare Soil and Litter. Describe overlap and separation by round,
#    report stress, and distinguish a visual pattern from a tested effect.
# B. Change selected_guild to "Saprotroph", "EcM", and then "Yeast".
#    Re-run from section 1. How do sample numbers and patterns change?
# C. Change nmds_dimensions from 3 to 2. Compare stress and the Shepard plots.
#    With k = 3, also inspect NMDS1 against NMDS3 in the saved coordinates.
# D. Replace Year_group with Sampling_round in the categorical PERMANOVA.
#    Explain how the scientific question changes when 2017a and 2017b separate.
# E. Set run_turnover_extension to TRUE. Why is an ordinary regression on
#    all sample pairs an unsuitable way to obtain independent-replicate p-values?
# F. For an extreme EcM point, inspect rowSums(results$Litter$presence_absence)
#    before deciding whether it is an error or a sample with few detected OTUs.


# WORKBOOK PREPARATION AND COLUMN DICTIONARY -------------------------------
# Source: the "OTU table" sheet of Krtiny_ITS2_2013_to_2023.xlsx.
# OTU IDs have their pipe characters removed. Counts are not renormalised.
# Sample metadata were decoded once from the original sample names.
# Year = numeric calendar year; Sampling_round keeps 2017a and 2017b separate.
# Site = repeated sampling location; Substrate = Soil or Litter.
#
# Preparation follows the source Rmd, in this order:
#   1. Remove samples with fewer than 3,000 reads in the source count table.
#   2. Within each substrate, retain sites with all six sampling rounds.
#   3. Retain OTUs with more than one read across ALL retained samples.
#      This means more than one read, not occurrence in more than one sample.
#   4. Copy the matching taxonomy and ecological annotations to Taxonomy.
# The exact preparation totals and excluded sites are listed below.
# Original data: 310 samples and 24,927 OTUs.
# Prepared data: 264 samples and 9,488 OTUs.
# Soil: 23 sites x 6 rounds = 138 samples.
# Litter: 21 sites x 6 rounds = 126 samples.
# Excluded Soil sites: 19, 43, 63. Excluded Litter sites: 19, 20, 43, 46, 49.
# Seven source samples were below 3,000 reads:
#   2016_L49 (2146), 2017_L19 (1513), 2017_L20 (442), 2017_L46 (2192),
#   2020_L43 (1959), 2020_S43 (1327), 2023_S63 (2677).
# Soil site 19 also lacks the 2017b sample. After depth filtering, all other
# excluded site series lack at least one required round. In the retained
# subset, 2,576 zero-total OTUs and 12,863 one-read OTUs were removed.
# Retained reads: 792,000 before OTU filtering; 779,137 afterwards.
# Retained guild assignments: 3,143 Saprotroph, 929 EcM, 571 Yeast OTUs.
#
# Reads_before_OTU_filter = source count total, used for the 3,000-read rule.
# Reads_in_OTU_table = total after step 3. It can be below 3,000. Do not apply
# the depth threshold to it again. No further rarefaction was performed.
# This reproduces a teaching subset, not independent proof that low-total
# samples originally failed rarefaction: upstream OTU filtering may reduce sums.
#
# Taxonomy preserves Kingdom, Phylum, Class, Order, Family, Genus, Species,
# and Primary_lifestyle. Source values of 0 in these text fields were
# replaced with blank cells, which mean missing annotation.
# Growth_form copies the source column labelled "19" (a column NAME,
# not the nineteenth physical column), interpreted as in the source Rmd.
# Guild flags reproduce the original text-matching rules:
#   Yeast      : Growth_form contains "yeast", excluding "non-yeast".
#   EcM        : Primary_lifestyle contains "ectomycorrhizal".
#   Saprotroph : Primary_lifestyle contains "saprotroph", excluding yeasts.
# Matching ignores letter case and surrounding spaces. A 0 means not assigned
# by these annotations, not proof of ecological absence. The three flags
# need not cover all OTUs, and EcM is not forced to be mutually exclusive.
# The workbook is a fixed prepared dataset. If raw samples change, repeat
# the preparation instead of treating these preprocessing decisions as universal.
#
# This example covers ordination, PERMANOVA and dispersion diagnostics.
# Soil chemistry, vegetation predictors and mixed-effects decay models from
# the later parts of the Rmd are not required and are not included here.
#
# Further reading (official vegan documentation):
# https://vegandevs.github.io/vegan/reference/metaMDS.html
# https://vegandevs.github.io/vegan/reference/adonis.html
# https://vegandevs.github.io/vegan/reference/betadisper.html
# https://vegandevs.github.io/vegan/reference/permutations.html
