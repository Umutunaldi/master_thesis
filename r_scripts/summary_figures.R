library(ggplot2)

# Run this script from the repository root directory

# Reads the complete comparison results as a tab-separated table [1]
summary_df <- read.delim(
  "summary_results/all_results.tsv",
  check.names = FALSE
)

# Converts the experimental results into short labels
# These labels are written beside their matching isoform names
summary_df$Y2H_1 <- ifelse(
  summary_df$Interaction_Found_1 == "positive",
  "Y2H+",
  "Y2H-"
)

summary_df$Y2H_2 <- ifelse(
  summary_df$Interaction_Found_2 == "positive",
  "Y2H+",
  "Y2H-"
)

# Makes one readable name for every structural comparison
# The template and directed chain pair separate repeated comparisons
summary_df$Comparison_Label <- paste0(
  summary_df$Iso_1,
  " (", summary_df$Y2H_1, ") vs ",
  summary_df$Iso_2,
  " (", summary_df$Y2H_2, ") | ",
  summary_df$Inter_ID,
  " | ",
  summary_df$Base_Tpl,
  " ",
  summary_df$Tpl_Iso_Chain,
  "-",
  summary_df$Tpl_Inter_Chain
)

# Calculates the two values that are only needed for plotting
summary_df$Mapped_Contact_Dissimilarity <- (
  1 - summary_df$Mapped_Contact_Jaccard
)

summary_df$Absolute_FoldX_Energy_Difference <- abs(
  summary_df$FoldX_Energy_Difference_Positive_Minus_Negative
)

# Labels comparisons using the final interactor Cα RMSD threshold
summary_df$RMSD_Group <- ifelse(
  summary_df$Inter_CA_RMSD <= 2,
  "At or below 2 Å",
  "Above 2 Å"
)

# Sets the order used in the RMSD legend [3]
summary_df$RMSD_Group <- factor(
  summary_df$RMSD_Group,
  levels = c("At or below 2 Å", "Above 2 Å")
)

# Keeps only comparisons that passed the final structural threshold
# These rows are used in every figure except the RMSD threshold figure
reliable_df <- summary_df[
  summary_df$Inter_CA_RMSD <= 2,
]

# Output folder for the summary figures
summary_dir <- "figures/summary"
dir.create(
  summary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Uses one simple paper style for every figure [2, 4]
# Axis titles and tick labels are bold to improve readability
# Legend titles such as FoldX direction and Contact type are also bold
paper_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.position = "bottom"
  )

# Uses the same result-direction colors in all FoldX figures [8]
direction_colors <- c(
  "yes" = "#009E73",
  "no" = "#D62728"
)


# Ranks comparisons from low to high interface similarity
ranked_df <- reliable_df[
  order(reliable_df$Inter_Interface_Jaccard),
]

# Keeps the ranked order on the vertical axis [3]
ranked_df$Comparison_Label <- factor(
  ranked_df$Comparison_Label,
  levels = ranked_df$Comparison_Label
)

# Draws one point for every reliable structural comparison [2, 5]
# labs gives the plot, axis and legend their readable names [9]
ranked_plot <- ggplot(
  ranked_df,
  aes(
    x = Inter_Interface_Jaccard,
    y = Comparison_Label,
    color = FoldX_Direction_Matches_Y2H
  )
) +
  geom_point(size = 2.3) +
  scale_color_manual(
    values = direction_colors,
    breaks = c("no", "yes"),
    labels = c(
      "no" = "Opposite to Y2H",
      "yes" = "Matches Y2H"
    )
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = "Interactor interface similarity",
    x = "Interactor interface Jaccard",
    y = NULL,
    color = "FoldX direction"
  ) +
  paper_theme +
  theme(
    axis.text.y = element_text(size = 6.5, face = "bold")
  )

# Saves the figure as a high-resolution PNG image [10]
ggsave(
  file.path(summary_dir, "ranked_interface_jaccard.png"),
  ranked_plot,
  width = 8.2,
  height = 9.5,
  dpi = 300,
  bg = "white"
)


# Puts the shared, lost and gained residue counts into one table
# rbind places the three contact types below each other
residue_df <- rbind(
  data.frame(
    Comparison_Label = reliable_df$Comparison_Label,
    Contact_Type = "Shared",
    Residue_Count = reliable_df$Shared_Inter_Interface_Residue_Count
  ),
  data.frame(
    Comparison_Label = reliable_df$Comparison_Label,
    Contact_Type = "Lost in positive",
    Residue_Count = reliable_df$Lost_Inter_Interface_Residue_Count
  ),
  data.frame(
    Comparison_Label = reliable_df$Comparison_Label,
    Contact_Type = "Gained in positive",
    Residue_Count = reliable_df$Gained_Inter_Interface_Residue_Count
  )
)

# Uses the same comparison order as the Jaccard figure [3]
residue_df$Comparison_Label <- factor(
  residue_df$Comparison_Label,
  levels = levels(ranked_df$Comparison_Label)
)

# Sets the contact order used in the legend and stacked bars [3]
residue_df$Contact_Type <- factor(
  residue_df$Contact_Type,
  levels = c(
    "Shared",
    "Lost in positive",
    "Gained in positive"
  )
)

# Uses gray for shared, red for lost and green for gained residues [8]
contact_colors <- c(
  "Shared" = "#B3B3B3",
  "Lost in positive" = "#D62728",
  "Gained in positive" = "#009E73"
)

# Draws the observed residue counts without another calculation [6]
residue_plot <- ggplot(
  residue_df,
  aes(
    x = Residue_Count,
    y = Comparison_Label,
    fill = Contact_Type
  )
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = contact_colors) +
  labs(
    title = "Interactor interface residue changes",
    x = "Interactor residue count",
    y = NULL,
    fill = "Contact type"
  ) +
  paper_theme +
  theme(
    axis.text.y = element_text(size = 6.5, face = "bold")
  )

# Saves the figure as a high-resolution PNG image [10]
ggsave(
  file.path(
    summary_dir,
    "interactor_interface_residue_changes.png"
  ),
  residue_plot,
  width = 8.2,
  height = 9.5,
  dpi = 300,
  bg = "white"
)


# Orders comparisons by the signed FoldX difference
foldx_df <- reliable_df[
  order(
    reliable_df$FoldX_Energy_Difference_Positive_Minus_Negative
  ),
]

# Keeps the FoldX ranking on the vertical axis [3]
foldx_df$Comparison_Label <- factor(
  foldx_df$Comparison_Label,
  levels = foldx_df$Comparison_Label
)

# Draws negative values left and positive values right of zero [5, 7]
foldx_plot <- ggplot(
  foldx_df,
  aes(
    x = FoldX_Energy_Difference_Positive_Minus_Negative,
    y = Comparison_Label,
    color = FoldX_Direction_Matches_Y2H
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#666666"
  ) +
  geom_point(size = 2.3) +
  scale_color_manual(
    values = direction_colors,
    breaks = c("no", "yes"),
    labels = c(
      "no" = "Opposite to Y2H",
      "yes" = "Matches Y2H"
    )
  ) +
  labs(
    title = "FoldX interaction-energy differences",
    x = "FoldX difference (Y2H+ minus Y2H-; kcal mol\u207B\u00B9)",
    y = NULL,
    color = "FoldX direction"
  ) +
  paper_theme +
  theme(
    axis.text.y = element_text(size = 6.5, face = "bold")
  )

# Saves the figure as a high-resolution PNG image [10]
ggsave(
  file.path(summary_dir, "foldx_energy_differences.png"),
  foldx_plot,
  width = 8.2,
  height = 9.5,
  dpi = 300,
  bg = "white"
)


# Compares mapped-contact dissimilarity with FoldX difference magnitude
# Only comparisons within the 2 Å threshold are used here
foldx_contact_plot <- ggplot(
  reliable_df,
  aes(
    x = Mapped_Contact_Dissimilarity,
    y = Absolute_FoldX_Energy_Difference,
    color = FoldX_Direction_Matches_Y2H
  )
) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_color_manual(
    values = direction_colors,
    breaks = c("no", "yes"),
    labels = c(
      "no" = "Opposite to Y2H",
      "yes" = "Matches Y2H"
    )
  ) +
  labs(
    title = "Contact dissimilarity and FoldX difference",
    x = "Mapped contact dissimilarity (1 - Jaccard)",
    y = "Absolute FoldX energy difference (kcal mol\u207B\u00B9)",
    color = "FoldX direction"
  ) +
  paper_theme

# Saves the figure as a high-resolution PNG image [10]
ggsave(
  file.path(summary_dir, "contact_dissimilarity_vs_foldx.png"),
  foldx_contact_plot,
  width = 6.2,
  height = 4.8,
  dpi = 300,
  bg = "white"
)


# Uses green for comparisons within the threshold and red above it [8]
rmsd_colors <- c(
  "At or below 2 Å" = "#009E73",
  "Above 2 Å" = "#D62728"
)

# Shows all comparisons so the 2 Å threshold remains visible [5, 7]
rmsd_plot <- ggplot(
  summary_df,
  aes(
    x = Inter_CA_RMSD,
    y = Inter_Interface_Jaccard,
    color = RMSD_Group
  )
) +
  geom_vline(
    xintercept = 2,
    linetype = "dashed",
    color = "#666666"
  ) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_color_manual(
    values = rmsd_colors,
    breaks = c("At or below 2 Å", "Above 2 Å")
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Interactor RMSD and interface similarity",
    x = "Interactor Cα RMSD (Å)",
    y = "Interactor interface Jaccard",
    color = "RMSD group"
  ) +
  paper_theme

# Saves the figure as a high-resolution PNG image [10]
ggsave(
  file.path(summary_dir, "rmsd_vs_interface_jaccard.png"),
  rmsd_plot,
  width = 6.2,
  height = 4.8,
  dpi = 300,
  bg = "white"
)


######### REFERENCES #########
# 1 - https://www.stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html
# 2 - https://ggplot2.tidyverse.org/reference/ggplot.html
# 3 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html
# 4 - https://ggplot2.tidyverse.org/reference/theme.html
# 5 - https://ggplot2.tidyverse.org/reference/geom_point.html
# 6 - https://ggplot2.tidyverse.org/reference/geom_bar.html
# 7 - https://ggplot2.tidyverse.org/reference/geom_abline.html
# 8 - https://ggplot2.tidyverse.org/reference/scale_manual.html
# 9 - https://ggplot2.tidyverse.org/reference/labs.html
# 10 - https://ggplot2.tidyverse.org/reference/ggsave.html
