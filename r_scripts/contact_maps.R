library(ggplot2)

# Reads the complete comparison results [1]
summary_df <- read.delim(
  "summary_results/all_results.tsv",
  check.names = FALSE
)

# Reads every mapped residue contact [1]
contact_df <- read.delim(
  "interface_analysis/interface_contact_changes.tsv",
  check.names = FALSE
)

# Applies the final structural comparison limit
# Rows above 2 A are not used in any final contact map
selected_df <- summary_df[
  summary_df$Inter_CA_RMSD <= 2,
]

# Uses the comparison IDs to keep their matching contact rows
contact_df <- contact_df[
  contact_df$Comparison_ID %in% selected_df$Comparison_ID,
]

# Residue IDs have the form chain:position:amino_acid
# strsplit separates these three parts [2]
# rbind puts the separated parts into one matrix
iso_parts <- do.call(
  rbind,
  strsplit(
    contact_df$Iso_Residue,
    ":",
    fixed = TRUE
  )
)

inter_parts <- do.call(
  rbind,
  strsplit(
    contact_df$Inter_Residue,
    ":",
    fixed = TRUE
  )
)

# The second part of each residue ID is its sequence position
# These positions become the two contact-map axes
contact_df$Iso_Position <- as.numeric(iso_parts[, 2])
contact_df$Inter_Position <- as.numeric(inter_parts[, 2])

# Output folder for the final contact maps
contact_dir <- "figures/final_2A/contact_maps"
dir.create(
  contact_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Draws one contact map for every final comparison
for (comparison_id in selected_df$Comparison_ID) {
  # Selects one comparison and its contacts
  result_row <- selected_df[
    selected_df$Comparison_ID == comparison_id,
  ]
  plot_df <- contact_df[
    contact_df$Comparison_ID == comparison_id,
  ]

  # Reads the names used in the title and axes
  iso_1 <- result_row$Iso_1[1]
  iso_2 <- result_row$Iso_2[1]
  inter_name <- result_row$Inter_ID[1]
  tpl_id <- result_row$Base_Tpl[1]
  tpl_iso_chain <- result_row$Tpl_Iso_Chain[1]
  tpl_inter_chain <- result_row$Tpl_Inter_Chain[1]
  inter_rmsd <- result_row$Inter_CA_RMSD[1]

  # Converts the experimental result into a short title label
  y2h_1 <- ifelse(
    result_row$Interaction_Found_1[1] == "positive",
    "Y2H+",
    "Y2H-"
  )
  y2h_2 <- ifelse(
    result_row$Interaction_Found_2[1] == "positive",
    "Y2H+",
    "Y2H-"
  )

  # Counts the three contact types
  # Model 1 is Y2H negative and model 2 is Y2H positive
  shared_count <- sum(plot_df$Change_Type == "shared")
  lost_count <- sum(
    plot_df$Change_Type == "lost_in_positive"
  )
  gained_count <- sum(
    plot_df$Change_Type == "gained_in_positive"
  )

  # Adds the counts to the legend text
  shared_label <- paste0(
    "Shared (n = ", shared_count, ")"
  )
  lost_label <- paste0(
    "Lost in ", iso_2, " (n = ", lost_count, ")"
  )
  gained_label <- paste0(
    "Gained in ", iso_2, " (n = ", gained_count, ")"
  )

  # factor fixes the order used for drawing and for the legend [3]
  plot_df$Contact_Type <- factor(
    plot_df$Change_Type,
    levels = c(
      "shared",
      "lost_in_positive",
      "gained_in_positive"
    ),
    labels = c(
      shared_label,
      lost_label,
      gained_label
    )
  )

  # Draws shared contacts first so changed contacts appear on top [4]
  plot_df <- plot_df[
    order(plot_df$Contact_Type),
  ]

  # Uses grey for shared, red for lost and green for gained
  contact_colors <- c(
    "#B3B3B3",
    "#D62728",
    "#009E73"
  )
  names(contact_colors) <- c(
    shared_label,
    lost_label,
    gained_label
  )

  # Makes shared points lighter so changed contacts remain visible
  contact_alpha <- c(0.30, 0.85, 0.85)
  names(contact_alpha) <- c(
    shared_label,
    lost_label,
    gained_label
  )

  # Makes an invisible point for every legend category
  # This keeps a zero-count category in the legend
  legend_df <- data.frame(
    Iso_Position = plot_df$Iso_Position[1],
    Inter_Position = plot_df$Inter_Position[1],
    Contact_Type = factor(
      c(shared_label, lost_label, gained_label),
      levels = c(shared_label, lost_label, gained_label)
    )
  )

  # Makes a title with the Y2H results, interactor, template and chains
  figure_title <- paste0(
    iso_1, " (", y2h_1, ") vs ",
    iso_2, " (", y2h_2, ")",
    " with ", inter_name,
    " (", tpl_id, ":",
    tpl_iso_chain, ":",
    tpl_inter_chain, ")"
  )

  # Each point represents one mapped isoform-interactor residue pair [5]
  contact_plot <- ggplot(
    plot_df,
    aes(
      x = Iso_Position,
      y = Inter_Position,
      color = Contact_Type,
      alpha = Contact_Type
    )
  ) +
    geom_point(size = 1.6) +
    # Adds the invisible points needed for zero-count legend entries
    geom_point(
      data = legend_df,
      aes(
        x = Iso_Position,
        y = Inter_Position,
        color = Contact_Type
      ),
      inherit.aes = FALSE,
      size = 0,
      alpha = 0,
      show.legend = TRUE
    ) +
    # Applies the fixed contact colors and transparency [6]
    scale_color_manual(
      values = contact_colors,
      drop = FALSE
    ) +
    scale_alpha_manual(
      values = contact_alpha,
      drop = FALSE,
      guide = "none"
    ) +
    # Keeps the three contact types in one legend row [7]
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2.5, alpha = 1)
      )
    ) +
    # Adds the comparison title, RMSD and residue axes
    labs(
      title = figure_title,
      subtitle = paste(
        "All mapped contacts | Inter C\u03B1 RMSD:",
        round(inter_rmsd, 3),
        "\u00C5"
      ),
      x = "Mapped isoform residue position",
      y = "Mapped interactor residue position",
      color = NULL
    ) +
    # Uses a simple background suitable for a scientific figure [8]
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      legend.position = "bottom"
    )

  # Makes a filename without spaces or punctuation
  figure_name <- paste0(
    iso_1, "_vs_", iso_2,
    "_with_", inter_name,
    "_", tpl_id,
    "_", tpl_iso_chain,
    "_", tpl_inter_chain
  )
  figure_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    figure_name
  )

  # Saves a high-resolution PNG image [9]
  ggsave(
    file.path(
      contact_dir,
      paste0(figure_name, "_contact_map.png")
    ),
    contact_plot,
    width = 7.2,
    height = 5.8,
    dpi = 300,
    bg = "white"
  )
}


######### REFERENCES #########
# 1 - https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html
# 2 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/strsplit.html
# 3 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html
# 4 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/order.html
# 5 - https://ggplot2.tidyverse.org/reference/geom_point.html
# 6 - https://ggplot2.tidyverse.org/reference/scale_manual.html
# 7 - https://ggplot2.tidyverse.org/reference/guide_legend.html
# 8 - https://ggplot2.tidyverse.org/reference/ggtheme.html
# 9 - https://ggplot2.tidyverse.org/reference/ggsave.html
