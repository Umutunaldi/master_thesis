library(ggplot2)

# Reads the complete comparison results as a tab-separated table [1]
summary_df <- read.delim(
  "summary_results/all_results.tsv",
  check.names = FALSE
)

# Reads every mapped interface contact [1]
contact_df <- read.delim(
  "interface_analysis/interface_contact_changes.tsv",
  check.names = FALSE
)

# Applies the final structural comparison limit
# Rows above 2 A are not used in any final network
selected_df <- summary_df[
  summary_df$Inter_CA_RMSD <= 2,
]

# Uses the comparison IDs to keep their matching contact rows
contact_df <- contact_df[
  contact_df$Comparison_ID %in% selected_df$Comparison_ID,
]

# Adds the protein side before the residue ID
# This prevents equal residue labels on the two chains from becoming one node
contact_df$Iso_Node_ID <- paste0(
  "iso__",
  contact_df$Iso_Residue
)
contact_df$Inter_Node_ID <- paste0(
  "inter__",
  contact_df$Inter_Residue
)

# Keeps one copy of every isoform residue node
iso_node_df <- unique(contact_df[c(
  "Comparison_ID",
  "Iso_Node_ID",
  "Iso_Residue"
)])

# Residue IDs have the form chain:position:amino_acid
# strsplit separates these three parts [2]
# rbind puts the separated parts into one matrix
iso_parts <- do.call(
  rbind,
  strsplit(
    iso_node_df$Iso_Residue,
    ":",
    fixed = TRUE
  )
)

iso_nodes <- data.frame(
  Comparison_ID = iso_node_df$Comparison_ID,
  Node_ID = iso_node_df$Iso_Node_ID,
  Protein_Side = "Iso",
  Residue_Position = as.numeric(iso_parts[, 2]),
  Residue_Label = paste0(iso_parts[, 3], iso_parts[, 2])
)

# Keeps one copy of every interactor residue node
inter_node_df <- unique(contact_df[c(
  "Comparison_ID",
  "Inter_Node_ID",
  "Inter_Residue"
)])

# Applies the same residue ID separation to the interactor side [2]
inter_parts <- do.call(
  rbind,
  strsplit(
    inter_node_df$Inter_Residue,
    ":",
    fixed = TRUE
  )
)

inter_nodes <- data.frame(
  Comparison_ID = inter_node_df$Comparison_ID,
  Node_ID = inter_node_df$Inter_Node_ID,
  Protein_Side = "Inter",
  Residue_Position = as.numeric(inter_parts[, 2]),
  Residue_Label = paste0(inter_parts[, 3], inter_parts[, 2])
)

# Combines both protein sides into one node table
node_df <- rbind(
  iso_nodes,
  inter_nodes
)

# Output folders for the complete and changed-only networks
all_network_dir <- paste0(
  "figures/final_2A/interaction_networks/",
  "all_contacts"
)
changed_network_dir <- paste0(
  "figures/final_2A/interaction_networks/",
  "changed_only"
)

dir.create(
  all_network_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  changed_network_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Draws two networks for every final comparison
# One contains all contacts and the other contains only changes
for (comparison_id in selected_df$Comparison_ID) {
  # Selects one comparison and its contacts
  result_row <- selected_df[
    selected_df$Comparison_ID == comparison_id,
  ]
  comparison_contact_df <- contact_df[
    contact_df$Comparison_ID == comparison_id,
  ]

  # Reads the names used in the title and legend
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
  shared_count <- sum(
    comparison_contact_df$Change_Type == "shared"
  )
  lost_count <- sum(
    comparison_contact_df$Change_Type == "lost_in_positive"
  )
  gained_count <- sum(
    comparison_contact_df$Change_Type == "gained_in_positive"
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

  # Makes the title used by both network versions
  figure_title <- paste0(
    iso_1, " (", y2h_1, ") vs ",
    iso_2, " (", y2h_2, ")",
    " with ", inter_name,
    " (", tpl_id, ":",
    tpl_iso_chain, ":",
    tpl_inter_chain, ")"
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

  # Makes the complete and changed-only network versions
  for (network_type in c("all_contacts", "changed_only")) {
    plot_contact_df <- comparison_contact_df

    # Removes shared contacts from the changed-only version
    if (network_type == "changed_only") {
      plot_contact_df <- plot_contact_df[
        plot_contact_df$Change_Type != "shared",
      ]

      contact_levels <- c(
        "lost_in_positive",
        "gained_in_positive"
      )
      contact_labels <- c(
        lost_label,
        gained_label
      )
      contact_colors <- c(
        "#D62728",
        "#009E73"
      )
      network_subtitle <- "Changed mapped contacts (shared contacts omitted)"
      output_dir <- changed_network_dir
      output_suffix <- "_changed_only_network.png"
    }

    # Keeps all contact types in the complete version
    if (network_type == "all_contacts") {
      contact_levels <- c(
        "shared",
        "lost_in_positive",
        "gained_in_positive"
      )
      contact_labels <- c(
        shared_label,
        lost_label,
        gained_label
      )
      contact_colors <- c(
        "#B3B3B3",
        "#D62728",
        "#009E73"
      )
      network_subtitle <- "All mapped contacts"
      output_dir <- all_network_dir
      output_suffix <- "_all_contacts_network.png"
    }

    # factor fixes the order used for drawing and for the legend [3]
    plot_contact_df$Contact_Type <- factor(
      plot_contact_df$Change_Type,
      levels = contact_levels,
      labels = contact_labels
    )
    names(contact_colors) <- contact_labels

    # Makes an invisible row for every legend category
    # This keeps a zero-count category in the legend
    legend_df <- data.frame(
      X = 0,
      Y = -1,
      X_End = 0,
      Y_End = -1,
      Contact_Type = factor(
        contact_labels,
        levels = contact_labels
      )
    )

    # Orders shared contacts first so changed contacts are drawn on top [4]
    plot_contact_df <- plot_contact_df[
      order(plot_contact_df$Contact_Type),
    ]

    # Keeps only residue nodes connected by the displayed contacts
    selected_nodes <- unique(c(
      plot_contact_df$Iso_Node_ID,
      plot_contact_df$Inter_Node_ID
    ))
    plot_node_df <- node_df[
      node_df$Comparison_ID == comparison_id &
        node_df$Node_ID %in% selected_nodes,
    ]

    # Orders isoform residues by sequence position on the left side [4]
    iso_plot_nodes <- plot_node_df[
      plot_node_df$Protein_Side == "Iso",
    ]
    iso_plot_nodes <- iso_plot_nodes[
      order(iso_plot_nodes$Residue_Position),
    ]
    iso_plot_nodes$X <- 0
    # Evenly spreads the residues from the top to the bottom [5]
    iso_plot_nodes$Y <- seq(
      1,
      0,
      length.out = nrow(iso_plot_nodes)
    )

    # Orders interactor residues by sequence position on the right side [4]
    inter_plot_nodes <- plot_node_df[
      plot_node_df$Protein_Side == "Inter",
    ]
    inter_plot_nodes <- inter_plot_nodes[
      order(inter_plot_nodes$Residue_Position),
    ]
    inter_plot_nodes$X <- 1
    # Uses the same top-to-bottom scale for the interactor side [5]
    inter_plot_nodes$Y <- seq(
      1,
      0,
      length.out = nrow(inter_plot_nodes)
    )

    # Combines both positioned node tables
    plot_node_df <- rbind(
      iso_plot_nodes,
      inter_plot_nodes
    )

    # match connects each residue ID with its node coordinates [6]
    # Isoform coordinates become the start and interactor coordinates become the end
    plot_contact_df$X <- plot_node_df$X[
      match(
        plot_contact_df$Iso_Node_ID,
        plot_node_df$Node_ID
      )
    ]
    plot_contact_df$Y <- plot_node_df$Y[
      match(
        plot_contact_df$Iso_Node_ID,
        plot_node_df$Node_ID
      )
    ]
    plot_contact_df$X_End <- plot_node_df$X[
      match(
        plot_contact_df$Inter_Node_ID,
        plot_node_df$Node_ID
      )
    ]
    plot_contact_df$Y_End <- plot_node_df$Y[
      match(
        plot_contact_df$Inter_Node_ID,
        plot_node_df$Node_ID
      )
    ]

    # Makes shared contacts lighter so lost and gained contacts remain visible
    plot_contact_df$Edge_Alpha <- 0.70
    plot_contact_df$Edge_Width <- 0.55
    plot_contact_df$Edge_Alpha[
      plot_contact_df$Change_Type == "shared"
    ] <- 0.18
    plot_contact_df$Edge_Width[
      plot_contact_df$Change_Type == "shared"
    ] <- 0.25

    # Gives crowded networks more vertical space
    # The larger protein side decides the minimum required height
    network_height <- max(
      7,
      nrow(iso_plot_nodes) * 0.13,
      nrow(inter_plot_nodes) * 0.13
    )

    # Draws every residue pair as a line between two nodes [7]
    network_plot <- ggplot() +
      geom_segment(
        data = legend_df,
        aes(
          x = X,
          y = Y,
          xend = X_End,
          yend = Y_End,
          color = Contact_Type
        ),
        linewidth = 0,
        alpha = 0,
        show.legend = TRUE
      ) +
      geom_segment(
        data = plot_contact_df,
        aes(
          x = X,
          y = Y,
          xend = X_End,
          yend = Y_End,
          color = Contact_Type,
          alpha = Edge_Alpha,
          linewidth = Edge_Width
        )
      ) +
      # Draws the residue nodes on both sides [8]
      geom_point(
        data = plot_node_df,
        aes(x = X, y = Y),
        shape = 21,
        size = 2.1,
        fill = "white",
        color = "black",
        stroke = 0.35
      ) +
      # Writes the isoform residue labels [9]
      geom_text(
        data = iso_plot_nodes,
        aes(x = X, y = Y, label = Residue_Label),
        hjust = 1.2,
        size = 2.2
      ) +
      # Writes the interactor residue labels [9]
      geom_text(
        data = inter_plot_nodes,
        aes(x = X, y = Y, label = Residue_Label),
        hjust = -0.2,
        size = 2.2
      ) +
      annotate(
        "text",
        x = 0,
        y = 1.05,
        label = "Mapped isoform residues",
        fontface = "bold",
        size = 3.1
      ) +
      annotate(
        "text",
        x = 1,
        y = 1.05,
        label = "Mapped interactor residues",
        fontface = "bold",
        size = 3.1
      ) +
      # Applies the fixed grey, red and green contact colors [10]
      scale_color_manual(
        values = contact_colors,
        drop = FALSE
      ) +
      scale_alpha_identity() +
      scale_linewidth_identity() +
      # Keeps the three contact types in one legend row [11]
      guides(
        color = guide_legend(
          nrow = 1,
          byrow = TRUE,
          override.aes = list(alpha = 1, linewidth = 1)
        )
      ) +
      # Adds horizontal room for residue labels without removing them [12]
      coord_cartesian(
        xlim = c(-0.35, 1.35),
        ylim = c(-0.03, 1.08),
        clip = "off"
      ) +
      labs(
        title = figure_title,
        subtitle = paste(
          network_subtitle,
          "| Inter C\u03B1 RMSD:",
          round(inter_rmsd, 3),
          "\u00C5"
        ),
        color = NULL
      ) +
      # Removes axes because the node positions are only a layout [13]
      theme_void(base_size = 10) +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.position = "bottom",
        plot.margin = margin(12, 45, 12, 45)
      )

    # Saves a high-resolution PNG image [14]
    ggsave(
      file.path(
        output_dir,
        paste0(figure_name, output_suffix)
      ),
      network_plot,
      width = 9,
      height = network_height,
      dpi = 300,
      bg = "white"
    )
  }
}


######### REFERENCES #########
# 1 - https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html
# 2 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/strsplit.html
# 3 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/factor.html
# 4 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/order.html
# 5 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/seq.html
# 6 - https://stat.ethz.ch/R-manual/R-devel/library/base/html/match.html
# 7 - https://ggplot2.tidyverse.org/reference/geom_segment.html
# 8 - https://ggplot2.tidyverse.org/reference/geom_point.html
# 9 - https://ggplot2.tidyverse.org/reference/geom_text.html
# 10 - https://ggplot2.tidyverse.org/reference/scale_manual.html
# 11 - https://ggplot2.tidyverse.org/reference/guide_legend.html
# 12 - https://ggplot2.tidyverse.org/reference/coord_cartesian.html
# 13 - https://ggplot2.tidyverse.org/reference/ggtheme.html
# 14 - https://ggplot2.tidyverse.org/reference/ggsave.html
