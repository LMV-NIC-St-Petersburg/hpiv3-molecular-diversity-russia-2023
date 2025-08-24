library(tidyverse)
library(tidytree)
library(ggtree)
library(ggtreeExtra)
library(ggstar)

dataframe_to_treedata <- function(
  tree, metadata
) {

  if (missing(tree)) {
    stop("Parameter 'tree' is mandatory and must be specified.")
  }

  if (missing(metadata)) {
    stop("Parameter 'metadata' is mandatory and must be specified.")
  }

  if (!inherits(metadata, "data.frame") && !inherits(metadata, "tbl_df")) {
    stop("Parameter 'metadata' must be of data.frame or tibble format.")
  }

  if (!inherits(tree, "phylo") && !inherits(tree, "treedata")) {
    stop("Parameter 'tree' must be of phylo or treedata format.")
  }

  if (inherits(tree, "treedata")) {
    phylo_tree <- tree@phylo
  } else {
    phylo_tree <- tree
  }

  meta_length <- metadata$label |> length()
  tree_length <- phylo_tree$tip.label |> length()

  tree_original_treedata <- full_join(
    tree |>
      as_tibble(),
    metadata |>
      filter(
        label %in% phylo_tree$tip.label
      ),
    by = "label"
  ) |> as.treedata()

  tree_data_length <- tree_original_treedata@phylo$tip.label |> length()

  if (tree_length == tree_data_length) {
    cat(crayon::bgGreen("[SUCCESS]"), crayon::green("---- Tips has been checked successfully ---\n"))
    cat(crayon::blue(glue::glue("Tree tips ({tree_length}) == Tree_data tips ({tree_data_length})\n")))
  } else {
    cat(crayon::bgYellow("[WARNING]"), crayon::yellow("---- Tips are not equal ---\n"))
    cat(crayon::blue(glue::glue("Tree tips ({tree_length}) != Tree_data tips ({tree_data_length})\n")))
  }

  tree_original_treedata
}

add_multiple_strips_flexible <- function(plot_object, annotations_df, 
                                       default_offset = 0.1, 
                                       default_offset.text = 0.01,
                                       default_align = TRUE, 
                                       default_barsize = 0.25,
                                       default_fontsize = 3.88
                                      ) {
  
  strip_layers <- purrr::pmap(annotations_df, function(taxa1, taxa2, label, colour,
                                                      offset = default_offset,
                                                      offset.text = default_offset.text,
                                                      align = default_align, fontsize = default_fontsize,
                                                      barsize = default_barsize, ...) {
    geom_strip(
      taxa1 = taxa1,
      taxa2 = taxa2, 
      label = label,
      offset = offset,
      offset.text = offset.text,
      align = align,
      barsize = barsize,
      color = colour
    )
  })
  
  plot_object + strip_layers
}

########################### Fig 3 ################################

metadata <- readxl::read_excel(
  "Supplementary table.xlsx", skip = 2,
  .name_repair = janitor::make_clean_names, sheet = 3
) |>
  mutate(
    year = as.numeric(year),
    label = accession_number
    # sub_cluster = cluster,
    # cluster = str_extract(cluster, "[A-Z]")
  )

tree_original <- treeio::read.raxml(
  "Fig3_data_boost.nwk"
)

tree_original@data$bootstrap <- tree_original@data$bootstrap |>
  enframe() |>
  mutate(
    bootstrap = case_when(
      value >= 70 & value < 90 ~ glue::glue("70 ≤ BS < 90"),
      value >= 90 ~ "90 ≤ BS ≤ 100"
    ),
    bootstrap = factor(bootstrap, levels = c("70 ≤ BS < 90", "90 ≤ BS ≤ 100"))
  ) |> pull(bootstrap)


cluster_colors <- c("A" = "#d62728", "B" = "#ff7f0e", "C" = "#1f77b4", "D" = "#2ca02c")

subcluster_colors <- c("C1" = "#e377c2", "C2" = "#ffd700", "C3" = "#17becf", 
                       "C4" = "#8c564b", "C5" = "#7f7f7f")

p <- ggtree(
  dataframe_to_treedata(tree_original, metadata),
  layout = "fan", open.angle = 20, size = 0.25
) +
  # geom_tippoint(
  #   aes(colour = geo_loc_name)
  # ) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = label, fill = cluster),
    width = 0.005,
    color = NA,
    offset = 0.2,
    pwidth = 0
  ) +
  scale_fill_manual(
    values = cluster_colors
  ) +
  labs(
    fill = "Cluster"
  ) +
  ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = label, fill = sub_cluster),
    width = 0.005,
    color = NA,
    offset = 0.046,
    pwidth = 0
  ) +
  scale_fill_manual(
    values = subcluster_colors,
    na.translate = FALSE
  ) +
  geom_nodepoint(
    aes(
      colour = bootstrap
    ),
    size = 0.25
  ) +
  scale_colour_manual(
    "Node bootstrap\nsupport values",
    values = c("#D89D6A", "#6D454C"),
    na.translate = FALSE
  ) +
  ggnewscale::new_scale_colour() +
  geom_treescale(
    fontsize = 3,
    linesize = 1,
    offset = 10,
    x = 0.077,
    # y = 100
  ) +
  labs(
    fill = "Subclusters"
  ) +
  theme_tree()


ggsave(
  "Fig_3.png",
  p,
  dpi = 800,
  scale = 1
)


########################### Fig 4 ################################


metadata <- readxl::read_excel(
  "Supplementary table.xlsx", skip = 2,
  .name_repair = janitor::make_clean_names, sheet = 4
) |>
  mutate(
    year = as.numeric(year),
    label = accession_number,
    label_text = glue::glue("{label}/{year}/{geo_loc_name}"),
    sub_data = ifelse(
      !is.na(sub_cluster_label),
      glue::glue("{sub_cluster}{sub_cluster_label}"),
      NA_character_
    ),
    Russia_check = ifelse(
      geo_loc_name == "Russia",
      "Russia",
      "Other"
    )
  )

tree_original <- treeio::read.raxml(
  "fig4_boot.tree"
)

tree_original@phylo$tip.label <- tree_original@phylo$tip.label |>
  str_replace_all("'", "") |>
  str_split("_") |> lapply(first) |> unlist()

tree_original@data$bootstrap <- tree_original@data$bootstrap |>
  enframe() |>
  mutate(
    bootstrap = case_when(
      value >= 70 & value < 90 ~ glue::glue("70 ≤ BS < 90"),
      value >= 90 ~ "90 ≤ BS ≤ 100"
    ),
    bootstrap = factor(bootstrap, levels = c("70 ≤ BS < 90", "90 ≤ BS ≤ 100"))
  ) |> pull(bootstrap)



p <- ggtree(
  tr = dataframe_to_treedata(tree_original, metadata),
  size = 0.25
) +
  # geom_nodelab(
  #   aes(
  #     x = branch,
  #     label = bootstrap
  #   ),
  #   size = 2,
  #   color = "black",
  #   nudge_y = 1
  # ) +
  geom_nodepoint(
    aes(
      colour = bootstrap
    ),
    size = 1
  ) +
  # geom_label(
  #   aes(label = bootstrap),
  #   nudge_y = 1,
  #   size = 1.5
  # ) +
  # geom_label2(
  #   aes(label = bootstrap)
  # )
  # geom_point2(
  #   aes(colour = bootstrap),
  #   size = 1
  # ) +
  scale_colour_manual(
    "Node bootstrap\nsupport values",
    values = c("#D89D6A", "#6D454C"),
    na.translate = FALSE
  ) +
  ggnewscale::new_scale_colour() +
  # geom_tippoint(
  #   aes(colour = geo_loc_name),
  #   size = 1
  # ) +
  geom_tiplab(
    aes(
      label = label_text,
      colour = Russia_check
    ),
    size = 1,
    show.legend = FALSE,
    align = TRUE,
    linesize = 0.1
  ) +
  scale_colour_manual(
    values = c(
      "black", "purple"
    )
  ) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = label, fill = cluster),
    width = 0.002,
    color = NA,
    offset = 0.16,
    pwidth = 0
  ) +
  scale_fill_manual(
    values = cluster_colors
  ) +
  labs(
    fill = "Cluster"
  ) +
  geom_treescale(fontsize = 3, linesize = 1, y = -2, x = 0.01)


subcluster_colors <- c("C1" = "#e377c2", "C2" = "#ffd700", "C3" = "#17becf", 
                       "C4" = "#8c564b", "C5" = "#7f7f7f")


clade_annotations <- tibble::tribble(
  ~taxa1,      ~taxa2,      ~label,      ~colour,
  "EU814623",  "EU814623",  "C4",  "#8c564b",
  "OR189470",  "MW654351",  "C3a",  "#17becf",
  "KT796417",  "KY973585",  "C3d",  "#17becf",
  "PP858578",  "KJ672565",  "C3e",  "#17becf",
  "OM240920",  "GU732154",  "C3f",  "#17becf",
  "KF530253",  "PP858551",  "C3b",  "#17becf",
  "GU732160",  "GU732149",  "C3c",  "#17becf",
  "PP858548",  "OR189441",  "C5b",  "#7f7f7f",
  "PP858552",  "KF687320",  "C5a",  "#7f7f7f",
  "KY973580",  "AB189960",  "C2",  "#ffd700",
  "JF416795",  "EU814624",  "C1a",  "#e377c2",
  "AB736166",  "HM460886",  "C1b",  "#e377c2",
  "MW654321",  "KY684754",  "C1c",  "#e377c2",
)


ppp <- add_multiple_strips_flexible(
  p, clade_annotations,
  default_offset = 0.021,
  default_offset.text = -0.008,
  default_barsize = 1
)

pppp <- ppp +
  geom_strip(
    "PV498575",
    "OR189479",
    "D",
    offset = 0.008,
    offset.text = 0.001,
    barsize = 1
  )

ggsave(
  "Fig_4.png",
  pppp,
  dpi = 600, scale = 1.2
)


