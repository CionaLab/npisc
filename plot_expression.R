library(tidyverse)
library(geojsonio)
library(infotheo)
library(ggdendro)

# Renaming function
f_rename <- function(.data) {
  .data %>%
    rename_with(~ str_to_lower(.) %>%
      str_replace_all("\\s", "_") %>%
      str_replace_all("-", "_"))
}

# List of stages to iterate over
stages <- c(
  "mid gastrula",
  "late gastrula",
  "early neurula",
  "mid neurula",
  "late neurula"
)

# Use map to iterate over stages
map(stages, function(stage_name) {
  file_stage <- str_replace_all(stage_name, " ", "_")
  # Read geojson file
  spdf <- geojson_read(paste0(file_stage, ".geojson"), what = "sp") %>%
    sf::st_as_sf()

  # Process data file
  df <- read_tsv("pass_02_full.tsv", na = "None") %>%
    f_rename() %>%
    filter(
      str_detect(territory, "[Aa]\\d+\\.\\d+ cell pair$") &
        str_detect(stage, stage_name)
    ) %>%
    mutate(territory = str_replace(territory, " .*", "")) %>%
    group_by(territory)
  # Join data and plot
  df_to_plot <- left_join(
    spdf,
    df %>% distinct(stage, gene) %>% count(),
    by = c("name" = "territory")
  )

  plot <- ggplot(data = df_to_plot) +
    geom_sf(aes(fill = n)) +
    geom_sf_text(aes(label = name), color = "white") +
    scale_fill_gradient(low = "blue", high = "red", limits = c(1, 20)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  file_stage <- paste0(file_stage, "_orig")
  ggsave(paste0("choropleth_", file_stage, ".png"))

  plot <- ggplot(df, aes(x = gene, y = territory)) +
    geom_raster() +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme(
      # Rotate the x-axis lables so they are legible
      axis.text.x = element_text(angle = 270, hjust = 0),
      # Force the plot into a square aspect ratio
      aspect.ratio = 1,
      # Hide the legend (optional)
      legend.position = "none"
    )

  ggsave(paste0("matrix_", file_stage, ".png"))

  m1 <- df %>%
    group_by(gene, territory) %>%
    count() %>%
    mutate_at(vars(n), ~ ifelse(. > 0, 1, 0)) %>%
    pivot_wider(names_from = gene, values_from = n, values_fill = 0) %>%
    column_to_rownames(var = "territory") %>%
    proxy::dist(method = "Jaccard")

  den1 <- hclust(m1) %>% as.dendrogram()
  plot <- ggdendrogram(data = den1, rotate = TRUE)
  ggsave(paste0("dendrogram_", file_stage, ".png"))

  df1 <- as.matrix(m1) %>%
    as.data.frame() %>%
    rownames_to_column("blastomere_1") %>%
    pivot_longer(
      -blastomere_1,
      names_to = "blastomere_2",
      values_to = "jaccard_distance"
    )

  plot <- ggplot(df1, aes(x = blastomere_1, y = blastomere_2)) +
    geom_tile(aes(fill = jaccard_distance)) +
    scale_fill_gradient2() +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE, limits = rev) +
    theme(
      # Rotate the x-axis lables so they are legible
      axis.text.x = element_text(angle = 270, hjust = 0),
      # Force the plot into a square aspect ratio
      aspect.ratio = 1,
    )

  ggsave(paste0("ji_", file_stage, ".png"))
})

# Use map to iterate over stages
map(stages, function(stage_name) {
  file_stage <- str_replace_all(stage_name, " ", "_")
  # Read geojson file
  spdf <- geojson_read(paste0(file_stage, ".geojson"), what = "sp") %>%
    sf::st_as_sf()

  # Process data file
  df <- read_tsv("pass_02.tsv", na = "None") %>%
    f_rename() %>%
    filter(
      str_detect(territory_eq, "[Aa]\\d+\\.\\d+$") &
        str_detect(stage, stage_name)
    ) %>%
    group_by(territory_eq)
  # Join data and plot
  df_to_plot <- left_join(
    spdf,
    df %>% distinct(stage, gene) %>% count(),
    by = c("name" = "territory_eq")
  )

  plot <- ggplot(data = df_to_plot) +
    geom_sf(aes(fill = n)) +
    geom_sf_text(aes(label = name)) +
    geom_sf_text(aes(label = name), color = "white") +
    scale_fill_gradient(low = "blue", high = "red", limits = c(1, 20)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  ggsave(paste0("choropleth_", file_stage, ".png"))

  plot <- ggplot(df, aes(x = gene, y = territory_eq)) +
    geom_raster() +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme(
      # Rotate the x-axis lables so they are legible
      axis.text.x = element_text(angle = 270, hjust = 0),
      # Force the plot into a square aspect ratio
      aspect.ratio = 1,
      # Hide the legend (optional)
      legend.position = "none"
    )

  ggsave(paste0("matrix_", file_stage, ".png"))

  m1 <- df %>%
    group_by(gene, territory_eq) %>%
    count() %>%
    mutate_at(vars(n), ~ ifelse(. > 0, 1, 0)) %>%
    pivot_wider(names_from = gene, values_from = n, values_fill = 0) %>%
    column_to_rownames(var = "territory_eq") %>%
    proxy::dist(method = "Jaccard")

  den1 <- hclust(m1) %>% as.dendrogram()
  plot <- ggdendrogram(data = den1, rotate = TRUE)
  ggsave(paste0("dendrogram_", file_stage, ".png"))

  df1 <- as.matrix(m1) %>%
    as.data.frame() %>%
    rownames_to_column("blastomere_1") %>%
    pivot_longer(
      -blastomere_1,
      names_to = "blastomere_2",
      values_to = "jaccard_distance"
    )

  plot <- ggplot(df1, aes(x = blastomere_1, y = blastomere_2)) +
    geom_tile(aes(fill = jaccard_distance)) +
    scale_fill_gradient2() +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE, limits = rev) +
    theme(
      # Rotate the x-axis lables so they are legible
      axis.text.x = element_text(angle = 270, hjust = 0),
      # Force the plot into a square aspect ratio
      aspect.ratio = 1,
    )

  ggsave(paste0("ji_", file_stage, ".png"))
})
