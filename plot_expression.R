library(tidyverse)
library(geojsonio)

# Renaming function
f_rename <- function(.data) {
  .data %>%
    rename_with(~ str_to_lower(.) %>%
      str_replace_all("\\s", "_") %>%
      str_replace_all("-", "_"))
}

# List of stages to iterate over
stages <- c("gastrula", "neurula")

# Use map to iterate over stages
map(stages, function(stage_name) {
  # Read geojson file
  spdf <- geojson_read(paste0(stage_name, ".geojson"), what = "sp") %>%
    sf::st_as_sf()

  # Process data file
  df <- read_tsv("output.txt", na = "None") %>%
    f_rename() %>%
    filter(str_detect(territory, "cell pair") & str_detect(stage, stage_name)) %>%
    mutate(territory = str_replace(territory, " .*", "")) %>%
    group_by(territory)
  # Join data and plot
  df_to_plot <- left_join(spdf, df %>% distinct(stage, gene) %>% count(), by = c("name" = "territory"))

  plot <- ggplot(data = df_to_plot) +
    geom_sf(aes(fill = n)) +
    geom_sf_text(aes(label = name)) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  ggsave(paste0("choropleth_", stage_name, ".png"))

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


  ggsave(paste0("matrix_", stage_name, ".png"))
})
