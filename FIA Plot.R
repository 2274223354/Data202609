##==============================================================================
## FIA Plot ---------
##==============================================================================
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggplot2)
library(ggsignif)

#### Figure 1a #############################################################


library(data.table)
library(dplyr)
library(sf)
library(ggplot2)
library(ggspatial)
library(grid)
library(rnaturalearth)



data <- fread(
  "/data/Data.csv"
)


data_plot <- data %>% dplyr::select(-V1) %>%
  filter(
    !is.na(LON),
    !is.na(LAT),
    !is.na(NA_L2NAME)
  )




plot_sf <- st_as_sf(
  data_plot,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)




eco <- st_read(
  "/data/us_eco_l3_state_boundaries.shp",
  quiet = TRUE
)




eco_l2 <- eco %>%
  group_by(NA_L2NAME) %>%
  dplyr::summarise()


eco_l2 <- st_simplify(
  eco_l2,
  dTolerance = 5000
)


eco_l2 <- st_transform(
  eco_l2,
  4326
)




sample_eco <- data_plot %>%
  distinct(NA_L2NAME) %>%
  pull(NA_L2NAME)

sample_eco <- sort(sample_eco)

print(sample_eco)
length(sample_eco)



eco_l2_plot <- eco_l2 %>%
  mutate(
    sampled_ecoregion = if_else(
      NA_L2NAME %in% sample_eco,
      NA_L2NAME,
      NA_character_
    )
  )





eco_names_shp <- unique(eco_l2$NA_L2NAME)

unmatched_eco <- setdiff(
  sample_eco,
  eco_names_shp
)



eco_colors <- c(
  
  "ATLANTIC HIGHLANDS" =
    "#b3dcff",
  
  "CENTRAL USA PLAINS" =
    "#fdd666",
  
  "COLD DESERTS" =
    "#F6E2AE",
  
  "EVERGLADES" =
    "#ebffac",
  
  "MARINE WEST COAST FOREST" =
    "#C7E4A2",
  
  "MEDITERRANEAN CALIFORNIA" =
    "#D2D99B",
  
  "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS" =
    "#CF6765",
  
  "MIXED WOOD PLAINS" =
    "#FFF0A6",
  
  "MIXED WOOD SHIELD" =
    "#AFC3E5",
  
  "OZARK/OUACHITA-APPALACHIAN FORESTS" =
    "#ffc2e5",
  
  "SOUTH CENTRAL SEMI-ARID PRAIRIES" =
    "#C7C7C7",
  
  "SOUTHEASTERN USA PLAINS" =
    "#c1f1fc",
  
  "TEMPERATE PRAIRIES" =
    "#D9D9D9",
  
  "TEXAS-LOUISIANA COASTAL PLAIN" =
    "#ebffac",
  
  "UPPER GILA MOUNTAINS" =
    "#7F95DB",
  
  "WEST-CENTRAL SEMI-ARID PRAIRIES" =
    "#D5CDA3",
  
  "WESTERN CORDILLERA" =
    "#76daff",
  
  "WESTERN SIERRA MADRE PIEDMONT" =
    "#C9DFA0"
)



eco_labels <- c(
  
  "ATLANTIC HIGHLANDS" =
    "Atlantic Highlands",
  
  "CENTRAL USA PLAINS" =
    "Central USA Plains",
  
  "COLD DESERTS" =
    "Cold Deserts",
  
  "EVERGLADES" =
    "Everglades",
  
  "MARINE WEST COAST FOREST" =
    "Marine West Coast Forest",
  
  "MEDITERRANEAN CALIFORNIA" =
    "Mediterranean California",
  
  "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS" =
    "Mississippi Alluvial and Southeast USA Coastal Plains",
  
  "MIXED WOOD PLAINS" =
    "Mixed Wood Plains",
  
  "MIXED WOOD SHIELD" =
    "Mixed Wood Shield",
  
  "OZARK/OUACHITA-APPALACHIAN FORESTS" =
    "Ozark/Ouachita-Appalachian Forests",
  
  "SOUTH CENTRAL SEMI-ARID PRAIRIES" =
    "South Central Semi-Arid Prairies",
  
  "SOUTHEASTERN USA PLAINS" =
    "Southeastern USA Plains",
  
  "TEMPERATE PRAIRIES" =
    "Temperate Prairies",
  
  "TEXAS-LOUISIANA COASTAL PLAIN" =
    "Texas-Louisiana Coastal Plain",
  
  "UPPER GILA MOUNTAINS" =
    "Upper Gila Mountains",
  
  "WEST-CENTRAL SEMI-ARID PRAIRIES" =
    "West-Central Semi-Arid Prairies",
  
  "WESTERN CORDILLERA" =
    "Western Cordillera",
  
  "WESTERN SIERRA MADRE PIEDMONT" =
    "Western Sierra Madre Piedmont"
)



### 10. 如果存在颜色表之外的生态区，自动补色


eco_present <- unique(
  eco_l2_plot$sampled_ecoregion[
    !is.na(eco_l2_plot$sampled_ecoregion)
  ]
)

missing_colors <- setdiff(
  eco_present,
  names(eco_colors)
)

if(length(missing_colors) > 0){
  
  extra_cols <- c(
    "#B7D5C4",
    "#E3B7A0",
    "#B8A7CE",
    "#A9C2B8",
    "#E8D59D",
    "#C6B0A0"
  )
  
  extra_cols <- rep(
    extra_cols,
    length.out = length(missing_colors)
  )
  
  names(extra_cols) <- missing_colors
  
  eco_colors <- c(
    eco_colors,
    extra_cols
  )
}



missing_labels <- setdiff(
  eco_present,
  names(eco_labels)
)

if(length(missing_labels) > 0){
  
  auto_labels <- tools::toTitleCase(
    tolower(missing_labels)
  )
  
 
  auto_labels <- gsub(
    "\\bUsa\\b",
    "USA",
    auto_labels
  )
  
  names(auto_labels) <- missing_labels
  
  eco_labels <- c(
    eco_labels,
    auto_labels
  )
}




eco_colors_use <- eco_colors[
  names(eco_colors) %in% eco_present
]

eco_labels_use <- eco_labels[
  names(eco_colors_use)
]



usa <- ne_states(
  country = "United States of America",
  returnclass = "sf"
)

usa <- st_transform(
  usa,
  4326
)


p <- ggplot() +
  
 
  geom_sf(
    data = usa,
    fill = "grey97",
    color = NA
  ) +
  
 
  geom_sf(
    data = eco_l2,
    fill = "grey93",
    color = "grey72",
    linewidth = 0.20
  ) +
  
  
  geom_sf(
    data = eco_l2_plot %>%
      filter(!is.na(sampled_ecoregion)),
    aes(fill = sampled_ecoregion),
    color = "grey55",
    linewidth = 0.22,
    alpha = 0.90
  ) +
  
  
  geom_sf(
    data = plot_sf,
    aes(color = "Forest inventory plots"),
    size = 0.1,
    alpha = 0.5
  ) +
  
 
  scale_fill_manual(
    name = "Ecoregion",
    values = eco_colors_use,
    labels = eco_labels_use,
    drop = TRUE,
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0,
      ncol = 1,
      keyheight = unit(0.34, "cm"),
      keywidth  = unit(0.50, "cm")
    )
  ) +
  
  
  scale_color_manual(
    name = "",
    values = c("Forest inventory plots" = "black"),
    guide = guide_legend(
      override.aes = list(
        size = 0.5,
        alpha = 1
      )
    )
  ) +
  
  
  coord_sf(
    xlim = c(-125, -65),
    ylim = c(24, 50),
    expand = FALSE
  ) +
  
 
  annotation_scale(
    location = "bl",
    width_hint = 0.18,
    text_cex = 0.75,
    text_family = "serif",
    line_width = 0.45
  ) +
  
 
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_minimal,
    height = unit(0.65, "cm"),
    width = unit(0.65, "cm")
  ) +
  
  theme_classic() +
  
  theme(
    
    axis.title = element_blank(),
    
    axis.text = element_text(
      size = 10,
      family = "serif",
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.30
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 10,
      family = "serif",
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8,
      family = "serif"
    ),
    
    legend.key = element_rect(
      color = NA
    ),
    
    legend.spacing.y = unit(0.03, "cm"),
    
    plot.margin = ggplot2::margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )


p



ggsave(
  filename = "FIA_sampling.png",
  path = "/Figure/",
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  plot = p
)

##==============================================================================
## Figure 1 b-d ---------
##==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(hexbin)
library(grid)



climate_plot_data <- fread(
  "/data/Data.csv"
) %>%
  as.data.frame()

climate_plot_data <- climate_plot_data %>%
  dplyr::select(-V1) %>%
  filter(Carbon_Mg_ha > 0)

climate_plot_data <- climate_plot_data %>%
  filter(!is.na(CWD_stand_age)) %>%
  filter(!is.na(MAT_stand_age)) %>%
  filter(!is.na(MAP_stand_age)) 

climate_plot_data <- climate_plot_data %>%
  transmute(
    MAT = MAT_stand_age,       # °C
    MAP = MAP_stand_age,     # mm yr-1
    CWD = CWD_stand_age,               # mm
    AGC = Carbon_Mg_ha            # Mg C ha-1
  ) %>%
  filter(
    complete.cases(MAT, MAP, CWD, AGC),
    AGC >= 0
  )


summary(climate_plot_data[, c("MAT", "MAP", "CWD", "AGC")])
range(climate_plot_data$MAT, na.rm = TRUE)
range(climate_plot_data$MAP, na.rm = TRUE)
range(climate_plot_data$CWD, na.rm = TRUE)
range(climate_plot_data$AGC, na.rm = TRUE)


AGC_low <- 0
AGC_high <-max(climate_plot_data$AGC)

AGC_high <- quantile(
  climate_plot_data$AGC,
  probs = 0.98,
  na.rm = TRUE
)
AGC_high



AGC_scale <- scale_fill_gradientn(
  colours = c(
    "#5B6CCF",   # low AGC
    "#71C9E8",
    "#BCE8D0",
    "#F6E7A1",
    "#F3A64A",
    "#D94A35"    # high AGC
  ),
  limits = c(AGC_low, AGC_high),
  oob = scales::squish,
  name = "Aboveground carbon storage (Mg C ha\u207b\u00b9)",
 
  breaks = pretty(c(AGC_low, AGC_high), n = 8)  
)



theme_climate <- theme_classic() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = grid::unit(0.15, "cm"),
    axis.text = element_text(size = 13, family = "serif", colour = "black"),
    axis.title = element_text(size = 14, family = "serif", colour = "black"),
    plot.title = element_text(size = 15, family = "serif", colour = "black", hjust = 0),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13, family = "serif", colour = "black"),
    legend.text = element_text(size = 12, family = "serif", colour = "black"),
    plot.margin = ggplot2::margin(t = 5, r = 8, b = 5, l = 8)
  )


P_MAT_MAP <- ggplot(
  climate_plot_data,
  aes(x = MAP, y = MAT, z = AGC)
) +
  stat_summary_hex(
    aes(fill = after_stat(value)),
    bins = 35,
    fun = mean,
    na.rm = TRUE
  ) +
  AGC_scale +                               
  labs(
    x = "Mean annual precipitation (mm)",
    y = "Mean annual temperature (\u00b0C)",
    title = "(b)"
  ) +
  theme_climate +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = grid::unit(8, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  )


  
P_MAT_CWD <- ggplot(
    climate_plot_data,
    aes(
      x = CWD,
      y = MAT,
      z = AGC
    )
  ) +
  
  stat_summary_hex(
    aes(
      fill = after_stat(value)
    ),
    bins = 35,
    fun = mean,
    na.rm = TRUE
  ) +
  
  AGC_scale +
  
  labs(
    x = "Climatic water deficit (mm)",
    y = "Mean annual temperature (\u00b0C)",
    title = "(c)"
  ) +
  
  theme_climate +
  
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = grid::unit(8, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  )
P_MAT_CWD


  
P_MAP_CWD <- ggplot(
    climate_plot_data,
    aes(
      x = MAP,
      y = CWD,
      z = AGC
    )
  ) +
  
  stat_summary_hex(
    aes(
      fill = after_stat(value)
    ),
    bins = 35,
    fun = mean,
    na.rm = TRUE
  ) +
  
  AGC_scale +
  
  labs(
    x = "Mean annual precipitation (mm)",
    y = "Climatic water deficit (mm)",
    title = "(d)"
  ) +
  
  theme_climate +
  
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = grid::unit(8, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  )
P_MAP_CWD

library(patchwork)
 
P_climate_hex <- (
    P_MAT_MAP |
      P_MAT_CWD |
      P_MAP_CWD
  ) +
  
  plot_layout(
    guides = "collect",
    widths = c(1, 1, 1)
  ) &
  
  theme(
    legend.position = "bottom"
  )


  
P_climate_hex


ggsave(
    filename = "Figure.png",
    plot = P_climate_hex,
    path = "/figure",
    width = 10,
    height = 4.5,
    units = "in",
    dpi = 600
  )

