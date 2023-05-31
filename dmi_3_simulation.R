source("sim_maxstable.R")

library(ggrepel)

library(tidyverse)
library(tikzDevice)

library(extRemes)
library(tailDepFun)

library(sf)
library(geojsonsf)
library(geosphere)

# Shapefile for Denmark
library(rnaturalearth)      # Provides ne_countries
library(rnaturalearthhires) # Hi-res shapefiles

load("BR_fit.RData")
load("dmi_block.RData")

standardize_evd <- function(value, x, y, par) {
  mu <- par["mu0"] + par["mu1"] * x + par["mu2"] * y
  
  if ("scale" %in% names(par)) {
    sig <- par["scale"]
  } else {
    sig <- par["sigma0"] + par["sigma1"] * x + par["sigma2"] * y
  }
  xi <- par["shape"]
  
  (1 + xi * (value - mu) / sig) ^ (1/xi)
}

inverse_standardize_evd <- function(z, x, y, par) {
  mu <- par["mu0"] + par["mu1"] * x + par["mu2"] * y
  if ("scale" %in% names(par)) {
    sig <- par["scale"]
  } else {
    sig <- par["sigma0"] + par["sigma1"] * x + par["sigma2"] * y
  }
  xi <- par["shape"]
  
  (z^xi - 1) / xi * sig + mu
}

### Simulate from estimated model ----
## Seed 3, length.out 30 for picture


## A) Simulate on dense grid for illustrations
set.seed(2)

coord <- expand.grid(
  x = seq(bbox$xmin, bbox$xmax, length.out = 60)+0.001,
  y = seq(bbox$ymin, bbox$ymax, length.out = 60)+0.001
) %>% 
  as.matrix()

coord_cart <- coord %>% 
  as_tibble() %>% 
  mutate(
    x_cart = map2_dbl(x,y, function(x,y) geosphere::distm(c(bbox$xmin, y), c(x, y), fun = distHaversine)/1000),
    y_cart = map2_dbl(x,y, function(x,y) geosphere::distm(c(x, bbox$ymin), c(x, y), fun = distHaversine)/1000)
  ) %>% 
  select(x_cart, y_cart) %>% 
  as.matrix() 

BR_model <- BR_factory(coord_cart, alpha = BR_fit$alpha, s = BR_fit$s)
Z <- replicate(
  9,
  sim_extremalfuns(BR_model, coord_cart, quietly = FALSE)
)
colnames(Z) <- paste0("sim", 1:ncol(Z))

coord %>% 
  as_tibble() %>% 
  bind_cols(as_tibble(Z_figure)) %>% 
  pivot_longer(-c(x,y), values_to = "Z", names_to = "sim") %>% 
  mutate(
    temp = inverse_standardize_evd(Z, x, y, fit2$results$par)
  ) %>%
  select(-Z) %>%
  mutate(sim = parse_number(sim)) %>% 
  filter(sim <= 9) %>% 
  group_by(sim) %>%
  nest() %>%
  mutate(data = purrr::map(data, interpolate_sheet, 200)) %>% 
  unnest(data) %>% 
  mutate(x1 = x, y1 = y) %>% 
  df_geojson(lon = "x", lat = "y") %>% 
  geojson_sf() %>% 
  st_filter(dk_map) %>% 
  as_tibble() %>%
  rename(y = y1, x = x1) %>% 
  group_by(sim) %>% 
  mutate(
    maxtemp = max(temp),
    maxtemp = if_else(temp == maxtemp, maxtemp, NA_real_)
  ) %>% 
  ggplot() +
  geom_sf(data = dk_map, alpha = 0, col = alpha("black", 0)) +
  geom_tile(aes(x, y, fill = temp)) +
  scale_fill_gradient(low = "blue", high = "red")+#, limits = c(15, 35))+
  geom_point(aes(x,y), data = function(x) subset(x, !is.na(maxtemp))) +
  geom_label_repel(aes(x,y, label = round(maxtemp, 1)), data = function(x) subset(x, !is.na(maxtemp)), alpha = 0.7) +
  facet_wrap(~ sim) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = alpha("black", 0)),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  ) +
  labs(fill = "°C") +
  xlab(NULL) + ylab(NULL)


## B) Simulate at station locations 
set.seed(1)

locations <- dmi_stations %>% 
  column_to_rownames("stationId") %>% 
  select(x,y) %>% 
  as.matrix()

locations_cart <- dmi_stations %>% 
  column_to_rownames("stationId") %>% 
  select(x_cart,y_cart) %>% 
  as.matrix()

x <- dmi_block %>% 
  arrange(block, stationId) %>% 
  group_by(block) %>% 
  pivot_wider(block, names_from = stationId, values_from = Z) %>% 
  ungroup() %>% 
  select(-block) %>% 
  as.matrix()

BR_model <- BR_factory(
  locations_cart,
  alpha = BR_fit$alpha,
  s = BR_fit$s
)

k <- 350000
Z <- replicate(
  k, 
  sim_extremalfuns(BR_model, locations_cart, quietly = TRUE)
)

## Extremal coefficients for k = 84 sample
theta <- SpatialExtremes::fitextcoeff(x, locations_cart, marge = "frech", plot = FALSE)$ext.coeff
theta_sim <- SpatialExtremes::fitextcoeff(t(Z[,1:84]), locations_cart, marge = "frech", plot = FALSE)$ext.coeff

## Compare to extremal coefficient from data
p <- left_join(
  as_tibble(theta),
  as_tibble(theta_sim), 
  by = "distance"
) %>% 
  rename(
    Data = ext.coeff.x,
    Simulation = ext.coeff.y
  ) %>% 
  pivot_longer(-distance) %>% 
  ggplot(aes(distance, value, col = name, group = name)) +
  geom_point(alpha = 0.8) +
  geom_function(fun = extcoeff_BR, args = list(alpha = BR_fit$alpha, s = BR_fit$s), col = "black") +
  facet_wrap(~ name) +
  theme_minimal() +
  xlab("Distance (km)") + ylab("Extremal coefficient") +
  theme(
    legend.position = "none",
    strip.background = element_rect(colour = NA, fill="lightgray")
  )

p

tikz('../tex/fig/dmi_sim_extcoeff.tex',width=6,height=3.5)
p
dev.off()

dmi_simulation <- dmi_stations %>% 
  bind_cols(as_tibble(Z)) %>% 
  pivot_longer(-c(stationId, name, x,y, x_cart, y_cart), values_to = "Z", names_to = "sim") %>% 
  mutate(
    temp = inverse_standardize_evd(Z, x, y, fit2$results$par)
  ) %>% 
  select(stationId, name, sim, Z, temp) %>% 
  mutate(sim = parse_number(sim)) 

save(dmi_simulation, file = "dmi_simulation.RData")

load("dmi_simulation.RData")

# Probability of exceeding `heatwave_limit` degrees
heatwave_limit <- 35

dmi_simulation %>% 
  group_by(sim) %>% 
  filter(temp[name == "Landbohøjskolen"] >= heatwave_limit) %>% 
  group_by(stationId) %>% 
  summarise(
    Probability = mean(temp >= heatwave_limit)
  ) %>% 
  left_join(dmi_stations, by = "stationId") %>% 
  ggplot() +
  geom_sf(data = dk_map, fill = "antiquewhite") +
  geom_point(aes(x,y, col = Probability)) +
  geom_label_repel(aes(x,y,label = round(Probability, 2), col = Probability)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = alpha("black", 0)),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  xlab(NULL) + ylab(NULL)

## Return period (years)
return_levels <- dmi_simulation %>% 
  group_by(sim) %>% 
  summarise(
    Copenhagen = temp[name == "Landbohøjskolen"],
    Esbjerg = temp[name == "Esbjerg Lufthavn"],
    Aarhus = temp[name == "Århus Syd"],
    Thyborøn = temp[name == "Thyborøn"],
    Denmark = max(temp)
  ) %>% 
  pivot_longer(-sim) %>% 
  group_by(name) %>% 
  summarise(
    Years = c(1, 10, 100, 1000),
    return_level = quantile(value, 1-1/Years/7)
  ) %>% 
  pivot_wider(values_from = return_level, names_from = name) %>% 
  relocate(Denmark, .after = last_col())

p <- return_levels %>% 
  pivot_longer(-Years) %>% 
  mutate(name = as_factor(name) %>% fct_relevel("Denmark", after = 0)) %>% 
  ggplot(aes(Years, value, col = name)) +
  geom_line() + geom_point() +
  scale_x_log10() +
  theme_minimal() +
  ylab("Return level ($^\\circ$C)") +
  labs(col = NULL)

p

tikz('../tex/fig/dmi_sim_returnlevel.tex',width=5,height=2.8)
p
dev.off()