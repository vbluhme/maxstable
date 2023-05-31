source("sim_extremalfuns_S3.R")
library(tikzDevice)

set.seed(3)

coord <- expand.grid(
  x = seq(0.01, 1, 0.01),
  y = seq(0.01, 1, 0.01)
)

## Illustration of extremal functions simulation
recalculate_L = TRUE
if(recalculate_L) {
  BR_obj1 <- BR_factory(coord, alpha = 0.5, s = 1)
  BR_obj2 <- BR_factory(coord, alpha = 1, s = 1)
  BR_obj3 <- BR_factory(coord, alpha = 1, s = 5)
  BR_obj4 <- BR_factory(coord, alpha = 1.5, s = 1)

  smith_obj1 <- smith_factory(type = "smith_sheet", sigma = 0.1)
  smith_obj2 <- smith_factory(type = "smith_sheet", sigma = 0.5)
  
  save(
    coord, BR_obj1, BR_obj2, BR_obj3,
    BR_obj4, smith_obj1, smith_obj2,
    file = "BR_obj_grid100_v2.RData")
}
load("BR_obj_grid100_v2.RData")

# This takes a few hours
recalculate_Z <- TRUE
if(recalculate_Z) {
  Z1 <- sim_extremalfuns(BR_obj1)
  Z2 <- sim_extremalfuns(BR_obj2)
  Z3 <- sim_extremalfuns(BR_obj3)
  Z4 <- sim_extremalfuns(BR_obj4)
  
  Z5 <- sim_extremalfuns(smith_obj1, coord)
  Z6 <- sim_extremalfuns(smith_obj2, coord)
  
  save(
    Z1, Z2, Z3, Z4, Z5, Z6,
    file = "Z_extremal_funs_seed3.RData"
  )
}
load("Z_extremal_funs_seed3.RData")

plot_fun <- function(df, label = NULL) {
  df %>% 
    interpolate_sheet(500) %>% 
    ggplot(aes(x,y)) +
    geom_raster(aes(fill = log(Z))) +
    scale_fill_gradient(low = "#BFE1B0", high = "#137177") +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab(label) + ylab(NULL)
}

p1 <- as_tibble(coord) %>% 
  bind_cols(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3,
    Z4 = Z4
  ) %>% 
  pivot_longer(-c(x,y), names_to = "sim", values_to = "Z") %>% 
  nest(data = c(x,y,Z)) %>% 
  mutate(
    label = c("(a) $\\alpha = 0.5$, $s = 1$",
              "(b) $\\alpha = 1$, $s = 1$",
              "(c) $\\alpha = 1$, $s = 5$",
              "(d) $\\alpha = 1.5$, $s = 1$")
  ) %>% 
  mutate(p = map2(data, label, plot_fun)) %>% 
  pull(p) %>% 
  cowplot::plot_grid(plotlist = ., ncol = 2)

p2 <- as_tibble(coord) %>% 
  bind_cols(
    Z5 = Z5,
    Z6 = Z6
  ) %>% 
  pivot_longer(-c(x,y), names_to = "sim", values_to = "Z") %>% 
  nest(data = c(x,y,Z)) %>% 
  mutate(
    label = c("(a) $\\sigma = 0.1$",
              "(b) $\\sigma = 0.5$")
  ) %>% 
  mutate(p = map2(data, label, plot_fun)) %>% 
  pull(p) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 1)

p1
p2

tikz('../tex/fig/extremal_functions_sheet_BR.tex',width=6,height=6)
p1
dev.off()

tikz('../tex/fig/extremal_functions_sheet_Smith.tex',width=6,height=3)
p2
dev.off()
