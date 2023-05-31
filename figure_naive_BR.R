library(tidyverse)
library(tikzDevice)

## Illustrate limitations of naive simulation approach
xmax <- 100
dx <- 1
N <- 50

set.seed(2)

Gamma <- rexp(N) %>% cumsum()
df <- expand.grid(
  x = seq(0,100,dx),
  sim = as_factor(1:N)
) %>% 
  as_tibble() %>% 
  mutate(Gamma = Gamma[sim]) %>% 
  group_by(sim) %>% 
  mutate(
    W = c(0, rnorm(xmax/dx, sd = sqrt(dx))) %>% cumsum(),
    V = exp(W - x/2),
    zV = V / Gamma
  ) %>% 
  select(sim, x, zV)

p2 <- ggplot() +
  geom_line(
    aes(x, zV, group = sim),
    data = filter(df),
    col = alpha("black", 0.2)
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab(NULL) + xlab("$x$")

p2

tikz('../tex/fig/threshold_stopping_naive.tex',width=5.5,height=3)
p2
dev.off()