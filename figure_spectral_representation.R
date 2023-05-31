library(tidyverse)
library(tikzDevice)

## Linear interpolation for df with $x and $zV
interpolate_approx <- function(df, npoints = 100) {
  xx <- seq(min(df$x), max(df$x), length.out = npoints)
  
  approx(x = df$x, y = df$zV, xout = xx) %>%
    as_tibble() %>% 
    rename(zV = y)
}


### Visualisation of spectral representation
### We simulate on 100 grid points
### The pointwise max is calculated on 500 grid points following linear interpolation.
### Just to make the plot nicer :)

N <- 50
dx <- 0.01

set.seed(9)

Gamma <- rexp(N) %>% cumsum()

df <- expand.grid(
  x = seq(0,1,dx),
  sim = as_factor(1:N)
) %>% 
  as_tibble() %>% 
  mutate(Gamma = Gamma[sim]) %>% 
  group_by(sim) %>% 
  mutate(
    W = c(0, rnorm(1/dx, sd = sqrt(dx))) %>% cumsum(),
    V = exp(W - x/2),
    zV = V / Gamma
  ) %>% 
  select(sim, x, zV)

# Get pointwise maximum. 
# Interpolate linearly before taking pointwise max for nice plot
df_max <- df %>% 
  nest(data = -sim) %>% 
  mutate(data = map(data, interpolate_approx, npoints = 300)) %>% 
  unnest(data) %>% 
  group_by(x) %>% 
  summarise(
    Z = max(zV),
    sim_max = which.max(zV)
  )

## Index of extremal functions 
sim_max_unique <- df_max %>% pull(sim_max) %>% unique()

p1 <- ggplot() +
  geom_line(
    aes(x, zV, group = sim),
    data = filter(df, !(sim %in% sim_max_unique)),
    col = alpha("black", 0.2)
  ) +
  geom_line(
    aes(x, zV, group = sim, col = sim), 
    data = filter(df, sim %in% sim_max_unique)
  ) +
  geom_line(
    aes(x, Z),
    data = df_max,
    linetype = "dashed",
    alpha = 0.8
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab(NULL) + xlab("$x$")

p1

tikz('../tex/fig/spectral_representation.tex',width=5.5,height=3)
p1
dev.off()



## Illustrate limitations of naive approachs
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


rY <- function() {
  
}
sim_naive <- function(tau, xmax, n) {
  dx <- xmax/n
  x <- seq(0,xmax, dx)
  Gamma <- rexp(1); Z <- numeric(n+1)
  while(tau/Gamma >= min(Z)) {
    W <- c(0, rnorm(n, sd = sqrt(dx))) %>% cumsum()
    V <- exp(W - 0.5 * x)
    Z <- pmax(Z, V/Gamma)
    Gamma <- Gamma + rexp(1)
  }
  Z
}

xmax <- 1
n <- 100

dx <- xmax/n
x <- seq(0,xmax, dx)
df <- replicate(100, sim_naive(0.0001, 100, 100)) %>% 
  as_tibble() %>% 
  mutate(x = x) %>% 
  pivot_longer(-x, names_to = "sim") %>% 
  mutate(sim = as_factor(parse_number(sim)))

df %>%
  ggplot(aes(x, value, group = sim, col = sim)) +
  geom_line() +
  theme(legend.position = "none")

xmax <- 100
dx <- 1
N <- 10000

set.seed(2)

Gamma <- rexp(N) %>% cumsum()
expand.grid(
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
  select(sim, x, zV) %>% 
  group_by(x) %>% 
  summarise(
    zV = max(zV)
  ) %>% 
  ggplot(aes(x, zV)) +
  geom_line()