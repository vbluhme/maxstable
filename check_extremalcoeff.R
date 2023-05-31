source("sim_extremalfuns_S3.R")
library(tikzDevice)

extcoeff_BR <- function(h, s, alpha) {
  gamma <- (h/s)^alpha
  2 * pnorm(sqrt(gamma)/2)
}

## Plot Brown-Resnick extremal coefficient function ----
p1 <- expand.grid(h = seq(0,5,0.05), s = 1, alpha = c(0.5, 1, 1.5, 2)) %>% 
  as_tibble() %>% 
  mutate(
    theta = extcoeff_BR(h, s, alpha),
  ) %>% 
  ggplot(aes(h, theta, col = as.factor(-alpha))) +
  geom_line() +
  geom_text(
    aes(label = ifelse(h == 5, paste0("$\\alpha = ", alpha, "$"),NA)),
    hjust = 1, nudge_y = 0.05, size = 3, nudge_x = 0.1
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0,5.1) + 
  scale_y_continuous(limits = c(1,2.1), breaks = seq(1,2,0.25)) +
  xlab("$h$") + ylab("$\\theta(h)$") 


p2 <- expand.grid(h = seq(0,5,0.05), s = c(0.5, 1, 2, 5), alpha = 1) %>% 
  as_tibble() %>% 
  mutate(
    theta = extcoeff_BR(h, s, alpha),
  ) %>% 
  ggplot(aes(h, theta, col = as.factor(s))) +
  geom_line() +
  geom_text(
    aes(label = ifelse(h == 5, paste0("$s = ", s, "$"),NA)),
    hjust = 1, nudge_y = 0.05, size = 3
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(0,5) + 
  scale_y_continuous(limits = c(1,2.1), breaks = seq(1,2,0.25)) +
  xlab("$h$") + ylab(NULL) 

cowplot::plot_grid(p1,p2, nrow = 1)

tikz('../tex/fig/BR_extcoeff.tex',width=5,height=2.5)
cowplot::plot_grid(p1,p2, nrow = 1)
dev.off()


## Show convergence from simulations
set.seed(1)

n <- 20
coord <- data.frame(x = runif(n), y = runif(n)) %>% as.matrix()

alpha <- 1.2
s <- 0.2

param <- list(alpha = alpha, s = s)

BR_model <- BR_factory(coord, alpha = alpha, s = s)

simulate_and_calculate_theta <- function(k, sim_algorithm, algorithm_label) {
  cat(algorithm_label, "\tk =", k, "\n")
  Z <- replicate(k, sim_algorithm(BR_model, quietly = TRUE))
  SpatialExtremes::fitextcoeff(t(Z), coord, marge = "frech", plot = FALSE)$ext.coeff %>% 
    as_tibble() %>% 
    mutate(k = k, algorithm = algorithm_label)
}


k <- c(50, 500, 5000)
theta_sim <- bind_rows(
  map(k, simulate_and_calculate_theta, sim_algorithm = sim_extremalfuns, algorithm_label = "Extremal functions") %>% bind_rows(),
  map(k, simulate_and_calculate_theta, sim_algorithm = sim_sumnorm, algorithm_label = "Sum-normalization") %>% bind_rows(),
)


p <- theta_sim %>% 
  mutate(
    k_lab = paste0("k = ", k) %>% fct_reorder(k),
    k = as_factor(k)
  ) %>% 
  ggplot(aes(distance, ext.coeff, col = algorithm)) +
  geom_point(alpha = 0.25) +
  geom_function(fun = extcoeff_BR, args = param, n = 30) +
  facet_grid(algorithm ~ k_lab) +
  #theme_minimal() +
  theme(legend.position = "none", panel.background = element_rect(fill = "white"), panel.grid = element_line(colour = rgb(0,0,0,0.1))) +
  xlab("$h$") + ylab("$\\theta(h)$") +
  scale_y_continuous(limits = c(1,2), breaks = c(1,1.5,2)) +
  scale_x_continuous(limits = c(0,1.5))

p

tikz('../tex/fig/BR_extcoeff_convergence.tex',width=5.5,height=4)
p
dev.off()


## Marginal distributions
set.seed(1)


n <- 10
coord <- data.frame(x = runif(n), y = runif(n)) %>% as.matrix()

alpha <- 1.2
s <- 0.2

BR_model <- BR_factory(coord, alpha = alpha, s = s)
Smith_model <- smith_factory(type = "smith_sheet", sigma = s)


k <- 10000
df <- bind_rows(
  replicate(k, sim_extremalfuns(BR_model, coord, quietly = TRUE)) %>% 
    as_tibble() %>% 
    bind_cols(as_tibble(coord)) %>% 
    pivot_longer(-c(x,y)) %>% 
    mutate(
      name = parse_number(name),
      model = "Brown-Resnick",
      algorithm = "Extremal functions"
    ),
  
  replicate(k, sim_sumnorm(BR_model, coord, quietly = TRUE)) %>% 
    as_tibble() %>% 
    bind_cols(as_tibble(coord)) %>% 
    pivot_longer(-c(x,y)) %>% 
    mutate(
      name = parse_number(name),
      model = "Brown-Resnick",
      algorithm = "Sum-normalization"
    ),
  
  replicate(k, sim_extremalfuns(Smith_model, coord, quietly = TRUE)) %>% 
    as_tibble() %>% 
    bind_cols(as_tibble(coord)) %>% 
    pivot_longer(-c(x,y)) %>% 
    mutate(
      name = parse_number(name),
      model = "Smith",
      algorithm = "Extremal functions"
    ),
  
  replicate(k, sim_sumnorm(Smith_model, coord, quietly = TRUE)) %>% 
    as_tibble() %>% 
    bind_cols(as_tibble(coord)) %>% 
    pivot_longer(-c(x,y)) %>% 
    mutate(
      name = parse_number(name),
      model = "Smith",
      algorithm = "Sum-normalization"
    )
)

qfrechet <- function(p) 1/log(1/p)

p <- df %>% 
  rename(sim = name) %>% 
  group_by(sim, model, algorithm) %>% 
  filter(row_number() == 5) %>% 
  group_by(model, algorithm) %>% 
  summarise(
    p = seq(0.01, 0.99, by = 0.01),
    qfrechet_emp = quantile(value, p),
    qfrechet = qfrechet(p)
  ) %>% 
  ggplot(aes(qfrechet_emp, qfrechet, col = interaction(model, algorithm))) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  geom_point() +
  facet_grid(model ~ algorithm) +
  scale_x_log10() + scale_y_log10() +
  xlab("Empirical quantiles") + ylab("Fréchet quantiles") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_rect(colour = NA, fill="lightgray")
  )

p

tikz('../tex/fig/check_simulation_qqplot.tex',width=5.5,height=4)
p
dev.off()
