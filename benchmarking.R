### THIS FILE INCLUDES:
### 1. Benchmark for evaluation of fBs covariance matrix
### 2. Benchmark of simulation algorithms for Brown-Resnick and Smith models

source("sim_extremalfuns_S3.R")
library(tikzDevice)

# To make \mathbb work in tikz().
options(tikzLatexPackages =c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))

## Evaluation of fractional Brownian motion covariance matrix ----
alpha <- 1
s <- 1

## R function to set up covariance matrix
## This is simply chol_fBs() without chol()
f <- function(coord, alpha, s) {
  coord_norm <- apply(coord, 1, norm, type = "2")^alpha
  n <- nrow(coord)
  scale <- 1/(2 * s^alpha)
  f <- function(i,j) apply(coord[i,] - coord[j,], 1, norm, type = "2")
  
  outer(1:n, 1:n, function(i,j) scale * (coord_norm[i] + coord_norm[j] - f(i,j)^alpha))
}

bench_df_cov <- bench::press(
  n = c(10, 20, 50, 100, 200, 500, 1000),
  {
    x <- data.frame(x = runif(n), y = runif(n)) %>% as.matrix()
    bench::mark(
      min_iterations = 10,
      "R" = {f(x, alpha, s)},
      "Rcpp" = {x %>% fBs_cov(alpha, s)}
    )
  }
) %>% bind_rows(bench::press(
  n = c(2000, 5000, 10000),
  {
    x <- data.frame(x = runif(n), y = runif(n)) %>% as.matrix()
    bench::mark(
      min_iterations = 10,
      #"R" = {f(x, alpha, s)},
      "Rcpp" = {x %>% fBs_cov(alpha, s)}
    )
  }
))

bench_df_cov_median <- bench_df_cov %>% 
  unnest(c(time, gc)) %>% 
  group_by(n, expression) %>% 
  summarise(median = median(time))

## We want to get rid of "1.67m" for 100 seconds
# Remove "60" and "m" to jump from seconds to hours.
fixInNamespace("time_units", "bench")

p <- bench_df_cov %>% 
  unnest(c(time, gc)) %>% 
  group_by(expression, n) %>% 
  slice_sample(n = 50) %>% 
  mutate(median = median(time)) %>% 
  ggplot(aes(n, time, col = expression)) +
  geom_jitter(height = 0, width = 0.04, alpha = 0.2) + 
  geom_line(aes(n, median, col = expression), data = bench_df_cov_median) +
  scale_x_log10(labels = scales::comma) +
  bench::scale_y_bench_time(breaks = 10^(-5:2)) +
  theme_minimal() +
  labs(colour = NULL) +
  ylab("Time elapsed") +
  xlab("$n$")

p

# Save plot
tikz('../tex/fig/bench_BR_covmat.tex',width=5.5,height=3.5)
p
dev.off()

## NEXT: Benchmark simulations of Brown-Resnick and Smith models ----
## using both algorithms
alpha <- 1
s <- 1

extremal_tracer <- tracer(
  c("N", "n", "nsim"),
  N = 0
)

bench_df1 <- bench::press(
  n = c(10, 20, 50, 100, 200, 500, 1000, 2000),
  {
    x <- data.frame(x = runif(n), y = runif(n)) %>% as.matrix()
    smith_model <- smith_factory(type = "smith_sheet", sigma = 1)
    BR_model <- BR_factory(x, alpha, s)
    
    bench::mark(
      memory = FALSE, check = FALSE, min_iterations = 100,
      "Brown-Resnick model, Extremal functions" = {sim_extremalfuns(BR_model, quietly = TRUE, cb_nsim = extremal_tracer$tracer)},
      "Brown-Resnick model, Spectral measure" = {sim_sumnorm(BR_model, cb_nsim = extremal_tracer$tracer)},
      "Smith model, Extremal functions" = {sim_extremalfuns(smith_model, x, quietly = TRUE, cb_nsim = extremal_tracer$tracer)},
      "Smith model, Spectral measure" = {sim_sumnorm(smith_model, x, cb_nsim = extremal_tracer$tracer)}
    )
  }
)

bench_df2 <- bench::press(
  n = c(5000, 10000),
  {
    x <- data.frame(x = runif(n), y = runif(n)) %>% as.matrix()
    smith_model <- smith_factory(type = "smith_sheet", sigma = 1)
    
    bench::mark(
      memory = FALSE, check = FALSE, min_iterations = 100,
      "Smith model, Extremal functions" = {sim_extremalfuns(smith_model, x, quietly = TRUE, cb_nsim = extremal_tracer$tracer)},
      "Smith model, Spectral measure" = {sim_sumnorm(smith_model, x, cb_nsim = extremal_tracer$tracer)}
    )
  }
)

bench_df <- bind_rows(bench_df1, bench_df2) %>% 
  mutate(expression = as.character(expression), tmp = expression) %>% 
  separate(tmp, c("model", "algorithm"), sep = ", ") %>% 
  unnest(c(time, gc)) %>% 
  mutate(
    nsim = summary(extremal_tracer)$nsim,
    time_nsim = time / nsim
  ) %>% 
  select(expression, model, algorithm, n, time, nsim, time_nsim)

bench_df_medians <- bench_df %>% 
  group_by(model, expression, algorithm, n) %>% 
  summarise(
    median_time = median(time),
    median_nsim = median(nsim),
    median_time_nsim = median(time_nsim)
  )

save(bench_df, bench_df_medians, file = "bench_df.RData")

load("bench_df.RData")


## We want to get rid of "1.67m" for 100 seconds
# Remove "60" and "m" to jump from seconds to hours.
fixInNamespace("time_units", "bench")

options(scipen = 1000)

p1 <- bench_df %>% 
  group_by(n, expression) %>%
  slice_head(n = 50) %>% 
  ggplot(aes(n, time, col = algorithm, linetype = model)) +
  geom_jitter(height = 0, width = 0.04, alpha = 0.07) + 
  geom_line(aes(n, median_time), data = bench_df_medians) +
  geom_point(aes(n, median_time), data = bench_df_medians) +
  scale_x_log10(labels = scales::comma) +
  bench::scale_y_bench_time(
    labels = scales::comma,
    breaks = 10^(-5:2)
  ) +
  theme_minimal() +
  ylab("Time elapsed") +
  xlab(NULL) +
  theme(legend.position = "none") +
  labs(linetype = "Model", col = "Algorithm")

p2 <- bench_df %>% 
  group_by(n, expression) %>%
  slice_head(n = 50) %>% 
  ggplot(aes(n, nsim, col = algorithm, linetype = model)) +
  geom_jitter(height = 0, width = 0.04, alpha = 0.07) +
  geom_line(aes(n, median_nsim), data = bench_df_medians) +
  geom_point(aes(n, median_nsim), data = bench_df_medians) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma, breaks = 10^(0:5)) +
  theme_minimal() +
  ylab("Samples from $\\mathbb P_{x_0}$") +
  xlab(NULL) +
  theme(legend.position = "none") +
  labs(linetype = "Model", col = "Algorithm")

p3 <- bench_df %>% 
  group_by(n, expression) %>% 
  slice_head(n = 50) %>% 
  ggplot(aes(n, time_nsim, col = algorithm, linetype = model)) +
  geom_jitter(height = 0, width = 0.05, alpha = 0.07) +
  geom_line(aes(n, median_time_nsim), data = bench_df_medians) +
  geom_point(aes(n, median_time_nsim), data = bench_df_medians) +
  scale_x_log10(labels = scales::comma) +
  bench::scale_y_bench_time() +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()) +
  labs(colour = "Algorithm", linetype = "Model") +
  ylab("Time per sample from $\\mathbb P_{x_0}$") +
  xlab("$n$")

cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1,1,1.5))

tikz('../tex/fig/bench_extremal_cowplot.tex',width=5.5,height=6.5)
cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1,1,1.4))
dev.off()


