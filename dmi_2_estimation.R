source("sim_extremalfuns_S3.R")

library(lubridate)
library(sf)
library(geojsonsf)
library(extRemes)
library(tailDepFun)

library(tikzDevice)

library(geosphere)

# Shapefile for Denmark
library(rnaturalearth)      # Provides ne_countries
library(rnaturalearthhires) # Hi-res shapefiles

# Denmark without Bornholm
dk_map <- ne_countries(country = "Denmark", scale = "large", returnclass = "sf") %>% 
  st_crop(xmin = 0, xmax = 14, ymin = 0, ymax = 180)

bbox <- st_bbox(dk_map)

# Load files:
# - dmi_stations
# - dmi (daily maximum temperatures, 2011-2022)
load("dmi.RData")

## Add cartesian coordinates in km for isotropic BR estimation.
dmi_stations <- dmi_stations %>%
  mutate(
    x_cart = map2_dbl(x,y, function(x,y) geosphere::distm(c(bbox$xmin, y), c(x, y), fun = distHaversine)/1000),
    y_cart = map2_dbl(x,y, function(x,y) geosphere::distm(c(x, bbox$ymin), c(x, y), fun = distHaversine)/1000)
  )

# Plot stations
ggplot() +
  geom_sf(data = dk_map, fill = "antiquewhite") +
  geom_point(aes(x,y), data = dmi_stations) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = alpha("black", 0))
  ) +
  xlab(NULL) + ylab(NULL)

# Plot maximal summer temperature by year
dmi %>% 
  group_by(stationId, year) %>% 
  summarise(temp = max(temp)) %>% 
  ungroup() %>% 
  left_join(dmi_stations) %>% 
  ggplot() +
  geom_sf(data = dk_map, fill = "antiquewhite") +
  geom_point(aes(x,y, col = temp)) +
  scale_color_gradient(low = "blue", high = "red") +
  facet_wrap(~ year) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = alpha("black", 0)),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  ) +
  xlab(NULL) + ylab(NULL)


#### Extreme value analysis ----
dmi_block <- dmi %>% 
  filter(
    # Filter for 1 June until 6 September each year.
    month(date) >= 6,
    make_date(year = 2000, month(date), day(date)) <= "2000-09-06"
  ) %>%
  mutate(
    summer_day = date - make_date(year = year(date), month = 6, day = 1),
    window_2week = 1 + as.numeric(summer_day) %/% 14,
    block = length(unique(window_2week))*(year-min(year)) + window_2week,
  ) %>% 
  group_by(stationId, block) %>%
  summarise(temp = max(temp), start_date = min(date)) %>% 
  ungroup() %>% 
  left_join(dmi_stations)

# tikz('../tex/fig/dmi_temp_2022july.tex',width=5,height=4)
# Figure post-processing:
# - insert 'fig/' in raster path
ggplot() +
  geom_sf(data = dk_map, fill = "antiquewhite") +
  ggrepel::geom_label_repel(
    aes(x,y, label = paste0(name, ", ", temp)),
    data = dmi_block %>% filter(temp == max(temp)),
    nudge_y = -0.2, alpha = 0.7, size = 3
  ) +
  geom_point(
    aes(x, y, col = temp),
    data = dmi_block %>% filter(start_date == "2022-07-13"),
    size = 2
  ) +
  scale_color_gradient(low = "blue", high = "red", limits = c(22, 38)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = alpha("black", 0)),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  ) +
  labs(col = "$^\\circ$C") +
  xlab(NULL) + ylab(NULL)

# dev.off()

### Do we see temporal dependence?
p2 <- dmi_block %>%
  group_by(stationId, name) %>% 
  nest() %>% 
  mutate(
    ac = map(data, ~ tibble(ac = acf(.x$temp, plot = FALSE)$acf %>% c(), lag = seq_along(ac)-1))
  ) %>% 
  unnest(ac) %>% 
  group_by(lag) %>% 
  mutate(mean = mean(ac)) %>% 
  ungroup() %>% 
  filter(lag >= 1, lag <= 15) %>% 
  ggplot(aes(x = factor(lag), y = ac, group = stationId)) +
  geom_hline(yintercept = c(1,-1)*1.96 / sqrt(84), linetype = "dashed", alpha = 0.7) +
  geom_line(alpha = 0.2) +
  geom_line(aes(factor(lag), mean), col = "red") +
  #geom_point(aes(factor(lag), mean), col = "red", size = 2) +
  theme_minimal() +
  xlab("Lag") + ylab("Autocorrelation coefficient")


tikz('../tex/fig/dmi_acf.tex',width=4,height=2.5)
p2
dev.off()

## Fitting marginal distributions
fit1 <- extRemes::fevd(
  temp, data = data.frame(dmi_block),
  location.fun = ~ x + y
)
fit1$results$par

fit2 <- extRemes::fevd(
  temp, data = data.frame(dmi_block),
  location.fun = ~ x + y,
  scale.fun = ~ x + y
)
fit2$results$par

## Likelihood ratio test at 0.05 level
# This is probably not meaningful as the likelihoods are misspecified.
lr.test(fit1, fit2)

# We can replicate the LR-test
lrtest <- function(fit1, fit2, alpha = 0.05) {
  df <- length(fit2$parnames)-length(fit1$parnames)
  lr.ratio <- -2 * (fit2$results$value - fit1$results$value)
  critical <- qchisq(1-alpha, df = df)
  p.value <- pchisq(lr.ratio, df = df, lower.tail = FALSE)
  
  structure(
    list(df = df, lr.ratio = lr.ratio, critical = critical, p.value = p.value, alpha = alpha),
    class = "lrtest"
  )
}
print.lrtest <- function(lrtest_object) {
  print(paste(
    "LR ratio =", round(lrtest_object$lr.ratio, 4),
    "| df =", lrtest_object$df,
    "| critical value =", round(lrtest_object$critical, 3),
    "| p-value =", round(lrtest_object$p.value, 8),
    "| alpha =", lrtest_object$alpha
  ))
}
lrtest(fit1,fit2)

## AIC
AIC1 <- 2 * length(fit1$results$par) + 2 * fit1$results$value
AIC2 <- 2 * length(fit2$results$par) + 2 * fit2$results$value
  
## Bootstrap for cluster-robust standard errors
boot_function <- function(df) {
  df <- df[sample(df$block, replace = TRUE),] %>% 
    unnest(data)
  
  list(
    fevd1 = extRemes::fevd(temp, data = data.frame(df), location.fun = ~ x + y)$results$par,
    fevd2 = extRemes::fevd(temp, data = data.frame(df), location.fun = ~ x + y, scale.fun = ~ x + y)$results$par,
    linear = lm(temp ~ x + y, data = df)$coefficients
  )
}

set.seed(2)
n_boot <- 1000

dmi_block_nest <- dmi_block %>% nest(data = -block)
boot <- replicate(n_boot, boot_function(dmi_block_nest)) %>% t()

se_fevd1 <- boot[,1] %>% do.call(rbind, .) %>% apply(2, sd)
se_fevd2 <- boot[,2] %>% do.call(rbind, .) %>% apply(2, sd)
se_lm <- boot[,3] %>% do.call(rbind, .) %>% apply(2, sd)

se_fevd1
se_fevd2
se_lm

rm(boot_function, boot, dmi_block_nest)

### Set up and export a nice table with estimation results.
{
  fake_data <- tibble(
    temp = rnorm(100), mu0 = rnorm(100), mu1 = rnorm(100), mu2 = rnorm(100), scale = rnorm(100), sigma1 = rnorm(100), sigma2 = rnorm(100), shape = rnorm(100)
  )
  fake_model <- lm(temp ~ mu0+mu1+mu2+scale+sigma1+sigma2+shape-1, data = fake_data)
  
  fit1_df <- tibble(
      var = c("Location, const.", "Location, lon.", "Location, lat.", "Scale", "Shape"),
      par = fit1$results$par, 
      se = summary(fit1)$se.theta,
      se_boot = se_fevd1
    ) %>% 
    mutate(
      t = par / se_boot,
      p = 2 * pnorm(abs(t), lower.tail = FALSE)
    ) %>% 
    add_row(.before = 5) %>%
    add_row(.before = 5)
  
  fit2_df <- tibble(
    var = c("Location, const.", "Location, lon.", "Location, lat.", "Scale, const.", "Scale, lon.", "Scale, lat.", "Shape"),
    par = fit2$results$par, 
    se = summary(fit2)$se.theta,
    se_boot = se_fevd2
  ) %>% 
    mutate(
      t = par / se_boot,
      p = 2 * pnorm(abs(t), lower.tail = FALSE)
    )
  
  fit_lm <- lm(temp ~ x + y, data = dmi_block)
  fit_lm_df <- summary(fit_lm)$coefficients %>% 
    `colnames<-`(c("par", "se", "t", "p")) %>% 
    as_tibble() %>% 
    mutate(
      se_boot = se_lm,
      t = par / se_boot,
      p = 2 * pnorm(abs(t), lower.tail = FALSE)
    ) %>% 
    rbind(NA) %>% rbind(NA) %>% rbind(NA) %>%  rbind(NA)
  
  texreg::texreg(
    list(fake_model, fake_model, fake_model),
    custom.header = list("Temperature (block-max)" = 1:3),
    custom.model.names = c("Linear regression", "GEV model (1)", "GEV model (2)"),
    override.coef = list(fit_lm_df$par, fit1_df$par, fit2_df$par),
    override.se = list(fit_lm_df$se_boot, fit1_df$se_boot, fit2_df$se_boot),
    override.pvalues = list(fit_lm_df$p, fit1_df$p, fit2_df$p),
    custom.coef.names = fit2_df$var,
    center = TRUE,
    include.rsquared = FALSE,
    include.adjrs = FALSE,
    include.nobs = FALSE,
    include.rmse = FALSE,
    custom.gof.rows = list(
      "Observations" = c("4,116", "4,116", "4,116"),
      "Res. standard error" = c(summary(fit_lm)$sigma, NA, NA),
      "Log-likelihood" = c("",format(c(-fit1$results$value, -fit2$results$value), big.mark = ",")),
      "AIC" = c("",format(c(2*length(fit1$results$par) + 2 * fit1$results$value, 2*length(fit2$results$par) + 2 * fit2$results$value), big.mark = ",")),
      "BIC" = c("",format(c(log(length(fit1$x))*length(fit1$results$par) + 2 * fit1$results$value, log(length(fit2$x))*length(fit2$results$par) + 2 * fit2$results$value), big.mark = ","))
    ),
    digits = 2,
    booktabs = TRUE,
    caption = "Parameter estimates for GEV maximum-likelihood estimation. Standard errors in parentheses obtained by bootstrapping from $1000$ samples with replacement from the $k = 84$ blocks. Linear regression included for comparison.",
    use.packages = FALSE,
    file = "../tex/fig/tab_GEV_estimates.tex"
  )
  
  rm(fit1_df, fit2_df, fit_lm_df, fake_model, fake_data, fit_lm)
}


## How does fevd work? We can replicate the results:
loglik <- function(theta) {
  mu <- theta[1] + theta[2] * dmi_block$x + theta[3] * dmi_block$y
  sig <- theta[4]; xi <- theta[5]
  z <- dmi_block$temp
  m <- length(z)

  y <- 1 + xi * (z-mu)/sig
  if(min(y) <= 0) return(10^6)
  m * log(sig) + (1+1/xi) * sum(log(y)) + sum(y^(-1/xi))
}

## Same optimum as fevd :)
init <- lm(temp~x + y, data = dmi_block)$coefficients %>% c(1,-1) %>% as.numeric()
fit1_optim <- optim(
  par = init,
  loglik,
  control = list(maxit = 100000, reltol = 1e-20),
  hessian = TRUE
)
rm(init)

tibble(
  parnames = fit1$parnames,
  fevd_par = fit1$results$par,
  fevd_se = fit1$results$hessian %>% solve() %>% diag() %>% sqrt(),
  optim_par = fit1_optim$par,
  optim_se = fit1_optim$hessian %>% solve() %>% diag() %>% sqrt(),
  boot_se = se_fevd1
)

## Plot mu(x) and sigma(x)
df_tmp <- expand.grid(
    x = seq(bbox$xmin, bbox$xmax, length.out = 500),
    y = seq(bbox$ymin, bbox$ymax, length.out = 500)
  ) %>% 
  mutate(
    mu = fit1$results$par[1] + fit1$results$par[2] * x + fit1$results$par[3] * y,
    sig = fit2$results$par[4] + fit2$results$par[5] * x + fit2$results$par[6] * y,
  ) %>% 
  mutate(x1 = x, y1 = y) %>% 
  df_geojson(lon = "x", lat = "y") %>% 
  geojson_sf() %>% 
  st_filter(dk_map) %>% 
  as_tibble() %>%
  rename(y = y1, x = x1)

p1 <- df_tmp %>% 
  ggplot() +
  geom_sf(data = dk_map, alpha = 0, col = alpha("black", 0)) +
  geom_tile(aes(x, y, fill = mu)) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = alpha("black", 0))
  ) +
  xlab(NULL) + ylab(NULL) +
  labs(fill="μ(x)")

p2 <- df_tmp %>% 
  ggplot() +
  geom_sf(data = dk_map, alpha = 0, col = alpha("black", 0)) +
  geom_tile(aes(x, y, fill = sig)) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = alpha("black", 0))
  ) +
  xlab(NULL) + ylab(NULL) +
  labs(fill="σ(x)")

p1
p2

##### Best result obtained by clicking Zoom, resizing and taking a screenshot. 
#### tikz creates a .tex file too big for the latex engine and cannot turn the map into a raster image.
# library(tikzDevice)
# #tikz('../tex/fig/dmi_mu.tex',width=4,height=4)
# pdf('../tex/fig/dmi_mu.pdf', width = 5, height = 5)
# p1
# dev.off()

rm(df_tmp, p1, p2)


### Back to our regular programming ....
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

dmi_block <- dmi_block %>% 
  mutate(
    Z = standardize_evd(temp, x, y, fit2$results$par)
  )

## QQ-plot, compare to Fréchet distribution
{
  #tikz(file = "../tex/fig/dmi_frechet", width = 6, height = 3, sanitize = TRUE)
  p <- par(mfrow = c(1,2), mar = c(4, 4, 0, 0) + 0.1)
  
  qfrechet <- function(p) 1/log(1/p)
  q <- seq(0.005, 0.995, by = 0.005)
  extRemes::qqplot(
    x = quantile(dmi_block$Z, q),
    y = qfrechet(q),
    regress = FALSE,
    log = "xy",
    xlab = "Empirical quantiles",
    ylab = "Fréchet quantiles"
  )
  
  break_fun <- function(x) seq(0, max(x)+1, by= 0.2)
  hist(
    dmi_block$Z,
    prob = TRUE,
    xlim = c(0,8),
    breaks = break_fun,
    main = NULL, xlab = "Z", ylab = NULL
  )
  curve(exp(-1/x)/x^2, add = TRUE, col = "red")
  par(p)
  #dev.off()
  
  rm(p, qfrechet, q, break_fun)
}

## Mean excess plot
## Linear plot -> Fréchet distribution
ME <- function(u, x) mean(x[x > u]- u)
u <- seq(0,200, length.out = 1000)
plot(
  u,
  map_dbl(u, ME, x = dmi_block$Z), 
  type = "l",
  xlab = "$u$",
  ylab = "$\\mathbb E[Z - u \\mid Z > u]$",
  main = "Mean excess plot of standardized observations"
)
rm(u, ME)

save(dmi_block, dmi_stations, dk_map, bbox, file = "dmi_block.RData")

#### FIT BROWN-RESNICK MODEL
load("dmi_block.RData")

set.seed(5)

n_pairs <- nrow(dmi_stations) * (nrow(dmi_stations)-1) / 2
keep_pairs1 <- sample(n_pairs, 100)
keep_pairs2 <- sample(n_pairs, 15)

indices <- selectGrid(cst = c(0,1), d = nrow(dmi_stations)) %>% 
  `colnames<-`(dmi_stations$stationId)

indices1 <- indices[keep_pairs1,]
indices2 <- indices[keep_pairs2,]

### Locations are cartesian coordinates,
#### km distance from corner of dk_map bbox
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

bind_rows(
  indices1 %>%
    as_tibble() %>%
    mutate(pair = keep_pairs1) %>%
    pivot_longer(-pair, names_to = "stationId") %>%
    filter(value == 1) %>% select(-value) %>% 
    mutate(indices = "Pairs for BR estimation with identity weight matrix"),
  indices2 %>%
    as_tibble() %>%
    mutate(pair = keep_pairs2) %>%
    pivot_longer(-pair, names_to = "stationId") %>%
    filter(value == 1) %>% select(-value) %>% 
    mutate(indices = "Pairs for BR estimation with optimal weight matrix")
  ) %>% 
  left_join(dmi_stations, by = "stationId") %>%
  ggplot() +
  geom_sf(data = dk_map, fill = "antiquewhite") +
  geom_point(aes(x,y), data = dmi_stations) +
  geom_line(aes(x,y,group = pair), alpha = 0.3) +
  theme_minimal() +
  facet_wrap(~ indices) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  ) +
  xlab(NULL) + ylab(NULL)

# M-estimator

### 100 pairs, no weight matrix iteration
BR_fit <- EstimationBR(
  x,
  locations_cart,
  indices1,
  k = 10,
  method = "Mestimator",
  startingValue = c(1, 200),
  isotropic = TRUE,
  iterate = FALSE, covMat = FALSE
)

### 15 pairs, optimal weight matrix + covariance matrix
BR_fit2 <- EstimationBR(
  x,
  locations_cart,
  indices2,
  k = 10,
  method = "Mestimator",
  startingValue = c(1, 200),
  isotropic = TRUE,
  iterate = TRUE, covMat = TRUE
)

## If fitting is too slow :)
# BR_fit <- list(theta = c(0.9, 9))

### Convert from (α, ρ) to (α, s)
### EstimationBR uses semivariogram parametrization.
BR_fit$alpha <- BR_fit$theta[1]
BR_fit$s <- BR_fit$theta[2] * 2^(1/BR_fit$alpha)

BR_fit2$alpha <- BR_fit2$theta[1]
BR_fit2$s <- BR_fit2$theta[2] * 2^(1/BR_fit2$alpha)

BR_bootstrap <- function(x) {
  x <- x[sample(1:nrow(x), replace = TRUE),]
  
  keep_pairs_boot <- sample(n_pairs, 100)
  indices_boot <- indices[keep_pairs_boot]
  
  BR_fit <- EstimationBR(
    x,
    locations_cart,
    indices_boot,
    k = 10,
    method = "Mestimator",
    startingValue = c(1, 200),
    isotropic = TRUE,
    iterate = FALSE, covMat = FALSE
  )
  BR_fit$theta[2] <- BR_fit$theta[2] * 2^(1/BR_fit$theta[1])
  BR_fit$theta
}
BR_fit1_boot <- replicate(20, BR_bootstrap(x))
BR_fit$se <- apply(BR_fit1_boot, 1, sd)

BR_fit2$se <- sqrt(diag(BR_fit2$covMatrix)) * c(1, 2^(1/BR_fit2$theta[1]))

save(fit1, fit2, BR_fit, BR_fit2, file = "BR_fit.RData")

## Let's create the LaTeX table manually :) :) :)

BR_table_function <- function(fit1, fit2) {
  paste0(
"\\begin{table}
\\begin{center}
\\begin{tabular}{l c c}
\\toprule
 & BR model 1 & BR model 2                                                      \\\\
\\midrule
$\\alpha$    		& $",round(fit1$alpha, 3),"$      & $",round(fit2$alpha,3),"$   \\\\
                & $(",round(fit1$se[1], 3),")$    & $(",round(fit2$se[1],3),")$ \\\\
$s$      			  & $",round(fit1$s, 1),"$          & $",round(fit2$s,1),"$       \\\\
                & $(",round(fit1$se[2], 1),")$    & $(",round(fit2$se[2],1),")$ \\\\
\\midrule
Blocks          & 84                              & 84                          \\\\
Locations       & 49                              & 49                          \\\\
Pairs           & 100                             & 15                          \\\\
Weight matrix   & Identity                        & Iterative                   \\\\
Standard errors & Bootstrap                       & Analytical                  \\\\
\\bottomrule
\\end{tabular}
\\caption{Brown-Resnick parameter estimation employing the M-estimator for spatial tail dependence. Standard errors (in parentheses) are obtained through bootstrapping in Model 1 and from the estimation routine in Model 2.}
\\label{table:BR_estimation}
\\end{center}
\\end{table}"
  ) %>% print()
}
BR_table_function(BR_fit, BR_fit2) %>% cat(file = "../tex/fig/tab_BR_estimates.tex")

file.show("../tex/fig/tab_BR_estimates.tex")

### Let's do the plot
theta <- SpatialExtremes::fitextcoeff(x, locations_cart, marge = "frech", plot = FALSE)$ext.coeff
p <- theta %>% 
  as_tibble() %>% 
  ggplot(aes(distance, ext.coeff)) +
  geom_point(alpha = 0.2) +
  geom_function(fun = extcoeff_BR, args = list(alpha = BR_fit$alpha, s = BR_fit$s), size = 1, aes(col = "Model 1")) +
  geom_function(fun = extcoeff_BR, args = list(alpha = BR_fit2$alpha, s = BR_fit2$s),size = 1, aes(col = "Model 2")) +
  #geom_line(data = df_extcoeff, aes(distance, ext.coeff, group = model, col = model), size = 2) +
  theme_minimal() +
  labs(col = NULL)+
  xlab("Distance (km)") + ylab("Extremal coefficient")
  
p

# tikz('../tex/fig/dmi_extcoeff_compare_BR.tex',width=5,height=3)
# p
# dev.off()

