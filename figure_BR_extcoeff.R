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