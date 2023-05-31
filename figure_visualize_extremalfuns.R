source("sim_extremalfuns_S3.R")

#set.seed(55223)
set.seed(446617)

x <- 1:4 
#x <- seq(0.1,1.2,0.1)
extremal_tracer <- tracer(
  c("k", "zY", "Z", "Gammainv", "accept"),
  expr = quote({Gammainv <- 1/Gamma})
)

BR_model <- BR_factory_univariate(x, alpha = 0.2, s = 1)

Z <- sim_extremalfuns(BR_model, cb = extremal_tracer$tracer)

df <- summary(extremal_tracer) %>% tibble()

## Adjustments for nicer plots
## These plots are for illustrative purposes only...
df <- df[-3,]
df[5, "zY.4"] <- 3.2


df <- df %>% 
  mutate(
    Y_index = cumsum(!is.na(zY.1)),
    Y_index = replace(Y_index, is.na(zY.1), NA)
  ) %>% 
  pivot_longer(
    -c(Y_index, k, Gammainv, accept, .time),
    names_to = c(".value", "x_index"),
    names_sep = "[.]"
  ) %>% 
  mutate(
    x_index = as.numeric(x_index),
    accept = as.logical(accept),
    #Y_index = as.numeric(Y_index),
    x = x[x_index],
    xk = x[k],
    nlab = paste0("$k = ", k, "$")
  )

p1 <- df %>% 
  #group_by(n) %>% filter(max(accept) == 1) %>% ungroup() %>% 
  ggplot(aes(x, zY, group = Y_index, col = accept, label = Y_index)) + 
  geom_line() +
  geom_point() +
  scale_color_manual(values=c(rgb(0,0,0,0.25), rgb(1,0,0,1))) +
  theme_minimal() +
  facet_wrap(~ nlab, nrow = 1) +
  geom_vline(aes(xintercept = xk), col = rgb(0,1,0,0.5)) +
  geom_point(aes(x = xk, y = Gammainv), col = rgb(0.5,1,0.5)) +
  geom_line(aes(x, Z), data = filter(df, is.na(Gammainv)), col = "black", linetype = "dashed") +
  #geom_line(aes(x, Z), data = filter(df, is.na(Gammainv)), col = "gray", linetype = "dashed") +
  #geom_line(aes(x, Z), data = filter(df, is.na(accept) & !is.na(Gammainv)), col = "black", linetype = "dashed") +
  geom_label_repel(data = filter(df, x_index == 1)) +
  scale_x_continuous(breaks = x, labels = ~ paste0("$x_", .x, "$")) +
  theme(legend.position = "none") +
  ylab("") + xlab("") +
  ylim(c(0,3.5))

p1

library(tikzDevice)
tikz('../tex/fig/extremal_functions_illustration.tex',width=6,height=2.5)
p1
dev.off()

stop("Stop en halv!")