#Load Required Packages
library(deSolve)
library(ggplot2)
library(gridExtra)
library(AICcmodavg)



# Step 1: Load Your CSV File
data <- read.csv("CovidNewCase.csv")  
days <- data[[1]]                  
# Compute estimated daily infectious population (I_{S_t}) from real data
D <- 9  # 9 days for COVID-19
I_est <- numeric(nrow(data)) #real case count for I_S

# Fill first D-1 days with cumulative sum
for (t in 1:(D - 1)) {
  I_est[t] <- sum(data$NewCase[1:t])
}

# Fill rest with rolling sum of last D days
for (t in D:nrow(data)) {
  I_est[t] <- sum(data$NewCase[(t - D + 1):t])
}

# Store in data
data$I_est <- I_est
real_data <- I_est
#Step 2: Define Constants
mu <- 1 / (80.3 * 365)
mu_S <- 0.00159
p <- 0.956
eta <- 5.5
epsilon <- 1 / 2.381
omega <- 1 / (eta - 1 / epsilon)

#Step 3: SEIR Model Function
seir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + E + IA + IS + R
    dS <- mu * (N - S) - ((beta_A * IA + beta_S * IS) / N) * S
    dE <- ((beta_A * IA + beta_S * IS) / N) * S - (epsilon + mu) * E
    dIA <- epsilon * E - (omega + mu) * IA
    dIS <- (1 - p) * omega * IA - (nu + mu_S) * IS
    dR <- p * omega * IA + nu * IS - mu * R
    list(c(dS, dE, dIA, dIS, dR))
  })
}

#Step 4: Initial Conditions
I_S0 <- 1             # use first real case count
I_A0 <- I_S0                      # assume I_A0 = I_S0
R0_estimate <- 2.5
E0 <- R0_estimate * (I_A0 + I_S0)
S0 <- 331002647
init <- c(S = S0, E = E0, IA = I_A0, IS = I_S0, R = 0)

#Step 5: Time Vector
times <- seq(0, length(real_data), by = 1)

#Step 6: Define Parameters
params_A <- c(beta_A = 0.55, beta_S = 0.55, nu = 0.05,
              epsilon = epsilon, omega = omega, mu = mu,
              mu_S = mu_S, p = p)

params_B <- c(beta_A = 0.5, beta_S = 0.45, nu = 0.0305,
              epsilon = epsilon, omega = omega, mu = mu,
              mu_S = mu_S, p = p)

#Step 7: Run Simulations
out_A <- ode(y = init, times = times, func = seir_model, parms = params_A)
out_B <- ode(y = init, times = times, func = seir_model, parms = params_B)

df_A <- as.data.frame(out_A)
df_B <- as.data.frame(out_B)

#Step 8: Model Predictions for Comparison
pred_A <- df_A$IS  
pred_B <- df_B$IS

# Trim to match real data length
min_len <- min(length(pred_A), length(real_data))
real_data <- real_data[1:min_len]
pred_A <- pred_A[1:min_len]
pred_B <- pred_B[1:min_len]
days <- days[1:min_len]

#Step 9: Fit Models and Compute AIC
lm_A <- lm(real_data ~ pred_A)
lm_B <- lm(real_data ~ pred_B)
aic_A <- AIC(lm_A)
aic_B <- AIC(lm_B)


y_max <- max(c(real_data, pred_A, pred_B), na.rm = TRUE)
y_min <- 0  # usually start infectious plots at 0

common_ylim <- c(y_min, y_max)

# === Define label x and individual label y for each plot ===
label_x <- min(days) + 0.2 * diff(range(days))
label_y <- y_max - 0.05 * diff(common_ylim) 

# === Updated theme with axis lines (boxed look) ===
boxed_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", size = 0.6),
    panel.border = element_rect(color = "black", fill = NA, size = 0.6)
  )

#Plot A
p1 <- ggplot() +
  geom_point(aes(x = days, y = real_data), shape = 4, size = 2.5, color = "black",stroke = 1.5) +
  geom_line(aes(x = days, y = pred_A), color = "forestgreen", size = 1) +
  labs(title = "High Transmission (USA)", x = "Days", y = "Number of People") +
  annotate("label", x = label_x, y = label_y,
           label = paste0("AIC = ", round(aic_A, 1)),
           size = 6.5, fontface = "bold", fill = "white", label.size = 0.5) +
  ylim(common_ylim) + 
  boxed_theme

#Plot B 
p2 <- ggplot() +
  geom_point(aes(x = days, y = real_data), shape = 4, size = 2.5, color = "black",stroke = 1.5) +
  geom_line(aes(x = days, y = pred_B), color = "firebrick", size = 1) +
  labs(title = "Low Recovery (USA)", x = "Days", y = "Number of People") +
  annotate("label", x = label_x, y = label_y,
           label = paste0("AIC = ", round(aic_B, 1)),
           size = 6.5, fontface = "bold", fill = "white", label.size = 0.5) +
  ylim(common_ylim)+
  boxed_theme

#Combine plots
grid.arrange(p1, p2, ncol = 2)
library(ggplotify)
combined_plot <- as.ggplot(~grid.arrange(p1, p2, ncol = 2))

# Save
ggsave("HipotheticalDisparities-Driven-Covid19.png", plot = combined_plot,
       width = 12, height = 5, dpi = 600, units = "in")
