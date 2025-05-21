# Install and load necessary packages
if (!require(deSolve)) install.packages("deSolve")
library(deSolve)

# Define the system of differential equations representing the memristive circuit dynamics
# Assuming corrected equations based on interpretation and typical forms
memristor_circuit_ode <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx <- y / C # Assuming dx/dt as current through C
    dy <- -(1/L) * (x + beta * (z^2 - 1) * y) # dy/dt as specified
    dz <- -y - alpha * z + y * z # dz/dt as specified
    list(c(dx, dy, dz))
  })
}

# Define parameters from the paper
params <- c(C = 1.2, L = 3.3, beta = 1.34)

# Set up time settings for simulation
t_start <- 0
t_end <- 500 # Long enough to capture dynamics
dt <- 0.01
times <- seq(t_start, t_end, by = dt)

# Initial conditions
initial_conditions <- c(x = 0.1, y = 0.1, z = 0.1)

# Function to find local maxima of x after discarding initial transients
get_local_maxima <- function(solution_df, transient_time_percentage = 0.5) {
  transient_index <- floor(nrow(solution_df) * transient_time_percentage)
  df_stable <- solution_df[transient_index:nrow(solution_df), ]
  
  # Iterate through the stable part to find local maxima
  x_vals <- df_stable$x
  maxima <- c()
  for (i in 2:(length(x_vals) - 1)) {
    if (x_vals[i] > x_vals[i-1] && x_vals[i] > x_vals[i+1]) {
      maxima <- c(maxima, x_vals[i])
    }
  }
  return(unique(maxima))
}

# Generate bifurcation diagram data
alpha_values <- seq(0.0, 1.2, by = 0.001)
x_maxima_list <- list()

# Progress bar for simulation feedback
pb <- txtProgressBar(min = 0, max = length(alpha_values), style = 3)

for (i in seq_along(alpha_values)) {
  alpha_val <- alpha_values[i]
  current_params <- c(params, alpha = alpha_val) # Add alpha to parameters
  
  # Simulate the system using 'lsoda' for solving ODEs
  sol <- ode(y = initial_conditions, times = times, func = memristor_circuit_ode, parms = current_params, method = "lsoda")
  solution_df <- as.data.frame(sol)
  
  # Discard initial transient dynamics and calculate local maxima of x
  transient_end_time <- 300
  transient_index <- which(solution_df$time >= transient_end_time)[1]
  df_stable <- if (is.na(transient_index)) solution_df else solution_df[transient_index:nrow(solution_df), ]
  
  # Detect upward crossings of the plane z=0
  poincare_x_vals <- c()
  z_vals <- df_stable$z
  x_vals_poincare <- df_stable$x
  
  for (j in 2:length(z_vals)) {
    if (z_vals[j-1] < 0 && z_vals[j] >= 0) { # Upward crossing
      fraction <- -z_vals[j-1] / (z_vals[j] - z_vals[j-1])
      interpolated_x <- x_vals_poincare[j-1] + (x_vals_poincare[j] - x_vals_poincare[j-1]) * fraction
      poincare_x_vals <- c(poincare_x_vals, interpolated_x)
    }
  }
  
  x_maxima_list[[i]] <- unique(poincare_x_vals)
  setTxtProgressBar(pb, i)
}
close(pb)

# Prepare data for plotting the bifurcation diagram
plot_data <- data.frame(alpha = numeric(), x_value = numeric())
for (i in seq_along(alpha_values)) {
  if (length(x_maxima_list[[i]]) > 0) {
    current_alpha <- rep(alpha_values[i], length(x_maxima_list[[i]]))
    current_x_values <- x_maxima_list[[i]]
    plot_data <- rbind(plot_data, data.frame(alpha = current_alpha, x_value = current_x_values))
  }
}

# Plot the bifurcation diagram
plot(plot_data$alpha, plot_data$x_value,
     pch = ".", 
     cex = 1,
     col = "black",
     xlab = expression(alpha),
     ylab = "Local maximum of x",
     main = "Bifurcation Diagram of Memristive Circuit (Replication of Figure 6a)"
)

# Fine analysis around critical alpha value
alpha_fine_range <- seq(1.041, 1.043, by = 0.00001)
fine_x_maxima_list <- list()

pb_fine <- txtProgressBar(min = 0, max = length(alpha_fine_range), style = 3)

for (i in seq_along(alpha_fine_range)) {
  alpha_val <- alpha_fine_range[i]
  current_params <- c(params, alpha = alpha_val)
  
  sol <- ode(y = initial_conditions, times = times, func = memristor_circuit_ode, parms = current_params, method = "lsoda")
  solution_df <- as.data.frame(sol)
  
  transient_end_time <- 300
  transient_index <- which(solution_df$time >= transient_end_time)[1]
  df_stable <- if (is.na(transient_index)) solution_df else solution_df[transient_index:nrow(solution_df), ]
  
  poincare_x_vals <- c()
  z_vals <- df_stable$z
  x_vals_poincare <- df_stable$x
  
  for (j in 2:length(z_vals)) {
    if (z_vals[j-1] < 0 && z_vals[j] >= 0) {
      fraction <- -z_vals[j-1] / (z_vals[j] - z_vals[j-1])
      interpolated_x <- x_vals_poincare[j-1] + (x_vals_poincare[j] - x_vals_poincare[j-1]) * fraction
      poincare_x_vals <- c(poincare_x_vals, interpolated_x)
    }
  }
  fine_x_maxima_list[[i]] <- unique(poincare_x_vals)
  setTxtProgressBar(pb_fine, i)
}
close(pb_fine)

plot_data_fine <- data.frame(alpha = numeric(), x_value = numeric())
for (i in seq_along(alpha_fine_range)) {
  if (length(fine_x_maxima_list[[i]]) > 0) {
    current_alpha <- rep(alpha_fine_range[i], length(fine_x_maxima_list[[i]]))
    current_x_values <- fine_x_maxima_list[[i]]
    plot_data_fine <- rbind(plot_data_fine, data.frame(alpha = current_alpha, x_value = current_x_values))
  }
}

# Plot fine-scale bifurcation diagram around alpha = 1.04
plot(plot_data_fine$alpha, plot_data_fine$x_value,
     pch = ".",
     cex = 1,
     col = "blue",
     xlab = expression(alpha),
     ylab = "Local maximum of x",
     main = "Bifurcation Diagram (Fine Scale around Alpha = 1.04)"
)
abline(v = 1.041574, col = "red", lty = 2) # Mark specified divergence point
abline(v = 1.041574 + 7e-4, col = "green", lty = 2) # Mark second specified point