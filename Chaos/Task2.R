library(deSolve)

# Defines the system of ordinary differential equations (ODEs) representing the circuit dynamics
circuit_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    # State variables:
    # x, y, z represent circuit variables (e.g., voltage, current, internal states)
    
    # Differential equations for x, y, and z
    dx <- y / C
    # Equation models the rate of change of y affected by x, y, z, and parameters L, beta
    dy <- -1/L * (x + beta * (z^2 - 1) * y)
    # Equation models the rate of change of z, a nonlinear function of y and z with parameter alpha
    dz <- -y - alpha * z + y * z
    
    # Return list of derivatives for use in ODE solver
    list(c(dx, dy, dz))
  })
}

# Parameters for the system: capacitance C, inductance L, nonlinear constants beta and alpha
params_a <- list(C = 1.2, L = 3.3, beta = 1.34, alpha = 1.15)  # Parameter set a (periodic case)
params_b <- list(C = 1.2, L = 3.3, beta = 1.34, alpha = 0.85)  # Parameter set b (chaotic case)

# Initial conditions for state variables x, y, and z
state <- c(x = 0.1, y = 0, z = 0.1)

# Time points to simulate over, from time 0 to 200 in steps of 0.01
time <- seq(0, 200, by = 0.01)

# Numerically integrate the ODE system for parameter set a
# 'ode' uses adaptive numerical methods (like Runge-Kutta) to approximate solutions
# rtol and atol control the solver's precision (relative and absolute tolerance)
out_a <- as.data.frame(ode(y = state, times = time, func = circuit_ode, parms = params_a, rtol=1e-8, atol=1e-8))

# Numerically integrate the ODE system for parameter set b (different alpha, chaotic dynamics)
out_b <- as.data.frame(ode(y = state, times = time, func = circuit_ode, parms = params_b, rtol=1e-8, atol=1e-8))

# Discard transient solution to focus on long-term behavior
transient_time <- 50
idx_a <- which(out_a$time > transient_time)
idx_b <- which(out_b$time > transient_time)

# Plot results:
# Phase portraits showing the behavior of state variables x vs y for both parameter sets
par(mfrow = c(2,1), mar = c(4,4,2,1))

# Plot for parameter set a: expected periodic behavior (stable limit cycle)
plot(out_a$x[idx_a], out_a$y[idx_a], type='l', col = "blue", lwd=1,
     xlab="x", ylab="y", xlim=c(-4,2), ylim=c(-2,2), main=expression(bold("(a) periodic ("*alpha*" = 1.15)")))

# Plot for parameter set b: expected chaotic behavior (strange attractor)
plot(out_b$x[idx_b], out_b$y[idx_b], type='l', col = "darkgreen", lwd=1,
     xlab="x", ylab="y", xlim=c(-4,2), ylim=c(-2,2), main=expression(bold("(b) chaotic ("*alpha*" = 0.85)")))