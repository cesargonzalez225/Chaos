library(deSolve)

params <- list(
  eta = 1,
  D = 1,          # nm
  R_OFF = 70,     # Ohms
  R_ON = 1,       # Ohms
  p = 10,
  v0 = 1,         # V
  T = 20          # s
)

# Initial condition for the state variable w (width of the doped region normalized)
state <- c(w = 0.6)

# Differential equation model describing the memristor dynamics
memristor_model <- function(t, state, params) {
  with(as.list(c(state, params)), {
    w <- as.numeric(w)
    
    # Ensure w remains physically valid between 0 and 1
    w <- min(max(w, 0), 1) 
    
    # Voltage waveform applied to the memristor: sinusoidal with amplitude v0 and period T
    V <- v0 * sin(2 * pi * t / T)
    
    # Window function f_wp that modulates the speed of state variable change
    # It depends on w and parameter p, shaping how w evolves near its bounds
    f_wp <- 1 - (2*w - 1)^(2*p)
    
    # Memristor resistance M as weighted average of ON and OFF resistances based on w
    M <- R_ON * w / D + R_OFF * (1 - w / D)
    
    # Current I through the memristor via Ohm's law
    I <- V / M
    
    # Differential equation for w
    # Rate of change of w depends on eta, window function, and current I
    dw_dt <- eta * f_wp * I
    
    # Return derivative and monitored variables for integration output
    list(c(dw_dt), c(V = V, w = w, M = M, I = I))
  })
}

# Time points at which to numerically solve the ODE (from 0 to 100 in steps of 0.01)
time <- seq(0, 100, by = 0.01)

# Numerical solution of the ODE using the 'ode' function from deSolve package
# This function applies adaptive Runge-Kutta based methods or others to solve ODEs
out <- ode(
  y = state,        # initial state vector
  times = time,     # time points for evaluation
  func = memristor_model,  # model describing derivatives
  parms = params    # parameters required by model
)

# Convert output to data frame for easier plotting and analysis
out <- as.data.frame(out)

# Plot the results: voltage (V), state variable (w), resistance (M), and current (I)
par(mfrow = c(4, 1), mar = c(2, 4, 1, 1))
plot(out$time, out$V, type = "l", col = "blue", ylab = expression(V), xlab = "", main = "(a)")
plot(out$time, out$w, type = "l", col = "blue", ylab = expression(w), xlab = "", main = "(b)")
plot(out$time, out$M, type = "l", col = "blue", ylab = expression(M), xlab = "", main = "(c)")
plot(out$time, out$I, type = "l", col = "blue", ylab = expression(I), xlab = "time (s)", main = "(d)")