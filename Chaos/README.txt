README Section


Task1

This script employs the `deSolve` package in R to model and simulate the dynamics of a memristor. The memristor's state evolution is described by differential equations that incorporate physical parameters, such as the width of the doped region and resistance states. Key parameters include `R_ON`, `R_OFF`, and a sinusoidal voltage waveform characterized by amplitude `v0` and period `T`. The simulation evaluates the state variable `w` over time, which influences the memristor's resistance and current flow. Results are plotted to show the voltage, state variable, resistance, and current over a simulated duration, providing insights into the memristor's dynamic responses under varying conditions.


Task2


This script models a circuit using ordinary differential equations (ODEs) to explore its dynamic behavior under different parameter settings. The `deSolve` package in R is used to define the ODE system, where state variables \(x\), \(y\), and \(z\) describe the circuit's voltage, current, and internal states. The script compares two cases: a periodic case with parameters inducing stable limit cycles (\(\alpha = 1.15\)) and a chaotic case resulting in a strange attractor (\(\alpha = 0.85\)). By simulating from time 0 to 200 with precise numerical integration settings, the script outputs phase portraits of state variables \(x\) vs \(y\), illustrating the long-term differences between periodic and chaotic dynamics after discarding initial transients. This visualization helps assess how parameter variations affect circuit behavior.

Task3

This script simulates the dynamics of a memristive circuit using a system of differential equations, with the goal of replicating a bifurcation diagram as a function of the parameter alpha. It employs the `deSolve` package in R to solve these equations and analyze system behavior. The script first defines the equations governing the circuit, then simulates its behavior over a set timeframe while varying \(\alpha\). It detects plane crossings to record significant state changes and generates a bifurcation diagram. Additionally, it conducts a fine-scale analysis around a critical \(\alpha\) value to observe subtle transitions and divergences. This detailed exploration helps in understanding the circuit's transition from periodic to chaotic behavior.