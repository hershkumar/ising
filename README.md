## 2D Ising Model Simulation
This program simulates a 2D lattice of sites with either spin up or down.
We have that each site has 4 nearest neighbors, and the lattice has periodic boundary conditions, so the nearest neighbors wrap around the edge of the lattice.

This implements the Metropolis algorithm, which starts with some random configuration, and then picks a random spot on the lattice. We then compute whether the change in energy of switching that specific spin is lower than 0. If it is lower than 0, we flip the spin, otherwise we probabilistic choose whether or not to flip the spin (via the Boltzmann factor).
We repeat this process N times, each time storing the changes we have in the observables (Average energy, average of energy squared, average magnetization, and average of magnetization squared).
