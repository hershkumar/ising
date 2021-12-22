// Ising Model Simulation
// Hersh Kumar
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

// struct to encode a site on the lattice
struct site {
	int x;
	int y;

	site(int a, int b) {
		x = a;
		y = b;
	}
};

// Lattice parameters
const int size = 4; // length of the square lattice
int total_spins = size * size; // total number of sites on the lattice
bool lattice[size - 1][size - 1]; // the actual lattice of sites
int J = 1; // J from the Hamiltonian
int seed = 1; // seed for the rng

// temperature parameters
float T_init = 10; // temperature of the lattice at the beginning
float dt = .1;
float T_min = 0.5;

// Monte Carlo parameters
int steps = 10000;
int trans_steps = 1000; // transient are the number of steps we ignore at the beginning
//printing symbols
char up = '+';
char down = '-';


// function to initialize the lattice with spins
bool init_state () {	
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			// set the spin to either 0 or 1
			lattice[i][j] = rand() % 2;
		}
	}
	return true;
}
// function to print the lattice
void print_state() {
	for (int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++) {
			if (lattice[i][j] == 0) {
				std::cout << up << " ";
			}
			else {
				std::cout << down << " ";
			}
		}
		std::cout << std::endl;
	}
	std::cout << "\n";
}

// function to get a random index in the lattice
// since its a square lattice, we can use this for both the x and y indices of the site
int rand_coord() {
	return rand() % size; // modulo restricts the random number to 0 - size-1 
}
// returns a random site on the lattice
site rand_site() {
	site ret(rand_coord(), rand_coord());
	return ret;
}
// returns the spin value of a certain site
int get_spin(site s) {
	// true corresponds to up spin, false to down spin
	return (lattice[s.x][s.y] ? 1 : -1); 
}

// given a site, returns the 4 nearest neighbor sites
// uses periodic boundary conditions to loop around the lattice if its at the edge
struct site* get_neighbors(site s) {
	// allocate memory for the return pointer
	struct site* neighbors = (site*) calloc(4, sizeof(struct site));
	// now look at the neighbors of s
	// first two slots in the array are for x neighbors
	if (s.x != 0) {
		if (s.x != size - 1) {
			neighbors[0] = site(s.x + 1, s.y);
			neighbors[1] = site(s.x - 1, s.y);
		} else {
			neighbors[0] = site(0, s.y);
			neighbors[1] = site(s.x - 1, s.y);
		}	
	} else {
		neighbors[0] = site(s.x + 1, s.y);
		neighbors[1] = site(size - 1, s.y);
	}

	// and the second two slots are for the y neighbors
	if (s.y != 0) {
		if (s.y != size - 1) {
			neighbors[2] = site(s.x, s.y + 1);
			neighbors[3] = site(s.x, s.y - 1);
		} else {
			neighbors[2] = site(s.x, 0);
			neighbors[3] = site(s.x, s.y - 1);
		}
	} else {
		neighbors[2] = site(s.x, s.y + 1);
		neighbors[3] = site(s.x, size - 1);
	}
	return neighbors;
}

// returns the energy at a certain site on the lattice
// uses the Hamiltonian to compute the energy
// H = -J sum_nn S_1 S_j
int energy(site s) {
	int e = 0;
	struct site* neighbors = get_neighbors(s);
	int spin_s = get_spin(s);
	
	for (int i = 0; i < 4; i++) {
		e += spin_s * get_spin(neighbors[i]);	
	}
	e *= -J;
	free(neighbors);
	return e;
}

// flips the spin at a certain site
void flip(site s) {
	lattice[s.x][s.y] = !lattice[s.x][s.y];
}

// function to compute the net total magnetization
int mag_tot() {
	int total = 0;
	for (int x = 0; x < size; x++) {
		for (int y = 0; y < size; y++){
			total += get_spin(site(x,y)); 
		}
	}
	return total;
}

// function that computes the total energy:
int e_tot () {
	int total = 0;
	for (int x = 0; x < size; x++) {
		for (int y = 0; y < size; y++){
			total += energy(site(x,y)); 
		}
	}
	return total;
}

// main loop
int main() {
	// declare vars for computing the observables
	double E = 0, Esq = 0, E_avg = 0, Esq_avg = 0, e_tot = 0, e_tot_sq = 0;
	double M = 0, Msq = 0, M_avg = 0, Msq_avg = 0, m_tot = 0, m_tot_sq = 0;

	// seed the random number generator
	srand(seed);
	// initialize the lattice to a random state
	init_state();
	// display the original lattice in the terminal
	print_state();
	// Now we have the simulation loop	
	for (float T = T_init; T >= T_min; T -= dt) {
		e_tot = 0;
		e_tot_sq = 0;
		m_tot = 0;
		m_tot_sq = 0;

		for (int step = 0; step <= steps; step++) {
			for (int b = 1; b <= total_spins; b++) {
				// pick a random site
				site s = rand_site();
				// check whether the site should be flipped based on the energy
				int de = -2 * J * energy(s);

				bool f;
				if (de < 0) {
					f = true;
				}
				else if (rand() % 2 < exp(-de / T)) {
					f = true;
				}
				else f = false;
				// if it should be flipped
				if (f) {
					// flip the spin
					flip(s);
					// compute a change in the observables
					E += 2 * J * de;
					M += 2 * lattice[s.x][s.y];

				}
			}	
		}
		print_state();	
	}
	return 1;
}
