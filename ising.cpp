// Ising Model Simulation
// Hersh Kumar
// December 2021
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <random>
#include <fstream>

// struct to encode a site on the lattice
struct site {
	int x;
	int y;
	// constructor takes in two positions on the lattice
	site(int a, int b) {
		x = a;
		y = b;
	}
};

// Lattice parameters
const int size = 25; // length of the square lattice
int total_spins = size * size; // total number of sites on the lattice
bool lattice[size - 1][size - 1]; // the actual lattice of sites

int seed = 436675; // seed for the rng
std::random_device rd;
// seed the random number generator
//std::mt19937 gen(seed);
// seeding the generator using hardware:
std::mt19937 gen(rd());


// NOTE: this is necessary because the default rand with a modulo
// does not produce a uniform distribution of values
// We create two distributions, one for the probability that a spin will be flipped
// and the other for picking coordinates on the grid
std::uniform_real_distribution<> prob_dist(0.0, 1.0);
std::uniform_int_distribution<> coord_dist(0, size);

// temperature parameters
float T_init = 10.0; // temperature of the lattice at the beginning
float dt = .1; // change in the temperature per simulation step
float T_min = 0.25; // stopping temperature

// Monte Carlo parameters
int steps = 10000;
int trans_steps = 1000; // transient are the number of steps we ignore at the beginning
//printing symbols
char up = '+';
char down = '-';

// data output files
std::string obs_out_file = "data.csv"; 
std::string lattice_out_file = "latt.txt";

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
	return coord_dist(gen); 
}

// returns a random float between 0 and 1
// Uses the uniform float distribution
float rand_prob() {
	return prob_dist(gen);
}


// returns a random site on the lattice
site rand_site() {
	site ret(rand_coord(), rand_coord());
	return ret;
}
// returns the spin value of a certain site
int get_spin(site s) {
	// true corresponds to up spin, false to down spin
	return lattice[s.x][s.y]; 
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
	e *= -1;
	free(neighbors);
	return e;
}

// flips the spin at a certain site
void flip(site s) {
	lattice[s.x][s.y] = -1 * lattice[s.x][s.y];
}

// function to compute the net total magnetization
int total_mag() {
	int total = 0;
	for (int x = 0; x < size; x++) {
		for (int y = 0; y < size; y++){
			total += get_spin(site(x,y)); 
		}
	}
	return total;
}

// function that computes the total energy:
int total_energy () {
	int total = 0;
	for (int x = 0; x < size; x++) {
		for (int y = 0; y < size; y++){
			total += energy(site(x,y)); 
		}
	}
	return total;
}

// takes the current lattice and writes it to the lattice output file for plotting
void serialize_lattice() {
	// open the file for writing
	std::ofstream lat_out;
	lat_out.open(lattice_out_file, std::fstream::app);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			int val = lattice[i][j];
			lat_out << val << ",";	
		}
		lat_out << "\n";
	}
	lat_out << "\n";
	lat_out.close();
}

// does the transient MC steps that are unstable
// Basically the same as the MC steps except we don't measure a change in the observables
void transient(float T) {
	for (int a = 1; a <= trans_steps; a++) {
		for(int b = 1; b < total_spins; b++) {
			// pick a random site
			site s = rand_site();
			// check whether the site should be flipped based on the energy
			int de = -energy(s);

			bool f = false;

			if (de < 0) {
				f = true;
			}
			else if (rand_prob() < exp(-de / T)) {
				f = true;
			}
			
			// if it should be flipped
			if (f) {
				// flip t he spin
				flip(s);
			}
		}
	}
}

// main loop
int main() {
	// set up the stream for outputting data
	std::ofstream obs_out;
	obs_out.open(obs_out_file);
	// add the data headers to the data file
	obs_out << "Temperature,Average Energy,Average Magnetization" << std::endl; 
	// declare vars for computing the observables
	double E_avg = 0, e_tot = 0;
	double M_avg = 0, m_tot = 0;

	// normalization factor for the observables
	float norm = (1.0/(float) (steps * total_spins));

	// initialize the lattice to a random state
	init_state();
	// display the original lattice in the terminal
	print_state();
	// write the initial state to the file
	serialize_lattice();
	// Now we have the simulation loop	
	for (float T = T_init; T >= T_min; T -= dt) {
		e_tot = 0;
		m_tot = 0;

		// run the transient steps, which are steps we don't care about initially
		transient(T);
		// Begin the main Monte Carlo steps
		for (int step = 1; step <= steps; step++) {

			for (int b = 1; b <= total_spins; b++) {
				// pick a random site
				site s = rand_site();
				// check whether the site should be flipped based on the energy
				int de = -energy(s);
				bool f = false;

				if (de < 0) {
					f = true;
				}
				else if (rand_prob() < exp(-de/T)) {
					f = true;
				}
				// if it should be flipped
				if (f) {
					// flip the spin
					flip(s);
				}				
			}
			e_tot += total_energy();
			m_tot += total_mag();
		}

		//print_state();	
		// scale the observablees by the norm
		E_avg = e_tot/2.0 * norm;
		M_avg = m_tot/2.0 * norm;
		// store the data in a csv file
		obs_out << T << "," << E_avg << "," << M_avg << std::endl;
		// write the lattice to the file
		serialize_lattice();	
	}
	// close the file object
	obs_out.close();
	return 1;
}
