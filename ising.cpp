// Ising Model Simulation
// Hersh Kumar
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

// Lattice parameters
const int size = 5; // length of the square lattice
int total_spins = size * size; // total number of sites on the lattice
float temperature = 0; // temperature of the lattice at the beginning
int seed = 1; // seed for the rng
bool lattice[size - 1][size - 1]; // the actual lattice of sites

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

// function to get nearest neighbors


// main loop
int main() {
	srand(seed);
	init_state();
	print_state();
	return 1;
}
