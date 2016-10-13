#include <iostream>
#include <thread>
#include <vector>
#include <cmath>

#include "gen.h"
#include "argparser.h"

using namespace std;

void targetEnergy(Generator g1, int n_configs) {
	int n_good_lattices(0), total_count(0), beta_count(0), local_stuck_count(0);
	double delta_E(0.0);
	vector<double> energies, distances;
	double g1_en = g1.totalCoulombEnergy();
	while (n_good_lattices < n_configs) {
		Generator g2 = g1;
		g2.modifySlightly();
		g2.printDopantAndDefectIndices();
		double g2_en = g2.totalCoulombEnergy();
		double e_prev = abs(g1_en - g1.etarget);
		double e_next = abs(g2_en - g1.etarget);

		if (e_next < e_prev) {
			g1 = g2;
			g1.writeLattice(n_good_lattices);
			g1_en = g2_en;
			n_good_lattices ++;
		} else {
			double uni = g1.rng.random();
			double val = exp( - 0.5 * abs(e_prev - e_next));
			if (val > uni) {
				g1 = g2;
				g1.writeLattice(n_good_lattices);
				g1_en = g2_en;
				n_good_lattices ++;
			} else {
				local_stuck_count ++;
				if (local_stuck_count > 99) {
					g1.placeDefectsAndDopants();
					g1_en = g1.totalCoulombEnergy();
					local_stuck_count = 0;
				}
			}
		}
	}

}

int main(int argc, char ** argv) {
	ArgParser parser(argc, argv);
	int n_atoms = 108, n_configs = parser["number"].as<int>();
	Generator g1(n_atoms, parser);
	targetEnergy(g1, n_configs);
	return 0;
}