#ifndef __GEN__
#define __GEN__

#include <vector>
#include <fstream>
#include <string> 
#include <iostream>
#include <random>

#include "rng.h"
#include "argparser.h"

using namespace std;

class Generator {
	private:
		vector<vector <double>> lattice;
		vector<vector<vector <double>>> md;
		vector<int> defect_indices, dopants_indices;

		Rng rng;

		int n_atoms, n_dopants, n_defects;
		double etarget;
	public:
		void readXyzFile() {
			ifstream xyzfile;
			xyzfile.open("static/Al_MD.xyz");
			string stuff;
			int time_steps = 0;
			int local_count = 0;
			vector<double> coords(3);
			vector<vector<double>> time_step;
			while ( !xyzfile.eof() ) {
				xyzfile >> stuff;
				if ( stuff == "Al") {
					xyzfile >> coords[0] >> coords[1] >> coords[2];
					time_step.push_back(coords);
					local_count ++;
					if (local_count % n_atoms == 0) {
						md.push_back(time_step);
						time_steps ++;
					}
				}
			}
		}

		void setLattice() {
			int rnd_time_step = this->rng.randint(0, this->md.size());
			this->lattice = md[rnd_time_step];
		}

		void printLattice() {
			for (auto& coords : this->lattice) {
				cout << coords[0] << "\t" << coords[1] << "\t" << coords[2] << "\n";
			}
		}

		void placeDefectsAndDopants() {

		}

		Generator(int n_atoms, ArgParser parser) {
			this->n_atoms = n_atoms;
			this->etarget = parser["etarget"].as<double>();
			this->rng = Rng(parser["seed"].as<int>());
			this->readXyzFile();
			this->setLattice();
			this->printLattice();
		}
};


#endif