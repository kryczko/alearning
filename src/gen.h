#ifndef __GEN__
#define __GEN__

#include <vector>
#include <fstream>
#include <string> 
#include <iostream>
#include <random>

#include "rng.h"
#include "argparser.h"
#include "helper.h"

using namespace std;

class Generator {
	private:
		vector<vector <double>> lattice;
		vector<vector<vector <double>>> md;
		vector<int> defect_indices, dopant_indices, z;
		vector<string> dopants;
		Rng rng;

		int n_atoms, n_dopants, n_defects, seed, n_configs;
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
			int total = this->n_dopants + this->n_defects;
			for (int i = 0; i < total; i ++) {
                int rand_index = -1;
				bool in_dopant_vector = true;
				bool in_defect_vector = true;
				while (in_dopant_vector && in_defect_vector) {
					rand_index = this->rng.randint(0, this->n_atoms);
					in_dopant_vector = find(this->dopant_indices.begin(), this->dopant_indices.end(), rand_index) != this->dopant_indices.end();
					in_defect_vector = find(this->defect_indices.begin(), this->defect_indices.end(), rand_index) != this->defect_indices.end();
				}
                if (i < this->n_dopants) {
                    this->dopant_indices.push_back(rand_index);
                } else {
                    this->defect_indices.push_back(rand_index);
                }
			}
		}

		void printConfig() {
			cout << "\nCONFIG INFO:\n";
			cout << "------------\n";
			cout << "Number of atoms: " << this->n_atoms << endl;
			cout << "Total number of dopant atoms: " << this->n_dopants << endl;
			for (int i = 0; i < this->dopants.size() / 2; i ++) {
				cout << "  Dopant " << i << " atom type: " << this->dopants[2*i] << endl;
				cout << "  Number of dopant atoms: " << this->dopants[2*i + 1] << endl;
			}
			cout << "Number of defects: " << this->n_defects << endl;
			cout << "Energy target: " << this->etarget << endl;
			cout << "Random seed: " << this->seed << endl;
			cout << "Number of configurations: " << this->n_configs << "\n\n";
		}

        double totalCoulombEnergy() {
            vector<int> z(this->n_atoms, 13); // 13 for Al
            int dopant_count = 0;
            for (int i = 0; i < this->dopants.size() / 2; i ++) {
                string name = this->dopants[2 * i];
                int count = stoi(this->dopants[2 * i + 1]);
                for (int j = 0; j < count; j ++) {
                    z[this->dopant_indices[dopant_count]] = coulomb_mapping.at(name);
                }
            }
            for (int index : this->defect_indices) {
                z[index] = 0;
            }

            double toten = 0.0;
            for (int i = 0; i < this->n_atoms; i ++) {
                for (int j = 0; j < this->n_atoms; j ++) {
                    if (i != j) {
                        vector<double> atom1 = this->lattice[i];
                        vector<double> atom2 = this->lattice[j];
                        double dist = r(atom1, atom2);
                        toten += (z[i] * z[j]) / dist;
                    }
                }
            }
            return toten / 2.0;
        }

		Generator(int n_atoms, ArgParser parser) {
			this->n_atoms = n_atoms;
			this->n_defects = parser["defects"].as<int>();
			this->seed = parser["seed"].as<int>();
			this->n_configs = parser["number"].as<int>();
			this->dopants = parser["dopants"].as<vector<string>>();
			this->n_dopants = 0;
			for (int i = 0; i < this->dopants.size() / 2; i ++) {
				this->n_dopants += stoi(this->dopants[2*i + 1]);
			}
			this->etarget = parser["etarget"].as<double>();
			this->rng = Rng(parser["seed"].as<int>());
            this->printConfig();

			this->readXyzFile();
			this->setLattice();
            this->placeDefectsAndDopants();
		} 
};


#endif