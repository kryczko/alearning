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
#include "boost/filesystem/operations.hpp"

namespace fs = boost::filesystem;

using namespace std;

class Generator {
    private:
        vector<vector <double>> lattice;
        vector<vector <double>> forces;
        vector<vector<vector <double>>> md;
        vector<int> defect_indices, dopant_indices, z;
        vector<string> dopants;

        int n_atoms, n_dopants, n_defects, seed;
    public:
        Rng rng;
        double etarget;

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
                    if (local_count % this->n_atoms == 0) {
                        this->md.push_back(time_step);
                        time_step.clear();
                        time_steps ++;
                    }
                }
            }
        }

        void setSeed(int seed) {
            this->seed = seed;
            this->rng = Rng(seed);

        }
        int getSeed() {
            return this->seed;
        }

        void setLattice() {
            int rnd_time_step = this->rng.randint(0, this->md.size());
            this->lattice = md[rnd_time_step];
        }

        void printDopantAndDefectIndices() {
            cout << "Defects: ";
            for (int index : this->defect_indices) {
                cout << index << "  ";
            }
            cout << "\nDopants: ";
            for (int index : this->dopant_indices) {
                cout << index << "  ";
            }
            cout << "\n";
        }

        void printLattice() {
            for (auto& coords : this->lattice) {
                cout << coords[0] << "\t" << coords[1] << "\t" << coords[2] << "\n";
            }
        }

        void placeDefectsAndDopants() {
            this->dopant_indices.clear();
            this->defect_indices.clear();
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
                        double dist = r(atom1, atom2) * 10e-10;
                        double numer = e * e * J_to_eV * k * z[i] * z[j];
                        toten += numer / dist;
                    }
                }
            }
            return toten / 2.0;
        }

        void modifySlightly() {
            int total = this->n_defects + this->n_dopants;
            if (!total) {
                return;
            }
            vector<vector<int>> split_dopants;
            int dopant_count = 0;
            for (int i = 0; i < this->dopants.size() / 2; i ++) {
                vector<int> dopant_indices;
                for (int j = dopant_count; j < stoi(this->dopants[i*2 + 1]) + dopant_count; j ++) {
                    dopant_indices.push_back(this->dopant_indices[j]);
                }
                dopant_count += stoi(this->dopants[i*2 + 1]);
                split_dopants.push_back(dopant_indices);
            }
            vector<int> all_indices;
            for (int val : this->defect_indices) {
                all_indices.push_back(val);
            }
            for (int val: this->dopant_indices) {
                all_indices.push_back(val);
            }
            int rand_int = this->rng.randint(0, total);
            int choice = all_indices[rand_int];
            string typeof = "null";
            vector<int> ignored_indices;
            if (rand_int <= this->defect_indices.size() && this->defect_indices.size() != 0) {
                typeof = "defect";
                ignored_indices = this->defect_indices;
            } else {
                typeof = "dopant";
                int new_rand_int = rand_int - this->defect_indices.size();
                for (int i = 0; i < this->dopants.size() / 2; i ++) {
                    if (new_rand_int - stoi(this->dopants[i * 2 + 1]) < 0) {
                        for (int j = new_rand_int; j < stoi(this->dopants[i * 2 + 1]); j ++ ) {
                            ignored_indices = split_dopants[i];
                        }
                        break;
                    }
                    new_rand_int -= stoi(this->dopants[i * 2 + 1]);
                }
            }
            bool valid_choice = false;
            int rand_atom = -1;
            while (!valid_choice) {
                rand_atom = this->rng.randint(0, this->n_atoms);
                valid_choice = find(ignored_indices.begin(), ignored_indices.end(), rand_atom) == ignored_indices.end();
            }

            if (typeof == "defect") {
                auto it = find(this->defect_indices.begin(), this->defect_indices.end(), choice);
                int index = distance(this->defect_indices.begin(), it);
                this->defect_indices[index] = rand_atom;
            } else if (typeof == "dopant") {
                auto it = find(this->dopant_indices.begin(), this->dopant_indices.end(), choice);
                int index = distance(this->dopant_indices.begin(), it);
                this->dopant_indices[index] = rand_atom;
            }
        }

        double areaOfDopantsAndDefects() {
            vector<int> all_indices;
            for (int val : this->defect_indices) {
                all_indices.push_back(val);
            }
            for (int val: this->dopant_indices) {
                all_indices.push_back(val);
            }
            vector<double> distances;
            for (int i = 0; i < all_indices.size(); i ++) {
                for (int j = 0; j < all_indices.size(); j ++) {
                    if (i != j) {
                        distances.push_back(r(this->lattice[all_indices[i]], this->lattice[all_indices[j]]));
                    }
                }
            }
            double area = distances[0];
            for (int i = 1; i < distances.size(); i ++) {
                area *= distances[i];
            }
            return sqrt(area);
        }

        void calculateForces() {
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
            for (int i = 0; i < this->lattice.size(); i ++) {
                double fx(0.0), fy(0.0), fz(0.0);
                for (int j = 0; j < this->lattice.size(); j++) {
                    if (i != j) {
                        double xdist = pbc_wrap(this->lattice[i][0] - this->lattice[j][0]) * 10e-10;
                        double ydist = pbc_wrap(this->lattice[i][1] - this->lattice[j][1]) * 10e-10;
                        double zdist = pbc_wrap(this->lattice[i][2] - this->lattice[j][2]) * 10e-10;
                        double dist = r(this->lattice[i], this->lattice[j]) * 10e-10;
                        double denom = dist * dist * dist;
                        double numer = e * e * N_to_eV_per_A * k * z[i] * z[j];

                        fx += numer * xdist / denom;
                        fy += numer * ydist / denom;
                        fz += numer * zdist / denom;                        
                    }
                }
                vector<double> force = {fx, fy, fz};
                this->forces.push_back(force);
             }
        }

        void writeXSF(int lattice_ind, double etarget, double toten) {
            this->calculateForces();
            if (!fs::exists("gen")) {
                fs::create_directory("gen");
            }
            if (!fs::exists("gen/xsf")) {
                fs::create_directory("gen/xsf");
            }
            ofstream output;
            output.open(("gen/xsf/gen_" + to_string(etarget) + "_" + to_string(lattice_ind) + ".xsf").c_str());
            output << "# total energy = " << toten << " eV\n\nCRYSTAL\nPRIMVEC\n";
            output << to_string(pbc) + " 0.000 0.000\n0.000 " + to_string(pbc) + " 0.000\n0.000 0.000 " + to_string(pbc) + "\n";
            string counts = to_string(this->n_atoms - this->n_defects - this->n_dopants) + "  ";
            output  << "PRIMCOORD\n";
            output << this->n_atoms - this->n_defects << "  1\n"; 
            for (int i = 0; i < this->lattice.size(); i ++) {
                bool in_dopants = find(this->dopant_indices.begin(), this->dopant_indices.end(), i) != this->dopant_indices.end();
                bool in_defects = find(this->defect_indices.begin(), this->defect_indices.end(), i) != this->defect_indices.end();
                if (!in_dopants && !in_defects) {
                    output << "Al  " << this->lattice[i][0] << "  " << this->lattice[i][1] << "  " << this->lattice[i][2] << "  " << this->forces[i][0] << "  " << this->forces[i][1] << "  " << this->forces[i][2] << "\n";
                }
            }
            int total = 0;
            for (int i = 0; i < this->dopants.size() / 2; i ++) {
                for (int j = total; j < stoi(this->dopants[2*i + 1]); j ++) {
                    output << this->dopants[2*i] << "  " << this->lattice[this->dopant_indices[j]][0] << "  " << this->lattice[this->dopant_indices[j]][1] << "  " << this->lattice[this->dopant_indices[j]][2] << "  " << this->forces[this->dopant_indices[j]][0] << "  " << this->forces[this->dopant_indices[j]][1] << "  " << this->forces[this->dopant_indices[j]][2] << endl;
                }
                total += stoi(this->dopants[2*i + 1]);
            }
            output.close();

        }

        void writePOSCAR(int lattice_ind, double etarget) {
            if (!fs::exists("gen/poscar")) {
                fs::create_directory("gen/poscar");
            }
            ofstream output;
            output.open(("gen//poscar/POSCAR_gen_" + to_string(etarget) + "_" + to_string(lattice_ind)).c_str());
            output << "generated-lattice\n1.0\n";
            output << to_string(pbc) + " 0.000 0.000\n0.000 " + to_string(pbc) + " 0.000\n0.000 0.000 " + to_string(pbc) + "\n";
            output << "Al";
            string names = "  ";
            string counts = to_string(this->n_atoms - this->n_defects - this->n_dopants) + "  ";
            for (int i = 0; i < this->dopants.size() / 2; i ++) {
                names += dopants[2*i] + "  ";
                counts += dopants[2*i + 1] + "  ";
            }
            output << names << "\n";
            output << counts << "\nCartesian\n";
            for (int i = 0; i < this->lattice.size(); i ++) {
                bool in_dopants = find(this->dopant_indices.begin(), this->dopant_indices.end(), i) != this->dopant_indices.end();
                bool in_defects = find(this->defect_indices.begin(), this->defect_indices.end(), i) != this->defect_indices.end();
                if (!in_dopants && !in_defects) {
                    output << this->lattice[i][0] << "  " << this->lattice[i][1] << "  " << this->lattice[i][1] << "\n";
                }
            }

            for (int i : this->dopant_indices) {
                output << this->lattice[i][0] << "  " << this->lattice[i][1] << "  " << this->lattice[i][1] << "\n"; 
            }
            output.close();

        }

        Generator(int n_atoms, ArgParser parser) {
            this->n_atoms = n_atoms;
            this->n_defects = parser["defects"].as<int>();
            this->seed = parser["seed"].as<int>();
            this->dopants = parser["dopants"].as<vector<string>>();
            this->n_dopants = 0;
            for (int i = 0; i < this->dopants.size() / 2; i ++) {
                this->n_dopants += stoi(this->dopants[2*i + 1]);
            }
            this->etarget = parser["etarget"].as<double>();
            this->rng = Rng(parser["seed"].as<int>());
            this->readXyzFile();
            this->setLattice();
            this->placeDefectsAndDopants();
            this->printConfig();

        } 
};


#endif
