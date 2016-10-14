#ifndef __RUNNER__
#define __RUNNER__

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "gen.h"

void targetEnergy(Generator g1, int n_configs) {
    int n_good_lattices(0), total_count(0), beta_count(0), local_stuck_count(0);
    double delta_E(0.0);
    vector<double> energies, distances;
    double g1_en = g1.totalCoulombEnergy();
    while (n_good_lattices < n_configs) {
        Generator g2(g1);
        g2.modifySlightly();
        // g2.printDopantAndDefectIndices();
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

vector<double> runningAverage(vector<double> stuff) {
    vector<double> ra;
    double running_average = 0.0;
    for (int i = 0; i < stuff.size(); i ++) {
        running_average *= i;
        running_average += stuff[i];
        running_average /= i + 1;
        ra.push_back(running_average);
    }
    return ra;
}

void equilTest(Generator g, int n_configs) {
    vector<double> beta_vals = {0, 0.001, 0.01, 0.1, 1, 10};

    vector<vector<double>> all_energies(beta_vals.size()), all_areas(beta_vals.size());

    #pragma omp parallel for
    for (int i = 0; i < beta_vals.size(); i ++) {
        Generator g1(g);
        g1.setSeed(g.getSeed() + i);
        g1.setLattice();
        g1.placeDefectsAndDopants();
        double beta = beta_vals[i];
        cout << "Working on beta val: " << beta << " (thread " << omp_get_thread_num() << ")\n";
        int n_good_lattices(0), total_count(0), beta_count(0), local_stuck_count(0);
        double delta_E(0.0);
        vector<double> energies, areas;
        double g1_en = g1.totalCoulombEnergy();
        while (n_good_lattices < n_configs) {
            Generator g2(g1);
            g2.modifySlightly();
            // g2.printDopantAndDefectIndices();
            double g2_en = g2.totalCoulombEnergy();
            double e_prev = g1_en;
            double e_next = g2_en;

            if (e_next < e_prev) {
                g1 = g2;
                // g1.writeLattice(n_good_lattices);
                g1_en = g2_en;
                energies.push_back(g1_en);
                areas.push_back(g1.areaOfDopantsAndDefects());
                n_good_lattices ++;
            } else {
                double uni = g1.rng.random();
                double val = exp( - beta * abs(e_prev - e_next));
                if (val > uni) {
                    g1 = g2;
                    // g1.writeLattice(n_good_lattices);
                    g1_en = g2_en;
                    energies.push_back(g1_en);
                    areas.push_back(g1.areaOfDopantsAndDefects());
                    n_good_lattices ++;
                } 
                // else {
                //     local_stuck_count ++;
                //     if (local_stuck_count > 999) {
                //         g1.placeDefectsAndDopants();
                //         g1_en = g1.totalCoulombEnergy();
                //         local_stuck_count = 0;
                //     }
                // }
            }
        }
        all_energies[i] = energies;
        all_areas[i] = areas;
    }
    if (!fs::exists("output")) {
        fs::create_directory("output");
    }
    ofstream output;
    vector<vector<double>> energy_ras(beta_vals.size()), area_ras(beta_vals.size());
    for (int i = 0; i < beta_vals.size(); i ++) {
        energy_ras[i] = runningAverage(all_energies[i]);
        area_ras[i] = runningAverage(all_areas[i]);
    }
    output.open("output/running_E_averages.dat");
    for (int i = 0; i < beta_vals.size(); i ++) {
        output << "# " << beta_vals[i] << "  ";
    }
    output << "\n";
    for (int i = 0; i < n_configs; i ++) {
        for (int j = 0; j < beta_vals.size(); j ++) {
            output << energy_ras[j][i] << "  "; 
        }
        output << "\n";
    }

    output.close();
}


void tempSweep(Generator g, int n_configs, bool flipped=false) {
    vector<double> beta_vals;
    int n_betas = 100;
    double start = 0.0, end = 10.0, incr = (end - start) / n_betas;
    for (int i = 0; i < n_betas; i ++) {
        beta_vals.push_back(start + incr*i);
    }

    vector<vector<double>> all_energies(n_betas), all_areas(n_betas);

    if (!flipped) {
        #pragma omp parallel for
        for (int i = 0; i < beta_vals.size(); i ++) {
            Generator g1(g);
            g1.setSeed(g.getSeed() + i);
            g1.setLattice();
            g1.placeDefectsAndDopants();
            double beta = beta_vals[i];
            cout << "Working on beta val: " << beta << " (thread " << omp_get_thread_num() << ")\n";
            int n_good_lattices(0), total_count(0), beta_count(0), local_stuck_count(0);
            double delta_E(0.0);
            vector<double> energies, areas;
            double g1_en = g1.totalCoulombEnergy();
            while (n_good_lattices < n_configs) {
                Generator g2(g1);
                g2.modifySlightly();
                // g2.printDopantAndDefectIndices();
                double g2_en = g2.totalCoulombEnergy();
                double e_prev = g1_en;
                double e_next = g2_en;

                if (e_next < e_prev) {
                    g1 = g2;
                    // g1.writeLattice(n_good_lattices);
                    g1_en = g2_en;
                    energies.push_back(g1_en);
                    areas.push_back(g1.areaOfDopantsAndDefects());
                    n_good_lattices ++;
                } else {
                    double uni = g1.rng.random();
                    double val = exp( - beta * abs(e_prev - e_next));
                    if (val > uni) {
                        g1 = g2;
                        // g1.writeLattice(n_good_lattices);
                        g1_en = g2_en;
                        energies.push_back(g1_en);
                        areas.push_back(g1.areaOfDopantsAndDefects());
                        n_good_lattices ++;
                    } 
                    // else {
                    //     local_stuck_count ++;
                    //     if (local_stuck_count > 999) {
                    //         g1.placeDefectsAndDopants();
                    //         g1_en = g1.totalCoulombEnergy();
                    //         local_stuck_count = 0;
                    //     }
                    // }
                }
            }
            all_energies[i] = energies;
            all_areas[i] = areas;
        }
    } else {
        #pragma omp parallel for
        for (int i = 0; i < beta_vals.size(); i ++) {
            Generator g1(g);
            g1.setSeed(g.getSeed() + i);
            g1.setLattice();
            g1.placeDefectsAndDopants();
            double beta = beta_vals[i];
            cout << "Working on beta val: " << beta << " (thread " << omp_get_thread_num() << ")\n";
            int n_good_lattices(0), total_count(0), beta_count(0), local_stuck_count(0);
            double delta_E(0.0);
            vector<double> energies, areas;
            double g1_en = g1.totalCoulombEnergy();
            while (n_good_lattices < n_configs) {
                Generator g2(g1);
                g2.modifySlightly();
                // g2.printDopantAndDefectIndices();
                double g2_en = g2.totalCoulombEnergy();
                double e_prev = g1_en;
                double e_next = g2_en;

                if (e_next > e_prev) {
                    g1 = g2;
                    // g1.writeLattice(n_good_lattices);
                    g1_en = g2_en;
                    energies.push_back(g1_en);
                    areas.push_back(g1.areaOfDopantsAndDefects());
                    n_good_lattices ++;
                } else {
                    double uni = g1.rng.random();
                    double val = exp( - beta * abs(e_prev - e_next));
                    if (val > uni) {
                        g1 = g2;
                        // g1.writeLattice(n_good_lattices);
                        g1_en = g2_en;
                        energies.push_back(g1_en);
                        areas.push_back(g1.areaOfDopantsAndDefects());
                        n_good_lattices ++;
                    } 
                    // else {
                    //     local_stuck_count ++;
                    //     if (local_stuck_count > 999) {
                    //         g1.placeDefectsAndDopants();
                    //         g1_en = g1.totalCoulombEnergy();
                    //         local_stuck_count = 0;
                    //         exit(-1);
                    //     }
                    // }
                }
            }
            all_energies[i] = energies;
            all_areas[i] = areas;
        }
    }

    if (!fs::exists("output")) {
        fs::create_directory("output");
    }
    ofstream output;
    output.open("output/beta_vs_E_vs_A.dat");
    for (int i = 0; i < beta_vals.size(); i ++) {
        double avg_en = accumulate(all_energies[i].begin(), all_energies[i].end(), 0.0) / all_energies[i].size();
        double avg_area = accumulate(all_areas[i].begin(), all_areas[i].end(), 0.0) / all_areas[i].size();
        output << beta_vals[i] << "\t" << avg_en << "\t" << avg_area << "\n";
    }
    output.close();

}

#endif