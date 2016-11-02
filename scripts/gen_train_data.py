#!/usr/bin/env python

import argparse
import numpy as np
import ase.io as aseio
import random
import os
import copy

# constants
pbc = 8.568
k = 1 / ( 4 * np.pi * 8.85e-12)
j_ev_conv = 1.6e-19
n_atoms = 108
beta =  0.11 # 300 K times the Boltzmann constant in eV

coulomb_mapping = {
    'H': 1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B': 5,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'Ar': 18,
    'K': 19,
    'Ca': 20,
    'Sc': 21,
    'Ti': 22,
    'V': 23,
    'Cr': 24,
    'Mn': 25,
    'Fe': 26,
    'Co': 27,
    'Ni': 28,
    'Cu': 29,
    'Zn': 30,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
    'Br': 35,
    'Kr': 36,
    'Rb': 37,
    'Sr': 38,
    'Y': 29,
    'Zr': 30,
    'Nb': 31,
    'Mo': 32
}


def parseArgs():
    parser = argparse.ArgumentParser(description='Generate FCC Al lattices of user defined sizes.\nExample usage: $ python gen_train_data.py -de 1 -do Ti 5 O 3 -n 100')
    parser.add_argument('--defects', '-de', help='Number of defects to introduce into the lattice.', dest='defects', default=0)
    parser.add_argument('--dopants', '-do', help='Dopants type followed by number of dopant atoms to insert into lattice.', dest='dopants', default=[], nargs='+')
    parser.add_argument('--number', '-n', help='Number of configurations to generate', dest='number', default=1)
    parser.add_argument('--etarget', '-et', help='Target energy in eV when sampling lattices.', dest='etarget', default=-1000)
    parser.add_argument('--seed', '-s', help='Random seed.', dest='seed', default=13)
    return parser.parse_args()


def formatDopantsFromArgs(args):
    dopants = []
    for i in range(len(args.dopants) // 2):
        dopants.append((args.dopants[2*i], int(args.dopants[2*i + 1])))
    return dopants

class Generator:
    def __init__(self, args):
        self.lattice = np.empty([n_atoms, 3])
        self.n = int(args.number)
        self.args = args
        self.defect_indices = []
        self.dopant_indices = []
        self.prev_en = None
        self.traj = aseio.read('static/Al_MD.xyz', index=':', format='xyz')
        self.n_dopants = sum([count for name, count in self.args.dopants])
        self.n_defects = int(self.args.defects)

        print 'INFO: There are', len(self.traj), 'frame(s) in the MD file:', 'static/Al_MD.xyz'

    def getLattice(self):
        # read an xyz file and pick a random traj
        index = random.randint(0, len(self.traj) - 1)
        print 'INFO: Choosing trajectory:', index

        atoms = self.traj[index]
        for i, atom in enumerate(atoms):
            self.lattice[i][0] = atom.position[0]
            self.lattice[i][1] = atom.position[1]
            self.lattice[i][2] = atom.position[2]

    def placeSite(self, count, used_indices=[]):
        indices = []
        for i in range(count):
            indices.append(random.choice(list(set([i for i in range(len(self.lattice))]) - set(used_indices))))
        return indices


    def placeDefectsAndDopants(self):
        self.defect_indices = []
        self.dopant_indices = []
        total =  self.n_defects + self.n_dopants
        defect_indices, dopant_indices = [], []
        if self.n_defects and self.n_dopants:
            if random.random() <= 0.5:
                defect_indices = self.placeSite(self.n_defects)
                dopant_indices = self.placeSite(self.n_dopants, used_indices=defect_indices)
            else:
                dopant_indices = self.placeSite(self.n_dopants)
                defect_indices = self.placeSite(self.n_defects, used_indices=dopant_indices)
        elif not self.n_defects and self.n_dopants:
            dopant_indices = self.placeSite(self.n_dopants)
        elif self.n_defects and not self.n_dopants:
            defect_indices = self.placeSite(self.n_defects)
        self.defect_indices = defect_indices
        self.dopant_indices = dopant_indices

    def pbc_wrap(self, array):
        for i, val in enumerate(array):
            if val > pbc / 2.0:
                array[i] -= pbc
        return array

    def calculateEnergy(self):
        z = np.empty(n_atoms)
        z.fill(13)
        total_dopant_count = 0
        for dopant, count in self.args.dopants:
            for i in range(count):
                z[self.dopant_indices[total_dopant_count]] = coulomb_mapping[dopant]
                total_dopant_count += 1
        for index in self.defect_indices:
            z[index] = 0.

        toten = 0.

        for i in range(n_atoms):
            for j in range(n_atoms):
                if i != j:
                    r = self.pbc_wrap(abs(self.lattice[i] - self.lattice[j]))
                    toten += ( z[i] * z[j] ) / np.linalg.norm(r, 1)
        return toten / 2.

    def writeLattice(self, lattice_ind):
        if not os.path.exists('gen'):
            os.mkdir('gen')
        f = open('gen/POSCAR_gen_' + str(lattice_ind).zfill(3), 'w')
        f.write('generated-lattice\n1.0\n')
        f.write(str(pbc) + ' 0.000 0.000\n0.000 ' + str(pbc) + ' 0.000\n0.000 0.000 ' + str(pbc) + '\n')
        f.write('Al  ' + '  '.join([name for name, count in self.args.dopants]) + '\n')
        f.write(str(n_atoms - int(self.args.defects) - sum([count for name, count in self.args.dopants])) + '  ' + '  '.join([str(count) for name, count in self.args.dopants]) + '\nCartesian\n')


        for i, row in enumerate(self.lattice):
            if i not in self.defect_indices + self.dopant_indices:
                f.write(' '.join([str(val) for val in self.lattice[i]]) + '\n')
        dopant_count = 0
        for name, count in self.args.dopants:
            for i in range(count):
                f.write(' '.join([str(val)
                 for val in self.lattice[self.dopant_indices[dopant_count]]]) + '\n')
                dopant_count += 1
        f.close()

    def avgDist(self):
        indices = self.defect_indices + self.dopant_indices
        total = 0.
        for i in indices:
            for j in indices:
                total += np.linalg.norm(self.pbc_wrap(abs(self.lattice[i] - self.lattice[j])), 2)
        return 0.5 * total / len(indices)

    def modifySlighty(self):
        total =  self.n_defects + self.n_dopants
        # select random dopant/defect
        randint = random.randint(0, total - 1)
        indices =  self.defect_indices + self.dopant_indices
        choice = indices[randint]
        # select another atom/defect of different type
        set_to_swap_with = set([i for i in range(n_atoms)])
        typeof = None
        if randint <= len(self.defect_indices) - 1:
            set_to_swap_with -= set(self.defect_indices)
            typeof = 'defect'
        else:
            typeof = 'dopant'
            new_randint = randint - (len(self.dopant_indices) - 1)
            for atom, count in self.args.dopants:
                if new_randint - (count - 1) <= 0:
                    set_to_swap_with -= set(self.dopant_indices[randint:randint + count  - 1])
                new_randint -= (count - 1)
        next_choice = random.choice(list(set_to_swap_with))
        while next_choice == choice:
            next_choice = random.choice(list(set_to_swap_with))
        # swap them
        if typeof == 'defect':
            self.defect_indices[self.defect_indices.index(choice)] = next_choice
        elif typeof == 'dopant':
            self.dopant_indices[self.dopant_indices.index(choice)] = next_choice

    def startTargetEnergy(self):
        etarget = -float(self.args.etarget)
        n_good_lattices = 0
        total_count = 0
        beta_count = 0
        delta_E = 0.
        energies = []
        avg_distances = []
        self.getLattice()
        self.placeDefectsAndDopants()
        self.prev_en = self.calculateEnergy()
        avg_distances.append(self.avgDist())
        energies.append(self.prev_en)
        while n_good_lattices < self.n:
            trial = copy.deepcopy(self)
            trial.modifySlighty()
            en = trial.calculateEnergy()
            e_next = abs(etarget - en)
            e_prev = abs(etarget - self.prev_en)
            if e_next < e_prev:
                self = copy.deepcopy(trial)
                self.writeLattice(n_good_lattices)
                avg_distances.append(self.avgDist())
                energies.append(en)
                print 'INFO: Total Coulomb energy of lattice: %.10E (iter %i)' % (en, n_good_lattices)
                print 'INFO: Found better energy, delta E = %.10E' % abs(e_next - e_prev)
                delta_E += abs(e_next - e_prev)
                self.prev_en = en
                n_good_lattices += 1
            elif e_next >= e_prev:
                total_count += 1
                uni = random.random()
                val = np.exp(- beta * abs(e_next - e_prev))
                if val > uni:
                    beta_count += 1
                    self = copy.deepcopy(trial)
                    self.writeLattice(n_good_lattices)
                    avg_distances.append(self.avgDist())
                    energies.append(en)
                    print 'INFO: Total Coulomb energy of lattice: %.10E (iter %i)' % (en, n_good_lattices)
                    print 'INFO: Hopped here, delta E = %.10E' % abs(e_next - e_prev)
                    delta_E += abs(e_next - e_prev)
                    n_good_lattices += 1
                    self.prev_en = en

        print
        print
        print "INFO: Jumped %0.1f%% of the time" % (100. * float(beta_count) / total_count )
        print "INFO: Average delta E: %.10E" % (delta_E / self.n  )        
        self.writeDistsVsEnergies(energies, avg_distances)

    def writeDistsVsEnergies(self, energies, avg_distances):
        if not os.path.exists('output'):
            os.mkdir('output')
        f = open('output/E_vs_avg_r.dat', 'w')
        for energy, dist in zip(energies, avg_distances):
            f.write('%.10E\t%.10E\n' % (dist, energy))
        f.close()

    def writeBetaStats(self, energies, avg_distances, betas):
        if not os.path.exists('output'):
            os.mkdir('output')
        f = open('output/beta_e_d.dat', 'w')
        for energy, dist, beta in zip(energies, avg_distances, betas):
            f.write('%.10E\t%.10E\t%.10E\n' % (beta, np.mean(dist), np.mean(energy)))
        f.close()

    def startTempSweep(self, flipped=False):
        all_energies = []
        all_dists = []
        for incr in np.arange(0,1,0.01):
            print 'INFO: beta =', incr
            n_good_lattices = 0
            total_count = 0
            beta_count = 0
            delta_E = 0.
            energies = []
            avg_distances = []
            self.getLattice()
            self.placeDefectsAndDopants()
            self.prev_en = self.calculateEnergy()
            avg_distances.append(self.avgDist())
            energies.append(self.prev_en)
            local_stuck_count = 0
            if not flipped:
                while n_good_lattices < self.n:
                    trial = copy.deepcopy(self)
                    trial.modifySlighty()
                    en = trial.calculateEnergy()
                    e_next = en
                    e_prev = self.prev_en
                    if e_next < e_prev:
                        local_stuck_count = 0
                        self = copy.deepcopy(trial)
                        self.writeLattice(n_good_lattices)
                        avg_distances.append(self.avgDist())
                        energies.append(en)
                        self.prev_en = en
                        n_good_lattices += 1
                    elif e_next >= e_prev:
                        total_count += 1
                        uni = random.random()
                        val = np.exp(- incr * abs(e_next - e_prev))
                        if val > uni:
                            local_stuck_count = 0
                            beta_count += 1
                            self = copy.deepcopy(trial)
                            self.writeLattice(n_good_lattices)
                            avg_distances.append(self.avgDist())
                            energies.append(en)
                            delta_E += abs(e_next - e_prev)
                            n_good_lattices += 1
                            self.prev_en = en
                        else:
                            local_stuck_count += 1
                            if local_stuck_count > 100:
                                self.placeDefectsAndDopants()
                                self.prev_en = self.calculateEnergy()
                                avg_distances.append(self.avgDist())
                                energies.append(self.prev_en)
                                local_stuck_count = 0
                                n_good_lattices += 1
            else:
                while n_good_lattices < self.n:
                    trial = copy.deepcopy(self)
                    trial.modifySlighty()
                    en = trial.calculateEnergy()
                    e_next = en
                    e_prev = self.prev_en
                    if e_next > e_prev:
                        self = copy.deepcopy(trial)
                        self.writeLattice(n_good_lattices)
                        avg_distances.append(self.avgDist())
                        energies.append(en)
                        self.prev_en = en
                        n_good_lattices += 1
                    elif e_next <= e_prev:
                        total_count += 1
                        uni = random.random()
                        val = np.exp(- incr * abs(e_next - e_prev))
                        if val > uni:
                            beta_count += 1
                            self = copy.deepcopy(trial)
                            self.writeLattice(n_good_lattices)
                            avg_distances.append(self.avgDist())
                            energies.append(en)
                            delta_E += abs(e_next - e_prev)
                            n_good_lattices += 1
                            self.prev_en = en
            print
            print
            print "INFO: Jumped %0.1f%% of the time" % (100. * float(beta_count) / total_count )
            print "INFO: Average delta E: %.10E" % (delta_E / self.n  )
            all_energies.append(energies)
            all_dists.append(avg_distances)
        self.writeBetaStats(all_energies, all_dists, np.arange(0,1,0.01))


def main():
    args = parseArgs()
    args.dopants = formatDopantsFromArgs(args)
    random.seed(args.seed)
    generator = Generator(args)
    generator.startTempSweep()

if __name__ == '__main__':
    main()