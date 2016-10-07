#!/usr/bin/env python

import argparse
import numpy as np
import ase.io as aseio
import random

random.seed(13)

pbc = 8.568

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
    return parser.parse_args()


def formatDopantsFromArgs(args):
    dopants = []
    for i in range(len(args.dopants) // 2):
        dopants.append((args.dopants[2*i], int(args.dopants[2*i + 1])))
    return dopants

class Generator:
    def __init__(self, args):
        self.lattice = np.empty([108, 3])
        self.n = args.number
        self.args = args
        self.defect_indices = []
        self.dopant_indices = []

    def getLattice(self):
        # read an xyz file and pick a random traj
        traj = aseio.read('static/Al_MD.xyz', index=':', format='xyz')
        print 'INFO: There are ', len(traj), ' in the MD file:', 'static/Al_MD.xyz'
        index = random.randint(0, len(traj))
        print 'INFO: Choosing trajectory:', index

        atoms = traj[index]
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
        n_dopants = sum([count for name, count in self.args.dopants])
        n_defects = int(self.args.defects)
        total =  n_defects + n_dopants
        defect_indices, dopant_indices = [], []
        if n_defects and n_dopants:
            if random.random() <= 0.5:
                defect_indices = self.placeSite(n_defects)
                dopant_indices = self.placeSite(n_dopants, used_indices=defect_indices)
            else:
                dopant_indices = self.placeSite(n_dopants)
                defect_indices = self.placeSite(n_defects, used_indices=dopant_indices)
        elif not n_defects and n_dopants:
            dopant_indices = self.placeSite(n_dopants)
        elif n_defects and not n_dopants:
            defect_indices = self.placeSite(n_defects)
        self.defect_indices = defect_indices
        self.dopant_indices = dopant_indices

    def pbc_wrap(self, array):
        for i, val in enumerate(array):
            if val > pbc / 2.0:
                array[i] -= pbc
        return array

    def calculateEnergy(self):
        z = np.empty(108)
        z.fill(13)
        total_dopant_count = 0
        for dopant, count in self.args.dopants:
            for i in range(count):
                z[self.dopant_indices[total_dopant_count]] = coulomb_mapping[dopant]
                total_dopant_count += 1
        for index in self.defect_indices:
            z[index] = 0.

        toten = 0.

        for i in range(108):
            for j in range(108):
                if i != j:
                    r = self.pbc_wrap(abs(self.lattice[i] - self.lattice[j]))
                    toten += - ( z[i] * z[j] ) / np.linalg.norm(r, 2)
        return toten / 2.

    def writeLattice(self):
        f = open('POSCAR_generated', 'w')
        f.write('generated-lattice\n1.0\n')
        f.write(str(pbc) + ' 0.000 0.000\n0.000 ' + str(pbc) + ' 0.000\n0.000 0.000 ' + str(pbc) + '\n')
        f.write('Al  ' + '  '.join([name for name, count in self.args.dopants]) + '\n')
        f.write(str(108 - int(self.args.defects) - sum([count for name, count in self.args.dopants])) + '  ' + '  '.join([str(count) for name, count in self.args.dopants]) + '\nCartesian\n')


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
                

    def start(self):
        for i in range(self.n):
            print "Generating configuration:", i
            self.getLattice()
            self.placeDefectsAndDopants()
            en = self.calculateEnergy()
            print 'Total energy of lattice:', en, 'J'
            self.writeLattice()



def main():
    args = parseArgs()
    args.dopants = formatDopantsFromArgs(args)
    generator = Generator(args)
    generator.start()

if __name__ == '__main__':
    main()