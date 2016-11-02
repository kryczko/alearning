#!/usr/bin/env python

from ase.lattice.cubic import FaceCenteredCubic as fcc
from ase.lattice.cubic import BodyCenteredCubic as bcc
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description='Generate FCC Al lattices of user defined sizes')
    parser.add_argument('--xyz', '-x', help='Generate xyz format, default is POSCAR', action='store_true')
    parser.add_argument('--size', '-s', help='Size of lattice to generate, default is unit cell', dest='size', nargs=3, default=[1,1,1])
    return parser.parse_args()

def main():
    args = parseArgs()
    args.size = tuple([int(i) for i in args.size])
    fcc_atoms = fcc(directions=[[1,0,0], [0,1,0], [0,0,1]], size=args.size, symbol='Al', pbc=(1,1,1), latticeconstant=2.856)
    if args.xyz:
        f = open('al_fcc_atoms.xyz', 'w')
        f.write(str(len(fcc_atoms)) + "\n\n")
        for atom in fcc_atoms:
            line = "Al  " + "  ".join([str(val) for val in atom.position])
            f.write(line + '\n')
        f.close()

    else:
        f = open('POSCAR_fcc_Al', 'w')
        f.write('Alsys\n1.0\n')
        for vec in fcc_atoms.get_cell():
            f.write("  ".join([str(val) for val in vec]) + "\n")
        f.write('  Al\n')
        f.write('  ' + str(len(fcc_atoms)) + '\n')
        f.write('Cartesian\n')
        for atom in fcc_atoms:
            f.write("  ".join([str(val) for val in atom.position]) + '\n')
        f.close()

    bcc_atoms = bcc(directions=[[1,0,0], [0,1,0], [0,0,1]], size=args.size, symbol='Al', pbc=(1,1,1), latticeconstant=2.789)
    if args.xyz:
        f = open('al_bcc_atoms.xyz', 'w')
        f.write(str(len(bcc_atoms)) + "\n\n")
        for atom in bcc_atoms:
            line = "Al  " + "  ".join([str(val) for val in atom.position])
            f.write(line + '\n')
        f.close()

    else:
        f = open('POSCAR_bcc_Al', 'w')
        f.write('Alsys\n1.0\n')
        for vec in bcc_atoms.get_cell():
            f.write("  ".join([str(val) for val in vec]) + "\n")
        f.write('  Al\n')
        f.write('  ' + str(len(bcc_atoms)) + '\n')
        f.write('Cartesian\n')
        for atom in bcc_atoms:
            f.write("  ".join([str(val) for val in atom.position]) + '\n')
        f.close()

if __name__ == '__main__':
    main()