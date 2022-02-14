"""
Python script to generate model MOFs

Requires MeXn_base.xyz and ligand_base.xyz in same folder

Usage:
$ python mof_builder.py mof_write_name.xyz [options]

e.g.
$ python mof_builder.py ZnO4_antroquinone.xyz --aq --metal Zn

"""
from ase.io import read
import argparse, sys, os

# Indices of changable atoms:
# Oxygens - 1, 2, 3, 4
# Metal - 0
# AQ model H's - 21, 22
# N for pyridine position - 5

# Argument parsing
parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('outputname', help='output filename .xyz')

parser.add_argument('--x_stretch', type=float, help='Cell stretch from default in X dir', default = 0)
parser.add_argument('--vac_y', type=float, help='Vacuum in y direction', default = 14)
parser.add_argument('--vac_z', type=float, help='Vacuum in z direction', default = 14)

parser.add_argument('--metal', help='Metal center', default = 'Ni')
parser.add_argument('--menx_xelem', help='Elements around metal, default O', default = 'O')
parser.add_argument('--pyridine', help='Change linker to pyridine type', action='store_true')
parser.add_argument('--halfaq', help='Change linker to half-antroquinone', action='store_true')
parser.add_argument('--aq', help='Change linker to antroquinone', action='store_true')

parser.add_argument('--nlayer', type=int, help='If >1 write slab instead of single layer; uses vac_z as spacing',
                    default = 1)
args=parser.parse_args()

# Fragment bases and default associated data
mexn = read("MeXn_base.xyz")
ligand = read("ligand_base.xyz")
Me_X_xdirlen = 1.4
X_linker_xdirlen = 1.2
cell_len_xdir = 10 + Me_X_xdirlen + X_linker_xdirlen
for atom in ligand:
    atom.position = atom.position+[Me_X_xdirlen+X_linker_xdirlen, -2.5, 0]
mof = mexn+ligand
mof.pbc = [1, 1, 1]
vacuum = 14
mof.cell = [cell_len_xdir+args.x_stretch,5+args.vac_y,args.vac_z*args.nlayer]

# Change metal center as specified in arguments
mof[0].symbol = args.metal
mof[1].symbol = args.menx_xelem
mof[2].symbol = args.menx_xelem
mof[3].symbol = args.menx_xelem
mof[4].symbol = args.menx_xelem

# Change linker as specified in arguments
if args.halfaq:
    mof[21].symbol = 'O'
    mof[21].position += [0, -0.3, 0]
    if args.aq:
        print("[-] Incompatible settings --aq and --halfaq! Aborting....")
        sys.exit(-1)
if args.aq:
    mof[21].symbol = 'O'
    mof[21].position += [0, -0.3, 0]
    mof[22].symbol = 'O'
    mof[22].position += [0, 0.3, 0]
if args.pyridine:
    mof[15].symbol = 'N'
    del mof[23]

if args.nlayer > 1:
    mof_base = mof.copy()
    for i in range(args.nlayer):
        mof_2 = mof_base.copy()
        for atom in mof_2:
            atom.position += [0, 0, i*args.vac_z]
        mof += mof_2

mof.center()

mof.write(args.outputname)



