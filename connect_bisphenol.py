import numpy as np

class mol(object):
    pass

bp = mol()

bp.box_dimensions = [10,10,10]

bp.coords = [
[ -6.57985     ,  1.40495  ,    -1.10509], # CA - p 90  1
[ -5.18100     ,  1.49216  ,    -1.12543], # CA 90
[ -4.39814     ,  0.46700  ,    -1.68126], # CA 90
[ -4.66028     ,  2.35152  ,    -0.70818], # HA 91      2
[ -4.99952     , -0.66058  ,    -2.24614], # C-O 108    3
[ -3.31596     ,  0.54998  ,    -1.67249], # HA 91
[ -6.37964     , -0.74993  ,    -2.25213], # CA 90
[ -7.16092     ,  0.27500  ,    -1.71781], # CA 90
[ -8.24049     ,  0.15422  ,    -1.74998], # HA 91
[ -6.89970     , -1.61960  ,    -2.63218], # HA 91
[-10.06834     ,  2.00455  ,    -0.24614], # CA 90
[ -8.74617     ,  1.80451  ,     0.19575], # CA - p 90
[ -9.56074     , -0.08329  ,     1.53081], # CA 90 
[-10.83222     ,  0.05323  ,     0.99099], # C - O 108  3
[-11.10533     ,  1.12710  ,     0.14072], # CA 90 
[-10.30420     ,  2.79294  ,    -0.95819], # HA 91
[-12.09514     ,  1.26014  ,    -0.28735], # HA 91
[ -7.50114     ,  2.48311  ,    -0.45994], # C -tert 84 4
[ -6.77671     ,  3.30965  ,     0.63192], # C - Me 80  5
[ -5.95503     ,  3.90567  ,     0.21572], # H 85       6
[ -6.34875     ,  2.67536  ,     1.41635], # H 85
[ -7.46390     ,  4.01556  ,     1.11359], # H 85
[ -7.84700     ,  3.46503  ,    -1.60652], # C - Me 80
[ -8.39080     ,  2.97188  ,    -2.42108], # H 85
[ -6.93682     ,  3.89169  ,    -2.05005], # H 85
[ -8.45581     ,  4.30432  ,    -1.25232], # H 85
[ -4.15384     , -1.63434  ,    -2.70118], # O 109      7
[-11.70690     , -0.98775  ,     1.20236], # O 109
[ -8.53610     ,  0.79194  ,     1.15705], # CA 90
[ -9.34881     , -0.91413  ,     2.19730], # HA 91
[ -7.53590     ,  0.59885  ,     1.54020], # HA 91
[-11.47410     , -1.68454  ,     1.78397], # HO 110     8
[ -3.22667     , -1.50372  ,    -2.66567]] # HO 110

bp.atom_labels =[1,1,1,2,3,2,1,1,2,2,1,1,1,3,1,2,2,4,5,6,6,6,5,6,6,6,7,7,1,2,2,8,8]

bp.vdw_defs = {1:90, 
               2:91,
               3:108,
               4:84,
               5:80,
               6:85,
               7:109,
               8:110}

bp.molecule_labels = [1] * len(bp.coords)

bp.bonds = np.array(([  
[1 ,2], 
[1 ,8],
[2 ,3], 
[2 ,4], 
[3 ,5], 
[3 ,6], 
[5 ,7], 
[7 ,10],
[8 ,7], 
[8 ,9], 
[11,12],
[11,15],
[11,16],
[13,14],
[15,14],
[15,17],
[12,18],
[18,19],
[19,20],
[19,21],
[19,22],
[18,23],
[23,24],
[23,25],
[23,26],
[18,1],
[5 ,27],
[14,28],
[29,12],
[29,13],
[13,30],
[29,31],
[27,33],
[28,32]]))





from connector import Connector
from write_coords import Writer
from params import Parameterise

connect = Connector()
bp.bond_types = connect.find_bond_types(bp.atom_labels,bp.bonds)
bp.bond_labels = connect.bond_labels(bp.atom_labels,bp.bonds,bp.bond_types)

bp.angles = connect.angles(bp.bonds)
bp.angle_types = connect.find_angle_types(bp.atom_labels,
                                           bp.angles)
bp.angle_labels = connect.angle_labels(bp.atom_labels,bp.angles,
                                        bp.angle_types)
bp.torsions = connect.torsions(bp.bonds)
bp.torsion_types = connect.find_torsion_types(bp.atom_labels,
                                           bp.torsions)
bp.torsion_labels = connect.torsion_labels(bp.atom_labels,bp.torsions,
                                        bp.torsion_types)
bp.impropers = connect.impropers(bp.bonds)
bp.improper_types = connect.find_improper_types(bp.atom_labels,
                                           bp.impropers)
bp.improper_labels = connect.improper_labels(bp.atom_labels,bp.impropers,
                                        bp.improper_types)

bp.masses = [0] * len(bp.coords)
#output.write_xyz('bp.xyz')

p = Parameterise(bp.vdw_defs)

bp.bond_coeffs = p.match_bonds(bp.bond_types)
bp.angle_coeffs = p.match_angles(bp.angle_types)
bp.torsion_coeffs = p.match_torsions(bp.torsion_types)
bp.improper_coeffs = p.match_impropers(bp.improper_types)
bp.pair_coeffs = p.match_pairs()
bp.masses = p.match_masses()
bp.charge_coeffs = p.match_charges()
bp.atom_charges = []

for atom in bp.atom_labels:
    index = bp.charge_coeffs['a'].index(atom)
    bp.atom_charges += [ bp.charge_coeffs['q'][index] ]


output = Writer(bp)
output.write_lammps('2.data')
output.write_xyz('2.xyz')

#a = Writer(paa)
#a.write_lammps()


