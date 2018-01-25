import numpy as np

class mol(object):
    pass

diamine = mol()

diamine.box_dimensions = [10,10,10]

diamine.coords = [
[-4.46220,       2.06632,      -0.28335],# N        1
[-3.33502,       1.28236,       0.26355],# C - N    2
[-3.51557,      -0.16282,      -0.19359],# C - Me   3
[-4.32267,      -0.63139,       0.38671],# H - Me   4
[-2.61242,      -0.75834,      -0.04337],# H - Me   4
[-3.80534,      -0.23082,      -1.24848],# H - Me   4
[-1.99319,       1.95161,      -0.13256],# C - O    5
[-3.40889,       1.28761,       1.35709],# H        4
[-0.77796,       1.25235,       0.23775],# O        6
[-1.97294,       2.92903,       0.36918],# H        7
[-1.97642,       2.12958,      -1.21515],# H        7
[ 0.34622,       2.10520,      -0.07618],# C - O    5
[ 1.77756,       1.56356,       0.19576],# C - N    2
[ 0.23849,       3.01922,       0.52320],# H        7
[ 0.27063,       2.39168,      -1.13304],# H        7
[ 2.03318,       0.11394,      -0.20380],# C - Me   3
[ 1.60051,      -0.56968,       0.53309],# H        4
[ 3.11140,      -0.10215,      -0.22416],# H        4
[ 1.63094,      -0.12372,      -1.19204],# H        4
[ 2.75124,       2.42686,      -0.50778],# N        1
[ 1.98488,       1.64063,       1.26995],# H        4
[-4.40957,       2.07236,      -1.30197],# HN       8
[-4.41148,       3.02351,       0.06541],# HN       8
[ 2.58510,       2.37768,      -1.51296],# HN       8
[ 3.69827,       2.10768,      -0.30368]]# HN       8

#diamine.atom_labels =[730,744,80,85,85,85,125,85,122,127,127,125,744,127,127,80,85,85,85,730,85,739,739,739,739]
diamine.atom_labels =[1,2,3,4,4,4,5,4,6,7,7,5,2,7,7,3,4,4,4,1,4,8,8,8,8]

diamine.vdw_defs = {1: 730,# N
                    2: 742,# C - N
                    3: 80,# C Me
                    4: 85,# H alkyl
                    5: 124,# C - O
                    6: 122,# O
                    7: 127,# H - CO
                    8: 739# H - N
                    } 

diamine.molecule_labels = [1] * len(diamine.coords)

diamine.bonds = np.array(([  
[1  ,2], 
[2  ,3], 
[3  ,4], 
[3  ,5], 
[3  ,6], 
[2  ,7], 
[2  ,8], 
[7  ,9], 
[7  ,10],
[7  ,11],
[9  ,12],
[12 ,13],
[12 ,14],
[12 ,15],
[13 ,16],
[16 ,17],
[16 ,18],
[16 ,19],
[13 ,20],
[13 ,21],
[1  ,22],
[1  ,23],
[20 ,24],
[20 ,25]]))


from connector import Connector
from write_coords import Writer
from params import Parameterise

connect = Connector()
diamine.bond_types = connect.find_bond_types(diamine.atom_labels,diamine.bonds)
diamine.bond_labels = connect.bond_labels(diamine.atom_labels,diamine.bonds,diamine.bond_types)

diamine.angles = connect.angles(diamine.bonds)
diamine.angle_types = connect.find_angle_types(diamine.atom_labels,
                                           diamine.angles)
diamine.angle_labels = connect.angle_labels(diamine.atom_labels,diamine.angles,
                                        diamine.angle_types)
diamine.torsions = connect.torsions(diamine.bonds)
diamine.torsion_types = connect.find_torsion_types(diamine.atom_labels,
                                           diamine.torsions)
diamine.torsion_labels = connect.torsion_labels(diamine.atom_labels,diamine.torsions,
                                        diamine.torsion_types)
diamine.impropers = connect.impropers(diamine.bonds)
diamine.improper_types = connect.find_improper_types(diamine.atom_labels,
                                           diamine.impropers)
diamine.improper_labels = connect.improper_labels(diamine.atom_labels,diamine.impropers,
                                        diamine.improper_types)

diamine.masses = [0] * len(diamine.coords)
#output.write_xyz('diamine.xyz')

p = Parameterise(diamine.vdw_defs)

diamine.bond_coeffs = p.match_bonds(diamine.bond_types)
diamine.angle_coeffs = p.match_angles(diamine.angle_types)
diamine.torsion_coeffs = p.match_torsions(diamine.torsion_types)
diamine.improper_coeffs = p.match_impropers(diamine.improper_types)
diamine.pair_coeffs = p.match_pairs()
diamine.masses = p.match_masses()
diamine.charge_coeffs = p.match_charges()
diamine.atom_charges = []

for atom in diamine.atom_labels:
    index = diamine.charge_coeffs['a'].index(atom)
    diamine.atom_charges += [ diamine.charge_coeffs['q'][index] ]

print diamine.impropers
output = Writer(diamine)
output.write_lammps('2.data')


#a = Writer(paa)
#a.write_lammps()


