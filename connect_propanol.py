import numpy as np

class mol(object):
    pass

propanol = mol()

propanol.box_dimensions = [10,10,10]

propanol.coords = [ 
[ -3.84584 ,      2.15936 ,     -0.21782], # C -Me 80   1
[ -2.47805 ,      2.21806 ,      0.46176], # C - O 100  2
[ -3.90609 ,      1.46855 ,     -1.06341], # H     89   3
[ -4.57868 ,      1.80371 ,      0.52101], # H     89   3
[ -1.29439 ,      1.51962 ,     -0.20421], #       80   1
[ -0.82426 ,      2.11761 ,     -0.99282], #       89   3
[ -1.53347 ,      0.52539 ,     -0.59142], #       89   3
[ -2.13223 ,      3.58802 ,      0.68331], #       96   4
[ -2.59747 ,      1.76017 ,      1.44653], #       89   3
[ -1.62465 ,      3.90104 ,     -0.07559], #       97   5
[ -4.20545 ,      3.14774 ,     -0.52766], #       89   3
[ -0.51525 ,      1.39178 ,      0.56033]] #       89   3



propanol.atom_labels =[1,2,3,3,1,3,3,4,3,5,3,3]
#propanol.atom_labels =[80,100,89,89,80,89,89,96,89,97,89,89]

propanol.vdw_defs = {1:80, 2:100, 3:85, 4:96, 5:97}
#propanol.vdw_defs = {80:80, 100:100, 89:89, 96:96, 97:97}

propanol.molecule_labels = [1] * len(propanol.coords)

propanol.bonds = np.array(([  
[1, 2],
[1, 3],
[1, 4],
[2, 5],
[5, 6],
[5, 7],
[2, 8],
[2, 9],
[8, 10],
[1, 11],
[5, 12]]))





from connector import Connector
from write_coords import Writer
from params import Parameterise

connect = Connector()
propanol.bond_types = connect.find_bond_types(propanol.atom_labels,propanol.bonds)
propanol.bond_labels = connect.bond_labels(propanol.atom_labels,propanol.bonds,propanol.bond_types)

propanol.angles = connect.angles(propanol.bonds)
propanol.angle_types = connect.find_angle_types(propanol.atom_labels,
                                           propanol.angles)
propanol.angle_labels = connect.angle_labels(propanol.atom_labels,propanol.angles,
                                        propanol.angle_types)
propanol.torsions = connect.torsions(propanol.bonds)
propanol.torsion_types = connect.find_torsion_types(propanol.atom_labels,
                                           propanol.torsions)
propanol.torsion_labels = connect.torsion_labels(propanol.atom_labels,propanol.torsions,
                                        propanol.torsion_types)
propanol.impropers = connect.impropers(propanol.bonds)
propanol.improper_types = connect.find_improper_types(propanol.atom_labels,
                                           propanol.impropers)
propanol.improper_labels = connect.improper_labels(propanol.atom_labels,propanol.impropers,
                                        propanol.improper_types)

propanol.masses = [0] * len(propanol.coords)
#output.write_xyz('propanol.xyz')

p = Parameterise(propanol.vdw_defs)

propanol.bond_coeffs = p.match_bonds(propanol.bond_types)
propanol.angle_coeffs = p.match_angles(propanol.angle_types)
propanol.torsion_coeffs = p.match_torsions(propanol.torsion_types)
propanol.improper_coeffs = p.match_impropers(propanol.improper_types)
propanol.pair_coeffs = p.match_pairs()
propanol.masses = p.match_masses()
propanol.charge_coeffs = p.match_charges()
propanol.atom_charges = []

for atom in propanol.atom_labels:
    index = propanol.charge_coeffs['a'].index(atom)
    propanol.atom_charges += [ propanol.charge_coeffs['q'][index] ]


output = Writer(propanol)
output.write_lammps('2.data')


#a = Writer(paa)
#a.write_lammps()


