import numpy as np

class mol(object):
    pass

epoxy = mol()

epoxy.box_dimensions = [10,10,10]

epoxy.coords = [ 

epoxy.atom_labels =[

epoxy.vdw_defs = {}

epoxy.molecule_labels = [1] * len(epoxy.coords)

epoxy.bonds = np.array(([  

from connector import Connector
from write_coords import Writer
from params import Parameterise

connect = Connector()
epoxy.bond_types = connect.find_bond_types(epoxy.atom_labels,epoxy.bonds)
epoxy.bond_labels = connect.bond_labels(epoxy.atom_labels,epoxy.bonds,epoxy.bond_types)

epoxy.angles = connect.angles(epoxy.bonds)
epoxy.angle_types = connect.find_angle_types(epoxy.atom_labels,
                                           epoxy.angles)
epoxy.angle_labels = connect.angle_labels(epoxy.atom_labels,epoxy.angles,
                                        epoxy.angle_types)
epoxy.torsions = connect.torsions(epoxy.bonds)
epoxy.torsion_types = connect.find_torsion_types(epoxy.atom_labels,
                                           epoxy.torsions)
epoxy.torsion_labels = connect.torsion_labels(epoxy.atom_labels,epoxy.torsions,
                                        epoxy.torsion_types)
epoxy.impropers = connect.impropers(epoxy.bonds)
epoxy.improper_types = connect.find_improper_types(epoxy.atom_labels,
                                           epoxy.impropers)
epoxy.improper_labels = connect.improper_labels(epoxy.atom_labels,epoxy.impropers,
                                        epoxy.improper_types)

epoxy.masses = [0] * len(epoxy.coords)
#output.write_xyz('epoxy.xyz')

p = Parameterise(epoxy.vdw_defs)

epoxy.bond_coeffs = p.match_bonds(epoxy.bond_types)
epoxy.angle_coeffs = p.match_angles(epoxy.angle_types)
epoxy.torsion_coeffs = p.match_torsions(epoxy.torsion_types)
epoxy.improper_coeffs = p.match_impropers(epoxy.improper_types)
epoxy.pair_coeffs = p.match_pairs()
epoxy.masses = p.match_masses()
epoxy.charge_coeffs = p.match_charges()
epoxy.atom_charges = []

for atom in epoxy.atom_labels:
    index = epoxy.charge_coeffs['a'].index(atom)
    epoxy.atom_charges += [ epoxy.charge_coeffs['q'][index] ]


output = Writer(epoxy)
output.write_lammps('2.data')


#a = Writer(paa)
#a.write_lammps()


