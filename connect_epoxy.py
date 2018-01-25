import numpy as np

class mol(object):
    pass

epoxy = mol()

epoxy.box_dimensions = [10,10,10]

epoxy.coords = [ 
       [ -6.57985,  1.40495, -1.10509], #p-aromatic C   90
       [ -5.18100,  1.49216, -1.12543], #aromatic C     90
       [ -4.39814,  0.46700, -1.68126], #ar C           90
       [ -4.66028,  2.35152, -0.70818], #ar H           91
       [ -4.99952, -0.66058, -2.24614], #Phenol C       413
       [ -3.31596,  0.54998, -1.67249], #ar H           91
       [ -6.37964, -0.74993, -2.25213], #ar C           90
       [ -7.16092,  0.27500, -1.71781], #ar C           90
       [ -8.24049,  0.15422, -1.74998], #ar H           91
       [ -6.89970, -1.61960, -2.63218], #ar H           91
       [-10.06834,  2.00455, -0.24614], #ar C           90
       [ -8.74617,  1.80451,  0.19575], #p-ar C         90
       [ -9.56074, -0.08329,  1.53081], #ar C           90 
       [-10.83222,  0.05323,  0.99099], #phenol C       413
       [-11.10533,  1.12710,  0.14072], #ar C           90
       [-10.30420,  2.79294, -0.95819], #ar H           91
       [-12.09514,  1.26014, -0.28735], #ar H           91
       [ -7.50114,  2.48311, -0.45994], #tert alkyl C   84
       [ -6.77671,  3.30965,  0.63192], #methyl c       80
       [ -5.95503,  3.90567,  0.21572], #alkyl H        85
       [ -6.34875,  2.67536,  1.41635], #alkyl H        85
       [ -7.46390,  4.01556,  1.11359], #alkyl H        85
       [ -7.84700,  3.46503, -1.60652], #methyl c       80
       [ -8.39080,  2.97188, -2.42108], #alkyl H        85
       [ -6.93682,  3.89169, -2.05005], #alkyl H        85
       [ -8.45581,  4.30432, -1.25232], #alkyl H        85
       [ -4.15384, -1.63434, -2.70118], #phenyl O       414
       [ -4.59802, -2.38279, -3.84326], # ether C       124
       [ -5.30272, -3.63806, -3.43229], # CH - epoxy    125
       [ -3.68916, -2.66600, -4.38674], #alkyl H        85
       [ -5.18884, -1.76282, -4.52922], #alkyl H        85
       [ -6.72870, -3.92671, -3.84433], #CH2 - epoxy    124
       [ -5.61964, -4.51409, -4.52873], #epoxy O        122
       [ -4.89395, -4.09586, -2.54222], #HC epoxy       85
       [ -7.29432, -3.16866, -4.37257], #H2C epoxy      85
       [ -7.32500, -4.56215, -3.20092], #H2C epoxy      85
       [-11.70690, -0.98775,  1.20236], #phenyl O       414 
       [ -8.53610,  0.79194,  1.15705], #ar C           90
       [ -9.34881, -0.91413,  2.19730], #ar H           91
       [ -7.53590,  0.59885,  1.54020], #ar H           91
       [-12.75654, -0.64314,  2.10695], #ether C        124
       [-13.12117, -1.75269,  3.06321], #CH epoxy       125
       [-13.64088, -0.39108,  1.51036], #alkyl H        85
       [-12.49898,  0.23562,  2.71112], #alkyl H        85
       [-12.53765, -3.14211,  3.29264], #CH2 epoxy      124
       [-12.04473, -2.00403,  3.99269], #epoxy O        122
       [-14.11200, -1.58990,  3.47061], #HC epoxy       85 
       [-11.91546, -3.61935,  2.54934], #H2C epoxy      85
       [-13.10046, -3.84446,  3.90484]] #H2C epoxy      85

epoxy.atom_labels =[
        90,
        90,
        90,
        91,
        413,
        91,
        90,
        90,
        91,
        91,
        90,
        90,
        90 ,
        413,
        90,
        91,
        91,
        84,
        80,
        85,
        85,
        85,
        80,
        85,
        85,
        85,
        414,
        124,
        125,
        85,
        85,
        124,
        122,
        85,
        85,
        85,
        414,
        90,
        91,
        91,
        124,
        125,
        85,
        85,
        124,
        122,
        85 ,
        85,
        85]


epoxy.vdw_defs = {}
for i in set(epoxy.atom_labels):
    epoxy.vdw_defs[i] = i

print epoxy.vdw_defs


epoxy.molecule_labels = [1] * len(epoxy.coords)
epoxy.bonds = np.array(([  
[1 ,2 ],
[1 ,8 ],
[2 ,3 ],
[2 ,4 ],
[3 ,5 ],
[3 ,6 ],
[5 ,7 ],
[7 ,10], 
[8 ,7 ],
[8 ,9 ],
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
[18,1 ], 
[5 ,27],
[27,28],
[28,29],
[28,30],
[28,31],
[29,32],
[29,33],
[29,34],
[33,32],
[32,35],
[32,36],
[14,37],
[38,12],
[38,13],
[13,39],
[38,40],
[37,41],
[41,42],
[41,43],
[41,44],
[42,45],
[42,46],
[42,47],
[45,46],
[45,48],
[45,49]]))


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


