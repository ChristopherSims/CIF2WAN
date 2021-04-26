from ase.io import read, write
from ase.optimize import BFGS
from ase.visualize import view
from ase.io import write as ASEwrite

import os
parentfolder = os.path.dirname(os.path.dirname( __file__ ))


crystal = read(parentfolder + '\CIFDIR\CIS2.cif')
#view(crystal)
#ASEwrite('image.png', crystal)
#opt = BFGS(crystal, trajectory='opt.traj', logfile='opt.log')
#opt.run(fmax=0.05)
print("Chemical formula \n")
print(crystal.get_chemical_formula())
print('Unit cell \n')
print(repr(crystal.get_cell()))
print("atom positions \n")
print(repr(crystal.get_scaled_positions()))
print(len(set(crystal.numbers)))
print(crystal.get_chemical_symbols())


#Base = Crystal.from_cif(parentfolder + '\CIFDIR\CSE2.cif') 
