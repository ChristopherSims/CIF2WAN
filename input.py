#INPUT FILE input.py
ncore = 50 # number of cores for calculation
ECUT = 60 # Energy cutoff in Ryberg
DFT = "ESPRESSO" #ESPRESSO, SIESTA, ABINIT, VASP
APIKEY = "4XK0y0eziIorGoHP0M" #APIKEY for mat
MATERIALID = "mp-20630" # input for pymatgen
SEEDNAME = "HMS" #your own name or default
EMAIL = "christopherssims95@gmail.com"
SOC = True # TRUE or FALSE
KMESH = [4,4,4] # KMESH for nscf/wan
NUMBANDS = 120 # Number of bands
MAGNETISM = False
conventional_cell = True # Select conventional cell #Ignore for pure CIF
Relax = False # Get the relaxed unit cell, ignore for pure CIF
