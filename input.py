#INPUT FILE input.py
ncore = 50 # number of cores for calculation
ECUT = 60 # Energy cutoff in Ryberg
DFT = "ESPRESSO" #ESPRESSO, SIESTA, ABINIT, VASP
APIKEY = "YourAPIKEY" #APIKEY for mat
MATERIALID = "mp-20630" # input for pymatgen
SEEDNAME = "HMS" #your own name or default
EMAIL = "yyy@mail.com"
SOC = True # TRUE or FALSE
KMESH = [4,4,4] # KMESH for nscf/wan
NUMBANDS = 120 # Number of bands
MAGNETISM = False
conventional_cell = True # Select conventional cell #Ignore for pure CIF
Relax = False # Get the relaxed unit cell, ignore for pure CIF
