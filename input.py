#INPUT FILE input.py
ncore = 50 # number of cores for calculation
ECUT = 60 # Energy cutoff in Ryberg
DFT = "ESPRESSO" #ESPRESSO, SIESTA, ABINIT, VASP
APIKEY = "4XK0y0eziIorGoHP0M" #APIKEY for mat
MATERIALID = "mp-1523" # imput for pymatgen
SEEDNAME = "ZP" #your own name or default
EMAIL = "Christopherssims95@gmail.com"
SOC = True # TRUE or FALSE
KMESH = [4,4,4] # KMESH for nscf/wan
NUMBANDS = 120 # Number of bands
MAGNETISM = False
conventional_cell = True
