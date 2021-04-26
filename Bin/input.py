# INPUT FILE input.py
##############
# This input file is read by Main.py in order to run the code.
##############
##############
# Ncore is forced to be 1 for Siesta due to a fatal error in the parallelization of the code
# Siesta runs very well on 1 core so there is no major issue.
##############
ncore = 50 # number of cores for calculation
ECUT = 60 # Energy cutoff in Ryberg
DFT = "ESPRESSO" # ESPRESSO, SIESTA, ABINIT, VASP
InputType = "Pymatgen " # Pymatgen or CIF
APIKEY = "qlGbypUfCtvT9INM" # APIKEY for Pymatgen
MATERIALID = "mp-1067587" # input for Pymatgen
SEEDNAME = "CIS" # your own name or default
ciffile = "D:\CIF2WANv0.5\CIFDIR\CIS.cif" # full path to CIF file
EMAIL = "christopherssims95@gmail.com"
SOC = False # TRUE or FALSE
KMESH = [4,4,4] # KMESH for nscf/wan
NUMBANDS = 120 # Number of bands
MAGNETISM = False # True or False; add magnetism
conventional_cell = False # True or False; Select conventional cell #Ignore for pure CIF
Relax = True # True or False; Get the relaxed unit cell, ignore for pure CIF
