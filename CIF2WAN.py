#MAIN
import sys, os
sys.path.append('Bin/')
#import input
from ESPRESSO import *
from ABINIT import *
from VASP import *
from SIESTA import *
import numpy as np
from datetime import datetime


# File reading
import re
delimiters = ' = | # |\n'

class CIF2WAN:
    def __init__(self,inputfile):
        input = open(inputfile,'r')
        for line in input:
            if "ncore" in line.lower():
                self.ncore = int(re.split(delimiters, line)[1])
                #print(self.ncore)
            if "ecut" in line.lower():
                self.ECUT = int(re.split(delimiters, line)[1])
            if "dft" in line.lower():
                self.DFT = str(re.split(delimiters, line)[1])
            if "inputtype" in line.lower():
                self.inputtype = str(re.split(delimiters, line)[1])
            if "apikey" in line.lower():
                self.APIKEY = str(re.split(delimiters, line)[1])
            if "materialid" in line.lower():
                self.MATERIALID = str(re.split(delimiters, line)[1])
            if "seedname" in line.lower():
                self.seedname = str(re.split(delimiters, line)[1])
            if "ciffile" in line.lower():
                self.CIF_file = str(re.split(delimiters, line)[1])
            if "email" in line.lower():
                self.Email = re.split(delimiters, line)[1]
            if "soc" in line.lower():
                self.SOC = bool(re.split(delimiters, line)[1])
            if "kmesh" in line.lower():
                self.KMESH = np.fromstring(re.split(delimiters, line)[1], dtype=int, sep=',')
                #print(self.KMESH)
            if "numbands" in line.lower():
                self.NUMBANDS = int(re.split(delimiters, line)[1])
            if "magnetism" in line.lower():
                self.Magnetism = bool(re.split(delimiters, line)[1])
            if "conventional_cell" in line.lower():
                self.conventional_cell = bool(re.split(delimiters, line)[1])
            if "relax" in line.lower():
                self.Relax = bool(re.split(delimiters, line)[1])
    def generate_log(self):
        with open("log.txt", 'w') as f:
            now = datetime.now()
            dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
            f.write("CIF2WAN V0.7 \n")
            f.write("CC Christopher Sims \n\n\n")
            f.write("Log File \n")
            f.write("Created on: " + dt_string + "\n")
            f.write("Generated input file for: " + C2W.DFT)
            #f.write
# if (input.DFT == 'ESPRESSO'):
#     if (input.InputType == 'Pymatgen'):
#         run = MAT2ESPRESSO()
#         run.PMG2WAN()
#     if (input.InputType == 'CIF'):
#         run = CIF2ESPRESSO(input.ciffile)
#         run.CIF2WAN()
# if (input.DFT == 'ABINIT'):
#     run = MAT2ABINIT()
#     run.PMG2WAN()
# if (input.DFT == 'VASP'):
#     run = MAT2VASP()
#     run.PMG2WAN()
# if (input.DFT == 'SIESTA'):
#     run = MAT2SIESTA()
#     run.PMG2WAN()

if __name__ == "__main__":
    C2W = CIF2WAN("input.in")
    #ESPRESSO # fully tested
    if (C2W.DFT.lower() == "espresso"):
        if (C2W.inputtype.lower() == "pymatgen"):
            run = CIF2ESPRESSO(C2W) # Tested
        if (C2W.inputtype.lower() == "cif"):
            run = PMG2ESPRESSO(C2W) # Tested
    #ABINIT # fully tested
    if (C2W.DFT.lower() == "abinit"):
        if (C2W.inputtype.lower() == "pymatgen"):
            run = PMG2ABINIT(C2W) # Tested
        if (C2W.inputtype.lower() == "cif"):
            run = CIF2ABINIT(C2W) # Tested
    #SIESTA # tested, may add for features to .fdf file
    if (C2W.DFT.lower() == "siesta"):
        if (C2W.inputtype.lower() == "pymatgen"):
            run = PMG2SIESTA(C2W) # Tested
        if (C2W.inputtype.lower() == "cif"):
          run = CIF2SIESTA(C2W) # tested
    if (C2W.DFT.lower() == "vasp"):
        if (C2W.inputtype.lower() == "pymatgen"):
            run = PMG2VASP(C2W) # Tested
        if (C2W.inputtype.lower() == "cif"):
            run = CIF2VASP(C2W) # NOT WELL TESTED/WORKING
    run.generate()


