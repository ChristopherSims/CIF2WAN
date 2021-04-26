### VASP ###
import sys, os
sys.path
parentfolder = os.path.dirname(os.path.dirname( __file__ ))
sys.path.append(parentfolder)
#### PyMATGEN
import pymatgen as mg
from pymatgen.core import Composition
from pymatgen.core import Lattice, Structure, Molecule, Species
from pymatgen.core.periodic_table import Element
from pymatgen.io.cif import CifParser,CifWriter
import pymatgen.ext.matproj as matproj
# PLOTTING
import matplotlib.pyplot as plt
import numpy as np
#FILE reading
from xml.dom import minidom # READ XML
from string import digits
from shutil import copyfile
import re
import glob as gb
import csv
#import input
import kmesh
import pandas as pd

class PMG2VASP:
    def __init__(self,input):
        self.CIFDIR = "cifdir/"
        self.key = input.APIKEY # Import API key from input
        self.MaterialID = input.MATERIALID # Set class variable MaterialID from input
        self.seedname = str(input.seedname) # set class variable SEEDNAME from input file
        self.Ecut = input.ECUT #Class variable; in Ryberg, transfer from input file
        self.kpts = input.KMESH ### KPTS from input file for a generator, class variable
        self.SEEDDIR = "OUTPUT/" + 'V-' + self.seedname + '/'  #set the seed directory
        self.EMAIL = input.Email #Class variable; set email for SLURM from input file
        self.NCORE = input.ncore #Class Variable; number of cores, taken from input file
        self.SOC = input.SOC
        self.NBND = input.NUMBANDS
        selectconventional = input.conventional_cell
        selectrelax = input.Relax
        # MAKE DIRECTORIES
        if not os.path.exists(self.SEEDDIR):
            os.makedirs(self.SEEDDIR)
        if not os.path.exists(self.SEEDDIR + "PP"):
            os.makedirs(self.SEEDDIR + "PP")
        if not os.path.exists(self.SEEDDIR + "WT"):
            os.makedirs(self.SEEDDIR + "WT")
        if not os.path.exists(self.SEEDDIR + 'PBE'):
            os.makedirs(self.SEEDDIR + "PBE")
        if not os.path.exists(self.SEEDDIR + "HSE06"):
            os.makedirs(self.SEEDDIR + "HSE06")
        if not os.path.exists(self.SEEDDIR + "W90"):
            os.makedirs(self.SEEDDIR + "W90")
        #import structure from PYMATGEN
        with matproj.MPRester(self.key) as m:
            struct = m.get_structure_by_material_id(self.MaterialID,final=selectrelax,conventional_unit_cell=selectconventional) # GET STRUCT FROM MATERIALS PROJECT
            self.struct = m.get_structure_by_material_id(self.MaterialID,final=selectrelax,conventional_unit_cell=selectconventional) # GET STRUCT FROM MATERIALS PROJECT
        C = Composition(str(struct.formula))
        self.lattice = struct.lattice
        self.a =  self.lattice.a
        self.b = self.lattice.b
        self.c =  self.lattice.c
        self.alpha = self.lattice.alpha
        self.beta = self.lattice.beta
        self.gamma = self.lattice.gamma
        self.UNIQUE_ATOMS = np.unique(struct.species)
        self.cord = struct.frac_coords
        self.ntypat = np.unique(struct.atomic_numbers).size
        self.natom = np.array(struct.atomic_numbers).size
        self.atomnum = np.unique(struct.atomic_numbers)[::1]
        self.typat = ""
        for x in struct.atomic_numbers:
            for at in self.atomnum:
                if at == x:
                        self.typat += ''.join(map(str,np.where(self.atomnum == at)
                            [0] + 1)) + " "
        # SOC CASE
                # SOC CASE
        self.SOC = input.SOC
        self.noncolin = '.FALSE.'
        self.lspinorb = '.FALSE.'
        if input.SOC:
            self.noncolin = '.TRUE.'
            self.lspinorb = '.TRUE.'
        num_bands = input.NUMBANDS
        num_wann = num_bands - (num_bands%2) # FORCE EVEN
        self.wan = num_wann
        self.WTwan = num_wann
        if self.SOC:
            self.WTwan = num_wann/2
        num_bands = 0
        #for x in struct.species:
           # PP = gb.glob(self.pseudodir + str(x) + ".*") 
            #readpp = minidom.parse(''.join(map(str,PP)))
            #items = readpp.getElementsByTagName('atom')
            #num_bands += float(items[0].attributes['valence'].value)
        num_wann = num_bands - (num_bands%2) # FORCE EVEN
    def generate(self):
        with open(self.SEEDDIR + 'pbe/' + "INCAR",'w',newline='\n') as f:
            f.write("ALGO = Fast\n")
            f.write("EDIFF = 1.0e-8\n")
            f.write("ENCUT = 500\n")
            f.write("ISIF = 2\n")
            f.write("ISMEAR = -5\n")
            f.write("KPAR = 2\n")
            f.write("LCHARG = .FALSE.\n")
            f.write("LREAL = .FALSE.\n")
            f.write("NBANDS = %i\n" %(self.NBND))
            f.write("NEDOS = 2000\n")
            f.write("NPAR = 1\n")
            f.write("PREC = Accurate\n")
            f.write("SYSTEM = %s\n" %(self.seedname))
        with open(self.SEEDDIR + 'pbe/' + "KPOINTS",'w',newline='\n') as f:
            f.write("Automatic %ix%ix%i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("0\n")
            f.write("Gamma\n")
            f.write("\t%i %i %i"%(self.kpts[0],self.kpts[1],self.kpts[2]))
        with open(self.SEEDDIR + 'pbe/' + "POSCAR",'w',newline='\n') as f:
            f.write("%s\n" %(self.seedname))
            f.write("1\n")
            f.write(str(self.struct.lattice) + "\n")
            f.write("\t" + ' '.join(map(str,self.UNIQUE_ATOMS)) + "\n")
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.struct.atomic_numbers == self.atomnum[i])).size
            f.write("\t" + ' '.join(map(str,typnum)) + "\n")
            f.write("Direct\n")
            for index,x in enumerate(self.struct.frac_coords):
                f.write("%f %f %f\n" %(x[0], x[1], x[2]))
    #def write_hse(self):
        with open(self.SEEDDIR + 'HSE06/' + "INCAR",'w',newline='\n') as f:
            f.write("ALGO = Fast\n")
            f.write("EDIFF = 1.0e-8\n")
            f.write("ENCUT = 500\n")
            f.write("ISIF = 2\n")
            f.write("ISMEAR = -5\n")
            f.write("KPAR = 2\n")
            f.write("LCHARG = .FALSE.\n")
            f.write("LREAL = .FALSE.\n")
            f.write("NBANDS = %i\n" %(self.NBND))
            f.write("NEDOS = 2000\n")
            f.write("NPAR = 1\n")
            f.write("PREC = Accurate\n")
            f.write("SYSTEM = %s\n" %(self.seedname))
            f.write("LHFCALC = .TRUE.\n")
            f.write("HFSCREEN = 0.2\n")
        with open(self.SEEDDIR + 'HSE06/' + "KPOINTS",'w',newline='\n') as f:
            f.write("Automatic %ix%ix%i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("0\n")
            f.write("Gamma\n")
            f.write("\t%i %i %i"%(self.kpts[0],self.kpts[1],self.kpts[2]))
        with open(self.SEEDDIR + 'HSE06/' + "POSCAR",'w',newline='\n') as f:
            f.write("%s\n" %(self.seedname))
            f.write("1\n")
            f.write(str(self.struct.lattice) + "\n")
            f.write("\t" + ' '.join(map(str,self.UNIQUE_ATOMS)) + "\n")
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.struct.atomic_numbers == self.atomnum[i])).size
            f.write("\t" + ' '.join(map(str,typnum)) + "\n")
            f.write("Direct\n")
            for index,x in enumerate(self.struct.frac_coords):
                f.write("%f %f %f\n" %(x[0], x[1], x[2]))
    #def write_w90(self):
        with open(self.SEEDDIR + 'W90/' + "INCAR",'w',newline='\n') as f:
            f.write("ALGO = Fast\n")
            f.write("EDIFF = 1.0e-8\n")
            f.write("ENCUT = 500\n")
            f.write("ISIF = 2\n")
            f.write("ISMEAR = -5\n")
            f.write("KPAR = 2\n")
            f.write("LCHARG = .FALSE.\n")
            f.write("LREAL = .FALSE.\n")
            f.write("NBANDS = %i\n" %(self.NBND))
            f.write("NEDOS = 2000\n")
            f.write("NPAR = 1\n")
            f.write("PREC = Accurate\n")
            f.write("SYSTEM = %s\n" %(self.seedname))
            f.write("LHFCALC = .TRUE.\n")
            f.write("HFSCREEN = 0.2\n")
            f.write("LWANNIER90 = .TRUE.\n")
        with open(self.SEEDDIR + 'W90/' + "KPOINTS",'w',newline='\n') as f:
            f.write("Automatic %ix%ix%i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("0\n")
            f.write("Gamma\n")
            f.write("\t%i %i %i"%(self.kpts[0],self.kpts[1],self.kpts[2]))
        with open(self.SEEDDIR + 'W90/' + "POSCAR",'w',newline='\n') as f:
            f.write("%s\n" %(self.seedname))
            f.write("1\n")
            f.write(str(self.struct.lattice) + "\n")
            f.write("\t" + ' '.join(map(str,self.UNIQUE_ATOMS)) + "\n")
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.struct.atomic_numbers == self.atomnum[i])).size
            f.write("\t" + ' '.join(map(str,typnum)) + "\n")
            f.write("Direct\n")
            for index,x in enumerate(self.struct.frac_coords):
                f.write("%f %f %f\n" %(x[0], x[1], x[2]))
        with open(self.SEEDDIR + 'W90/wannier90.win','w',newline='\n') as f:
            f.write("!write_hr = .TRUE.\n")
            f.write("!write_xyz = .TRUE.\n")
            f.write("guiding_centres= .TRUE.\n")
            f.write("!wannier_plot = .TRUE. \n")
            f.write("spinors = %s\n" %(self.lspinorb))
            f.write("num_wann = %i\n" %(self.wan))
            f.write("dis_num_iter=1000\n")
            f.write("num_iter = 2000\n\n\n\n")
            f.write("begin unit_cell_cart\n")
            f.write(str(self.struct.lattice) + "\n") ### atomic structure
            f.write("end unit_cell_cart\n\n\n")
            f.write("begin atoms_frac\n")
            for index,x in enumerate(self.struct.frac_coords):
                f.write("%s %f %f %f\n" %(str(self.struct.species[index]),x[0], x[1], x[2]))
            f.write("end atoms_frac\n\n\n")
            f.write("begin projections \n")
            f.write("random \n")
            f.write("end projections\n\n\n")
            f.write(kmesh.WANNIER(self.kpts[0],self.kpts[1],self.kpts[2]))

    #def write_wt(self):
        with open(self.SEEDDIR + "WT/wt.in",'w',newline='\n') as f:
            f.write("#### wt file ###########\n")
            f.write("&TB_FILE\n")
            f.write("Hrfile = '%s_hr.dat'\n" %(self.seedname))
            f.write("Package = 'ESPRESSO'\n")
            f.write("/\n\n\n")
            f.write("LATTICE\n")
            f.write("Angstrom\n")
            f.write(str(self.struct.lattice) + "\n\n\n") ### atomic structure
            f.write("ATOM_POSITIONS\n")
            f.write("%i !Number of atoms for projectors\n" %(len(self.struct.frac_coords)))
            f.write("Direct ! Direct or Cartisen coordinateS\n")
            for index,x in enumerate(self.struct.frac_coords):
                f.write("%s %f %f %f\n" %(str(self.struct.species[index]),x[0], x[1], x[2]))
            f.write("\n\n\n\n")
            f.write("PROJECTORS\n")
            #### Places prjectors and stoms into PRJCARD in wannier tools
            proj = int(self.WTwan/self.natom)
            proj_rem = self.WTwan % self.natom
            sumproj = 0
            PROJPARAM = ""
            projector = ["s", "pz","px", "py"]
            LATPARAM = str(self.struct.lattice) + "\n" ### atomic structur
            PROJCARD =""
            ATOMCARD = ""
            appen = 0 #### append number
            for index,atom in enumerate(self.struct.species):
                if proj_rem > 0:
                    appen = 1
                proj_rem = proj_rem - 1
                PROJPARAM += str(proj + appen) + " "
                PROJCARD += str(atom) + " " + ' '.join(map(str,projector)) + "\n"
                #print(str(atom) + str(proj + appen))
                sumproj += proj + appen
                #print(sumproj)
                appen = 0
            f.write(PROJPARAM + "\n")
            f.write(PROJCARD)
            f.write("\n\n\n")

            f.write("&CONTROL\n")
            f.write("! BULK BAND CALCULATIONS \n")
            f.write("BulkBand_calc       =  F\n")
            f.write("BulkBand_plane_calc =  F \n")
            f.write("BulkFS_calc         =  F \n")
            f.write("BulkFS_Plane_calc   =  F\n")
            f.write("SlabBand_calc       =  F\n")
            f.write("Dos_calc            =  F\n")
            f.write("! BULK GAP\n")
            f.write("BulkGap_cube_calc   =  F\n")
            f.write("BulkGap_plane_calc  =  F\n")
            f.write("! SURFACE STATES\n")
            f.write("SlabSS_calc         =  T\n")
            f.write("SlabArc_calc        =  T\n")
            f.write("SlabSpintexture_calc =  F\n")
            f.write("! TOPO INV\n")
            f.write("wanniercenter_calc   = F\n")
            f.write("BerryPhase_calc     =  F\n")
            f.write("BerryCurvature_calc =  F\n")
            f.write("BerryCurvature_slab_calc =  F\n")
            f.write("Z2_3D_calc          =  F\n")
            f.write("WeylChirality_calc  =  F\n")
            f.write("NLChirality_calc    =  F\n")
            f.write("Chern_3D_calc       =  F\n")
            f.write("MirrorChern_calc    =  F\n")
            f.write("!QUASIPARTICLE (STM)\n")
            f.write("JDos_calc           =  F\n")
            f.write("FindNodes_calc      =  F\n")
            f.write("EffectiveMass_calc  =  F\n")
            f.write("AHC_calc            =  F\n")
            f.write("Boltz_OHE_calc      =  F\n")
            f.write("LOTO_correction     =  F\n")
            f.write("OrbitalTexture_calc    =  F\n")
            f.write("OrbitalTexture_3D_calc =  F\n")
            f.write("LandauLevel_k_calc     =  F\n")
            f.write("LandauLevel_B_calc     =  F\n")
            f.write("LandauLevel_wavefunction_calc     =  F\n")
            f.write("Fit_kp_calc         =  F\n")
            f.write("DMFT_MAG_calc       =  F\n")
            f.write("Translate_to_WS_calc=  F\n")
            f.write("LandauLevel_kplane_calc = F\n")
            f.write("LandauLevel_k_dos_calc = F\n")
            f.write("LandauLevel_B_dos_calc = F \n/\n\n\n")

            f.write("&SYSTEM\n")
            f.write("NSLAB = 20                ! for thin film system\n")
            f.write("NumOccupied = %i        ! NumOccupied\n" %(self.WTwan/2))
            soc = 0 
            if self.SOC:
                soc = 1
            f.write("SOC = %i\n" %(soc))
            f.write("E_FERMI = 0.0\n")
            f.write("surf_onsite= 0.0\n")
            f.write("/\n\n\n")
            f.write("&PARAMETERS \n")
            f.write("Eta_Arc = 0.01     ! infinite small value, like brodening \n")
            f.write("E_arc = 0.0      ! energy level for contour plot of spectrum\n")
            f.write("OmegaNum = 401      ! omega number       \n")
            f.write("OmegaMin = -1.0   ! energy interval\n")
            f.write("OmegaMax =  1.0     ! energy interval\n")
            f.write("Nk1 =  101            ! number k points  odd number would be better\n")
            f.write("Nk2 = 101            ! number k points  odd number would be better\n")
            f.write("Nk3 = 11            ! number k points  odd number would be better\n")
            f.write("NP = 2              ! number of principle layers\n")
            f.write("Gap_threshold = 0.01 ! threshold for FindNodes_calc output\n")
            f.write("/\n\n\n")

            f.write("MILLER_INDEX\n")
            f.write("0 0 1\n\n\n")
            f.write("KPATH_BULK            ! k point path \n")
            f.write("4              ! number of k line only for bulk band\n")
            f.write("G 0.00000 0.00000 0.0000 Z 0.00000 0.00000 0.5000\n")
            f.write("Z 0.00000 0.00000 0.5000 F 0.50000 0.50000 0.0000\n")
            f.write("F 0.50000 0.50000 0.0000 G 0.00000 0.00000 0.0000\n ")
            f.write("G 0.00000 0.00000 0.0000 L 0.50000 0.00000 0.0000 \n\n\n\n")
            ###
            f.write("KPATH_SLAB \n")
            f.write("4        ! numker of k line for 2D case \n")
            f.write("Y 0.0 0.5 G 0.0 0.0\n")
            f.write("G 0.0 0.0 X 0.5 0.0\n")
            f.write("X 0.5 0.0 M 0.5 0.5\n ")
            f.write("M 0.5 0.5 G 0.0 0.0\n\n\n")
            ###
            f.write("KPLANE_SLAB \n")
            f.write("-0.5  -0.5      ! Original point for 2D k plane\n")
            f.write("1.0  0.0      ! The first vector to define 2D k plane \n")
            f.write("0.0  1.0      ! The second vector to define 2D k plane  for arc plots\n\n\n")
            ####
            f.write("KCUBE_BULK\n")
            f.write("0.50  0.50  0.50   ! Original point for 3D k plane \n")
            f.write("1.00  0.00  0.00   ! The first vector to define 3d k space plane\n")
            f.write("0.00  1.00  0.00   ! The second vector to define 3d k space plane\n ")
            f.write("0.00  0.00  1.00   ! The third vector to define 3d k cube\n\n")

            f.write("KPLANE_BULK\n")
            f.write("0.00  0.00  0.00   ! Original point for 3D k plane \n")
            f.write("1.00  0.00  0.00   ! The first vector to define 3d k space plane\n")
            f.write("0.00  1.00  0.00   ! The second vector to define 3d k space plane\n\n\n")
            
        with open(self.SEEDDIR + "WT/wt.slurm",'w',newline= '\n') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --time=5-12:00:00\n")
            f.write("#SBATCH --error=job.%J.err\n")
            f.write("#SBATCH --output=job.%J.out\n")
            f.write("#SBATCH --job-name=%s_WT\n" %(self.seedname))
            f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
            f.write("#SBATCH --mail-type=ALL\n\n\n")
            f.write('#Select how logs get stored\n')
            f.write("mkdir $SLURM_JOB_ID\n")
            f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
            f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
            f.write('#Load Modules\n')
            f.write('ml load wannier_tools\n\n\n')
            f.write("cd $SLURM_SUBMIT_DIR\n")
            f.write('# Create Log File\n')
            f.write('echo $SLURM_SUBMIT_DIR')
            f.write('echo "JobID: $SLURM_JOB_ID" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NODELIST" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NNODES nodes." >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NPROCS processors." >> $debug_logs\n')
            f.write('echo  "Current working directory is `pwd`" >> $debug_logs\n\n\n')
            f.write('# Module debugging\n')
            f.write('module list >> $debug_logs\n')
            f.write('which mpirun >> $debug_logs\n\n\n')
            f.write('#Start Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l: " >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('# run file\n')
            f.write('mpirun -np %i wt.x\n' %(self.NCORE))
            f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
            f.write('#End Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l" >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('#Cleanup\n')
            f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
            f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
            f.write('rm -rf $SLURM_JOB_ID\n')
#def write_slurm(self):
        with open(self.SEEDDIR + "vsp.slurm",'w',newline= '\n') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --time=10:00:00\n")
            f.write("#SBATCH --error=job.%J.err\n")
            f.write("#SBATCH --output=job.%J.out\n")
            f.write("#SBATCH --job-name=%s_NSCF\n" %(self.seedname))
            f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
            f.write("#SBATCH --mail-type=ALL\n\n\n")
            f.write('#Select how logs get stored\n')
            f.write("mkdir $SLURM_JOB_ID\n")
            f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
            f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
            f.write('#Load Modules\n')
            f.write('ml load espresso\n\n\n')
            f.write("cd $SLURM_SUBMIT_DIR\n")
            f.write('# Create Log File\n')
            f.write('echo $SLURM_SUBMIT_DIR')
            f.write('echo "JobID: $SLURM_JOB_ID" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NODELIST" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NNODES nodes." >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NPROCS processors." >> $debug_logs\n')
            f.write('echo  "Current working directory is `pwd`" >> $debug_logs\n\n\n')
            f.write('# Module debugging\n')
            f.write('module list >> $debug_logs\n')
            f.write('which mpirun >> $debug_logs\n\n\n')
            f.write('#Start Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l: " >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('# run file\n')
            f.write('mpirun -n %i vasp_std\n' %(self.NCORE))
            f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
            f.write('#End Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l" >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('#Cleanup\n')
            f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
            f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
            f.write('rm -rf $SLURM_JOB_ID\n')
            f.write('#sbatch p2w.slurm\n')
##############################################
#
#
#
#
#############################################
from ase.io import read, write
from ase.optimize import BFGS
from ase.visualize import view
from ase.io import write as ASEwrite
class CIF2VASP:
    def __init__(self,input):
        crystal = read(input.CIF_file)
        self.CIFDIR = "cifdir/"
        self.key = input.APIKEY # Import API key from input
        self.MaterialID = input.MATERIALID # Set class variable MaterialID from input
        self.seedname = str(input.seedname) # set class variable SEEDNAME from input file
        self.Ecut = input.ECUT #Class variable; in Ryberg, transfer from input file
        self.kpts = input.KMESH ### KPTS from input file for a generator, class variable
        self.SEEDDIR = "OUTPUT/" + 'V-' + self.seedname + '/'  #set the seed directory
        self.EMAIL = input.Email #Class variable; set email for SLURM from input file
        self.NCORE = input.ncore #Class Variable; number of cores, taken from input file
        self.SOC = input.SOC
        self.NBND = input.NUMBANDS
        selectconventional = input.conventional_cell
        selectrelax = input.Relax
        # MAKE DIRECTORIES
        if not os.path.exists(self.SEEDDIR):
            os.makedirs(self.SEEDDIR)
        if not os.path.exists(self.SEEDDIR + "PP"):
            os.makedirs(self.SEEDDIR + "PP")
        if not os.path.exists(self.SEEDDIR + "WT"):
            os.makedirs(self.SEEDDIR + "WT")
        if not os.path.exists(self.SEEDDIR + 'PBE'):
            os.makedirs(self.SEEDDIR + "PBE")
        if not os.path.exists(self.SEEDDIR + "HSE06"):
            os.makedirs(self.SEEDDIR + "HSE06")
        if not os.path.exists(self.SEEDDIR + "W90"):
            os.makedirs(self.SEEDDIR + "W90")
        ################################
        # Import structure from CIF file
        # ##############################
        self.cell = crystal.get_cell()
        self.cord = crystal.get_scaled_positions(wrap=True)
        self.species = crystal.get_chemical_symbols()
        self.UNIQUE_ATOMS = np.unique(self.species)
        self.atomic_numbers = crystal.get_atomic_numbers()
        self.ntypat = np.unique(self.atomic_numbers).size
        self.natom = np.array(self.atomic_numbers).size
        self.atomnum = np.unique(self.atomic_numbers)[::1]
        lattice_old = pd.DataFrame(np.array(self.cell))
        self.lattice = lattice_old.to_string(header=False, index=False)
        ##############
        angles = self.cell.angles()
        abc = np.array(self.cell.lengths())
        self.a = abc[0]
        self.b = abc[1]
        self.c = abc[2]
        self.alpha = angles[0]
        self.beta = angles[1]
        self.gamma = angles[2]
        self.typat = ""
        for x in self.atomic_numbers:
            for at in self.atomnum:
                if at == x:
                        self.typat += ''.join(map(str,np.where(self.atomnum == at)
                            [0] + 1)) + " "
        # SOC CASE
                # SOC CASE
        self.SOC = input.SOC
        self.noncolin = '.FALSE.'
        self.lspinorb = '.FALSE.'
        if input.SOC:
            self.noncolin = '.TRUE.'
            self.lspinorb = '.TRUE.'
        num_bands = input.NUMBANDS
        num_wann = num_bands - (num_bands%2) # FORCE EVEN
        self.wan = num_wann
        self.WTwan = num_wann
        if self.SOC:
            self.WTwan = num_wann/2
        num_bands = 0
        #for x in struct.species:
           # PP = gb.glob(self.pseudodir + str(x) + ".*") 
            #readpp = minidom.parse(''.join(map(str,PP)))
            #items = readpp.getElementsByTagName('atom')
            #num_bands += float(items[0].attributes['valence'].value)
        num_wann = num_bands - (num_bands%2) # FORCE EVEN

    def generate(self):
        with open(self.SEEDDIR + 'pbe/' + "INCAR",'w',newline='\n') as f:
            f.write("ALGO = Fast\n")
            f.write("EDIFF = 1.0e-8\n")
            f.write("ENCUT = 500\n")
            f.write("ISIF = 2\n")
            f.write("ISMEAR = -5\n")
            f.write("KPAR = 2\n")
            f.write("LCHARG = .FALSE.\n")
            f.write("LREAL = .FALSE.\n")
            f.write("NBANDS = %i\n" %(self.NBND))
            f.write("NEDOS = 2000\n")
            f.write("NPAR = 1\n")
            f.write("PREC = Accurate\n")
            f.write("SYSTEM = %s\n" %(self.seedname))
        with open(self.SEEDDIR + 'pbe/' + "KPOINTS",'w',newline='\n') as f:
            f.write("Automatic %ix%ix%i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("0\n")
            f.write("Gamma\n")
            f.write("\t%i %i %i"%(self.kpts[0],self.kpts[1],self.kpts[2]))
        with open(self.SEEDDIR + 'pbe/' + "POSCAR",'w',newline='\n') as f:
            f.write("%s\n" %(self.seedname))
            f.write("1\n")
            f.write(str(self.lattice) + "\n")
            f.write("\t" + ' '.join(map(str,self.UNIQUE_ATOMS)) + "\n")
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.atomic_numbers == self.atomnum[i])).size
            f.write("\t" + ' '.join(map(str,typnum)) + "\n")
            f.write("Direct\n")
            for index,x in enumerate(self.cord):
                f.write("%f %f %f\n" %(x[0], x[1], x[2]))
    #def write_hse(self):
        with open(self.SEEDDIR + 'HSE06/' + "INCAR",'w',newline='\n') as f:
            f.write("ALGO = Fast\n")
            f.write("EDIFF = 1.0e-8\n")
            f.write("ENCUT = 500\n")
            f.write("ISIF = 2\n")
            f.write("ISMEAR = -5\n")
            f.write("KPAR = 2\n")
            f.write("LCHARG = .FALSE.\n")
            f.write("LREAL = .FALSE.\n")
            f.write("NBANDS = %i\n" %(self.NBND))
            f.write("NEDOS = 2000\n")
            f.write("NPAR = 1\n")
            f.write("PREC = Accurate\n")
            f.write("SYSTEM = %s\n" %(self.seedname))
            f.write("LHFCALC = .TRUE.\n")
            f.write("HFSCREEN = 0.2\n")
        with open(self.SEEDDIR + 'HSE06/' + "KPOINTS",'w',newline='\n') as f:
            f.write("Automatic %ix%ix%i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("0\n")
            f.write("Gamma\n")
            f.write("\t%i %i %i"%(self.kpts[0],self.kpts[1],self.kpts[2]))
        with open(self.SEEDDIR + 'HSE06/' + "POSCAR",'w',newline='\n') as f:
            f.write("%s\n" %(self.seedname))
            f.write("1\n")
            f.write(str(self.lattice) + "\n")
            f.write("\t" + ' '.join(map(str,self.UNIQUE_ATOMS)) + "\n")
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.atomic_numbers == self.atomnum[i])).size
            f.write("\t" + ' '.join(map(str,typnum)) + "\n")
            f.write("Direct\n")
            for index,x in enumerate(self.cord):
                f.write("%f %f %f\n" %(x[0], x[1], x[2]))
    #def write_w90(self):
        with open(self.SEEDDIR + 'W90/' + "INCAR",'w',newline='\n') as f:
            f.write("ALGO = Fast\n")
            f.write("EDIFF = 1.0e-8\n")
            f.write("ENCUT = 500\n")
            f.write("ISIF = 2\n")
            f.write("ISMEAR = -5\n")
            f.write("KPAR = 2\n")
            f.write("LCHARG = .FALSE.\n")
            f.write("LREAL = .FALSE.\n")
            f.write("NBANDS = %i\n" %(self.NBND))
            f.write("NEDOS = 2000\n")
            f.write("NPAR = 1\n")
            f.write("PREC = Accurate\n")
            f.write("SYSTEM = %s\n" %(self.seedname))
            f.write("LHFCALC = .TRUE.\n")
            f.write("HFSCREEN = 0.2\n")
            f.write("LWANNIER90 = .TRUE.\n")
        with open(self.SEEDDIR + 'W90/' + "KPOINTS",'w',newline='\n') as f:
            f.write("Automatic %ix%ix%i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("0\n")
            f.write("Gamma\n")
            f.write("\t%i %i %i"%(self.kpts[0],self.kpts[1],self.kpts[2]))
        with open(self.SEEDDIR + 'W90/' + "POSCAR",'w',newline='\n') as f:
            f.write("%s\n" %(self.seedname))
            f.write("1\n")
            f.write(str(self.lattice) + "\n")
            f.write("\t" + ' '.join(map(str,self.UNIQUE_ATOMS)) + "\n")
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.atomic_numbers == self.atomnum[i])).size
            f.write("\t" + ' '.join(map(str,typnum)) + "\n")
            f.write("Direct\n")
            for index,x in enumerate(self.cord):
                f.write("%f %f %f\n" %(x[0], x[1], x[2]))
        with open(self.SEEDDIR + 'W90/wannier90.win','w',newline='\n') as f:
            f.write("!write_hr = .TRUE.\n")
            f.write("!write_xyz = .TRUE.\n")
            f.write("guiding_centres= .TRUE.\n")
            f.write("!wannier_plot = .TRUE. \n")
            f.write("spinors = %s\n" %(self.lspinorb))
            f.write("num_wann = %i\n" %(self.wan))
            f.write("dis_num_iter=1000\n")
            f.write("num_iter = 2000\n\n\n\n")
            f.write("begin unit_cell_cart\n")
            f.write(str(self.lattice) + "\n") ### atomic structure
            f.write("end unit_cell_cart\n\n\n")
            f.write("begin atoms_frac\n")
            for index,x in enumerate(self.cord):
                f.write("%s %f %f %f\n" %(str(self.species[index]),x[0], x[1], x[2]))
            f.write("end atoms_frac\n\n\n")
            f.write("begin projections \n")
            f.write("random \n")
            f.write("end projections\n\n\n")
            f.write(kmesh.WANNIER(self.kpts[0],self.kpts[1],self.kpts[2]))

    #def write_wt(self):
        with open(self.SEEDDIR + "WT/wt.in",'w',newline='\n') as f:
            f.write("#### wt file ###########\n")
            f.write("&TB_FILE\n")
            f.write("Hrfile = '%s_hr.dat'\n" %(self.seedname))
            f.write("Package = 'ESPRESSO'\n")
            f.write("/\n\n\n")
            f.write("LATTICE\n")
            f.write("Angstrom\n")
            f.write(str(self.lattice) + "\n\n\n") ### atomic structure
            f.write("ATOM_POSITIONS\n")
            f.write("%i !Number of atoms for projectors\n" %(len(self.cord)))
            f.write("Direct ! Direct or Cartisen coordinateS\n")
            for index,x in enumerate(self.cord):
                f.write("%s %f %f %f\n" %(str(self.species[index]),x[0], x[1], x[2]))
            f.write("\n\n\n\n")
            f.write("PROJECTORS\n")
            #### Places prjectors and stoms into PRJCARD in wannier tools
            proj = int(self.WTwan/self.natom)
            proj_rem = self.WTwan % self.natom
            sumproj = 0
            PROJPARAM = ""
            projector = ["s", "pz","px", "py"]
            LATPARAM = str(self.lattice) + "\n" ### atomic structur
            PROJCARD =""
            ATOMCARD = ""
            appen = 0 #### append number
            for index,atom in enumerate(self.species):
                if proj_rem > 0:
                    appen = 1
                proj_rem = proj_rem - 1
                PROJPARAM += str(proj + appen) + " "
                PROJCARD += str(atom) + " " + ' '.join(map(str,projector)) + "\n"
                #print(str(atom) + str(proj + appen))
                sumproj += proj + appen
                #print(sumproj)
                appen = 0
            f.write(PROJPARAM + "\n")
            f.write(PROJCARD)
            f.write("\n\n\n")

            f.write("&CONTROL\n")
            f.write("! BULK BAND CALCULATIONS \n")
            f.write("BulkBand_calc       =  F\n")
            f.write("BulkBand_plane_calc =  F \n")
            f.write("BulkFS_calc         =  F \n")
            f.write("BulkFS_Plane_calc   =  F\n")
            f.write("SlabBand_calc       =  F\n")
            f.write("Dos_calc            =  F\n")
            f.write("! BULK GAP\n")
            f.write("BulkGap_cube_calc   =  F\n")
            f.write("BulkGap_plane_calc  =  F\n")
            f.write("! SURFACE STATES\n")
            f.write("SlabSS_calc         =  T\n")
            f.write("SlabArc_calc        =  T\n")
            f.write("SlabSpintexture_calc =  F\n")
            f.write("! TOPO INV\n")
            f.write("wanniercenter_calc   = F\n")
            f.write("BerryPhase_calc     =  F\n")
            f.write("BerryCurvature_calc =  F\n")
            f.write("BerryCurvature_slab_calc =  F\n")
            f.write("Z2_3D_calc          =  F\n")
            f.write("WeylChirality_calc  =  F\n")
            f.write("NLChirality_calc    =  F\n")
            f.write("Chern_3D_calc       =  F\n")
            f.write("MirrorChern_calc    =  F\n")
            f.write("!QUASIPARTICLE (STM)\n")
            f.write("JDos_calc           =  F\n")
            f.write("FindNodes_calc      =  F\n")
            f.write("EffectiveMass_calc  =  F\n")
            f.write("AHC_calc            =  F\n")
            f.write("Boltz_OHE_calc      =  F\n")
            f.write("LOTO_correction     =  F\n")
            f.write("OrbitalTexture_calc    =  F\n")
            f.write("OrbitalTexture_3D_calc =  F\n")
            f.write("LandauLevel_k_calc     =  F\n")
            f.write("LandauLevel_B_calc     =  F\n")
            f.write("LandauLevel_wavefunction_calc     =  F\n")
            f.write("Fit_kp_calc         =  F\n")
            f.write("DMFT_MAG_calc       =  F\n")
            f.write("Translate_to_WS_calc=  F\n")
            f.write("LandauLevel_kplane_calc = F\n")
            f.write("LandauLevel_k_dos_calc = F\n")
            f.write("LandauLevel_B_dos_calc = F \n/\n\n\n")

            f.write("&SYSTEM\n")
            f.write("NSLAB = 20                ! for thin film system\n")
            f.write("NumOccupied = %i        ! NumOccupied\n" %(self.WTwan/2))
            soc = 0 
            if self.SOC:
                soc = 1
            f.write("SOC = %i\n" %(soc))
            f.write("E_FERMI = 0.0\n")
            f.write("surf_onsite= 0.0\n")
            f.write("/\n\n\n")
            f.write("&PARAMETERS \n")
            f.write("Eta_Arc = 0.01     ! infinite small value, like brodening \n")
            f.write("E_arc = 0.0      ! energy level for contour plot of spectrum\n")
            f.write("OmegaNum = 401      ! omega number       \n")
            f.write("OmegaMin = -1.0   ! energy interval\n")
            f.write("OmegaMax =  1.0     ! energy interval\n")
            f.write("Nk1 =  101            ! number k points  odd number would be better\n")
            f.write("Nk2 = 101            ! number k points  odd number would be better\n")
            f.write("Nk3 = 11            ! number k points  odd number would be better\n")
            f.write("NP = 2              ! number of principle layers\n")
            f.write("Gap_threshold = 0.01 ! threshold for FindNodes_calc output\n")
            f.write("/\n\n\n")

            f.write("MILLER_INDEX\n")
            f.write("0 0 1\n\n\n")
            f.write("KPATH_BULK            ! k point path \n")
            f.write("4              ! number of k line only for bulk band\n")
            f.write("G 0.00000 0.00000 0.0000 Z 0.00000 0.00000 0.5000\n")
            f.write("Z 0.00000 0.00000 0.5000 F 0.50000 0.50000 0.0000\n")
            f.write("F 0.50000 0.50000 0.0000 G 0.00000 0.00000 0.0000\n ")
            f.write("G 0.00000 0.00000 0.0000 L 0.50000 0.00000 0.0000 \n\n\n\n")
            ###
            f.write("KPATH_SLAB \n")
            f.write("4        ! numker of k line for 2D case \n")
            f.write("Y 0.0 0.5 G 0.0 0.0\n")
            f.write("G 0.0 0.0 X 0.5 0.0\n")
            f.write("X 0.5 0.0 M 0.5 0.5\n ")
            f.write("M 0.5 0.5 G 0.0 0.0\n\n\n")
            ###
            f.write("KPLANE_SLAB \n")
            f.write("-0.5  -0.5      ! Original point for 2D k plane\n")
            f.write("1.0  0.0      ! The first vector to define 2D k plane \n")
            f.write("0.0  1.0      ! The second vector to define 2D k plane  for arc plots\n\n\n")
            ####
            f.write("KCUBE_BULK\n")
            f.write("0.50  0.50  0.50   ! Original point for 3D k plane \n")
            f.write("1.00  0.00  0.00   ! The first vector to define 3d k space plane\n")
            f.write("0.00  1.00  0.00   ! The second vector to define 3d k space plane\n ")
            f.write("0.00  0.00  1.00   ! The third vector to define 3d k cube\n\n")

            f.write("KPLANE_BULK\n")
            f.write("0.00  0.00  0.00   ! Original point for 3D k plane \n")
            f.write("1.00  0.00  0.00   ! The first vector to define 3d k space plane\n")
            f.write("0.00  1.00  0.00   ! The second vector to define 3d k space plane\n\n\n")
            
        with open(self.SEEDDIR + "WT/wt.slurm",'w',newline= '\n') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --time=5-12:00:00\n")
            f.write("#SBATCH --error=job.%J.err\n")
            f.write("#SBATCH --output=job.%J.out\n")
            f.write("#SBATCH --job-name=%s_WT\n" %(self.seedname))
            f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
            f.write("#SBATCH --mail-type=ALL\n\n\n")
            f.write('#Select how logs get stored\n')
            f.write("mkdir $SLURM_JOB_ID\n")
            f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
            f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
            f.write('#Load Modules\n')
            f.write('ml load wannier_tools\n\n\n')
            f.write("cd $SLURM_SUBMIT_DIR\n")
            f.write('# Create Log File\n')
            f.write('echo $SLURM_SUBMIT_DIR')
            f.write('echo "JobID: $SLURM_JOB_ID" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NODELIST" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NNODES nodes." >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NPROCS processors." >> $debug_logs\n')
            f.write('echo  "Current working directory is `pwd`" >> $debug_logs\n\n\n')
            f.write('# Module debugging\n')
            f.write('module list >> $debug_logs\n')
            f.write('which mpirun >> $debug_logs\n\n\n')
            f.write('#Start Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l: " >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('# run file\n')
            f.write('mpirun -np %i wt.x\n' %(self.NCORE))
            f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
            f.write('#End Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l" >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('#Cleanup\n')
            f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
            f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
            f.write('rm -rf $SLURM_JOB_ID\n')
#def write_slurm(self):
        with open(self.SEEDDIR + "vsp.slurm",'w',newline= '\n') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --time=10:00:00\n")
            f.write("#SBATCH --error=job.%J.err\n")
            f.write("#SBATCH --output=job.%J.out\n")
            f.write("#SBATCH --job-name=%s_NSCF\n" %(self.seedname))
            f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
            f.write("#SBATCH --mail-type=ALL\n\n\n")
            f.write('#Select how logs get stored\n')
            f.write("mkdir $SLURM_JOB_ID\n")
            f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
            f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
            f.write('#Load Modules\n')
            f.write('ml load espresso\n\n\n')
            f.write("cd $SLURM_SUBMIT_DIR\n")
            f.write('# Create Log File\n')
            f.write('echo $SLURM_SUBMIT_DIR')
            f.write('echo "JobID: $SLURM_JOB_ID" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NODELIST" >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NNODES nodes." >> $debug_logs\n')
            f.write('echo "Running on $SLURM_JOB_NPROCS processors." >> $debug_logs\n')
            f.write('echo  "Current working directory is `pwd`" >> $debug_logs\n\n\n')
            f.write('# Module debugging\n')
            f.write('module list >> $debug_logs\n')
            f.write('which mpirun >> $debug_logs\n\n\n')
            f.write('#Start Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l: " >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('# run file\n')
            f.write('mpirun -n %i vasp_std\n' %(self.NCORE))
            f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
            f.write('#End Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l" >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('#Cleanup\n')
            f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
            f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
            f.write('rm -rf $SLURM_JOB_ID\n')
            f.write('#sbatch p2w.slurm\n')