### VASP ###

import os
#### PyMATGEN
import pymatgen as mg
from pymatgen import Composition
from pymatgen import Lattice, Structure, Molecule, Specie
import pymatgen.core.periodic_table as Element
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
import input
import kmesh

class MAT2VASP:
    def __init__(self):
        self.CIFDIR = "cifdir/"
        self.key = input.APIKEY # Import API key from input
        self.MaterialID = input.MATERIALID # Set class variable MaterialID from input
        self.seedname = str(input.SEEDNAME) # set class variable SEEDNAME from input file
        self.Ecut = input.ECUT #Class variable; in Ryberg, transfer from input file
        self.kpts = input.KMESH ### KPTS from input file for a generator, class variable
        self.SEEDDIR = "OUTPUT/" + 'VSP-' + self.seedname + '/'  #set the seed directory
        self.pseudodir = "ABINIT/ABINITPP/ATOMICDATA/" # Set the directory for pseudopotentials
        self.EMAIL = input.EMAIL #Class variable; set email for SLURM from input file
        self.NCORE = input.ncore #Class Variable; number of cores, taken from input file
        self.SOC = input.SOC
        self.NBND = input.NUMBANDS
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
            struct = m.get_structure_by_material_id(self.MaterialID,final=True,conventional_unit_cell=False) # GET STRUCT FROM MATERIALS PROJECT
            self.struct = m.get_structure_by_material_id(self.MaterialID,final=True,conventional_unit_cell=False) # GET STRUCT FROM MATERIALS PROJECT
        C = Composition(str(struct.formula))
        ########## MAKE A CIF FILE ###
        cif = CifWriter(struct) # Start generator
        cif.write_file(self.CIFDIR + self.seedname + ".cif") #write cif file to directory
        ciffile = self.CIFDIR + self.seedname + ".cif" # define ciffile
        cif = open(ciffile) # open cif file
        parser = CifParser(ciffile) #Start Parser
        ####### PARSE VARIABLES FROM CIF ##########
        chk_a = "_cell_length_a"
        chk_b = "_cell_length_b"
        chk_c = "_cell_length_c"
        chk_alpha = "_cell_angle_alpha"
        chk_beta = "_cell_angle_beta"
        chk_gamma = "_cell_angle_gamma"
        chk_name = "_chemical_formula_structural"
        for x in cif:
            if chk_a in x:
                self.a = x.replace(chk_a, '').replace('\n', '')
            if chk_b in x:
                self.b = x.replace(chk_b, '').replace('\n', '')
            if chk_c in x:
                self.c = x.replace(chk_c, '').replace('\n', '')
            if chk_alpha in x:
                self.alpha = x.replace(chk_alpha, '').replace('\n', '')
            if chk_beta in x:
                self.beta = x.replace(chk_beta, '').replace('\n', '')
            if chk_gamma in x:
                self.gamma = x.replace(chk_gamma, '').replace('\n', '')
            #if chk_name in x:
            #    seedname = x.replace(chk_name, '')
            #    seedname = seedname.replace('\n', '').translate(
            #    {ord(i): None for i in ' '})
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
        if self.SOC:
            self.WTwan = self.wan/2
        num_bands = 0
        for x in struct.species:
            PP = gb.glob(self.pseudodir + str(x) + ".*") 
            readpp = minidom.parse(''.join(map(str,PP)))
            items = readpp.getElementsByTagName('atom')
            num_bands += float(items[0].attributes['valence'].value)
        num_wann = num_bands - (num_bands%2) # FORCE EVEN
    def PMG2WAN(self):
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



#run = MAT2VASP()
