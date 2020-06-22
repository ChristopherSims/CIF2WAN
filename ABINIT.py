#ABI
import os
#### PyMATGEN
import math
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
class MAT2ABINIT:
    def __init__(self):
        self.CIFDIR = "cifdir/"
        self.key = input.APIKEY # Import API key from input
        self.MaterialID = input.MATERIALID # Set class variable MaterialID from input
        self.seedname = str(input.SEEDNAME) # set class variable SEEDNAME from input file
        self.Ecut = input.ECUT #Class variable; in Ryberg, transfer from input file
        self.kpts = input.KMESH ### KPTS from input file for a generator, class variable
        self.SEEDDIR = "OUTPUT/" + 'Abi-' + self.seedname + '/'  #set the seed directory
        self.pseudodir = "ABINIT/ABINITPP/ATOMICDATA/" # Set the directory for pseudopotentials
        self.EMAIL = input.EMAIL #Class variable; set email for SLURM from input file
        self.NCORE = input.ncore #Class Variable; number of cores, taken from input file
        self.SOC = input.SOC
        
        # MAKE DIRECTORIES
        if not os.path.exists(self.SEEDDIR):
            os.makedirs(self.SEEDDIR)
        if not os.path.exists(self.SEEDDIR + "PP"):
            os.makedirs(self.SEEDDIR + "PP")
        if not os.path.exists(self.SEEDDIR + "WT"):
            os.makedirs(self.SEEDDIR + "WT")
        #import structure from PYMATGEN
        with matproj.MPRester(self.key) as m:
            struct = m.get_structure_by_material_id(self.MaterialID,final=False,conventional_unit_cell=False) # GET STRUCT FROM MATERIALS PROJECT
            self.struct = m.get_structure_by_material_id(self.MaterialID,final=False,conventional_unit_cell=False) # GET STRUCT FROM MATERIALS PROJECT
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
        #print(float(self.a)*1.88973)
        #print((float(self.b)/float(self.a))*1.88973)
        #print((float(self.c)/float(self.a))*1.88973)
        #print(math.cos(float(self.alpha)))
        #print(math.cos(float(self.beta)))
        #print(math.cos(float(self.gamma)))
        

        for x in struct.atomic_numbers:
            for at in self.atomnum:
                if at == x:
                        self.typat += ''.join(map(str,np.where(self.atomnum == at)
                            [0] + 1)) + " "
        # SOC CASE
        num_bands = 0
        for x in struct.species:
            PP = gb.glob(self.pseudodir + str(x) + ".*") 
            readpp = minidom.parse(''.join(map(str,PP)))
            items = readpp.getElementsByTagName('atom')
            num_bands += float(items[0].attributes['valence'].value)
        num_wann = num_bands - (num_bands%2) # FORCE EVEN
        self.wan = num_wann
        if self.SOC:
            self.WTwan = self.wan/2

    def PMG2WAN(self):
        with open(self.SEEDDIR + self.seedname + ".in",'w',newline='\n') as file:
            file.write("ndtset %i \n" % (2)) # for wanniertools
            file.write("enunit 1 #change to EV\n")
            file.write("prtvol 2 #print all kpoints\n")
            file.write("##### " + self.seedname + " stucture ######\n") # header
            file.write("acell" + self.a + self.b + self.c + " angstrom \n") # cell parameters
            file.write("angdeg" + self.alpha + self.beta + self.gamma + "\n") #angle param
            file.write("natom  " + str(self.natom) + "\n") #number of atoms
            #### write atom positions ####
            file.write("xred\n") # Fractional coordinates
            for x in self.cord:
                file.write("%12.8f%12.8f%12.8f\n" %(x[0],x[1],x[2]))
            ##### end write atom positions ####
            file.write("ntypat " + repr(self.ntypat) + " #number of unique atoms\n")

            file.write("typat " + self.typat + "#identity of atoms\n")
            #### Write Z of atoms ####
            file.write("znucl ")
            for item in np.unique(self.struct.atomic_numbers):
                file.write(str(item) + " ")
            file.write(" # atomic mass of atoms\n")
            #### END of Z write ###
            file.write("optforces 0 \n") # enforce calculation of forces at each SCF step
            ############ Plane wave  ########
            file.write("########### PLANE WAVE INFO ##########\n")
            file.write("ecut 20  ry \n")
            file.write("pawecutdg 22 ry\n")
            file.write("pawovlp  15\n")
            file.write("nstep 100\n")
            #### USE DFT+U ####
            is_HF = False
            csvfile = open('HF.csv', newline='')
            HF = csv.reader(csvfile, delimiter=',')
            DFTUmat = np.full(len(self.UNIQUE_ATOMS),-1)
            Umat = np.zeros(len(self.UNIQUE_ATOMS))
            Jmat = np.zeros(len(self.UNIQUE_ATOMS))
            for index,atom in enumerate(self.UNIQUE_ATOMS):
                for el in HF:
                    if str(atom) in str(el[0]) and len(str(atom)) == len(el[0]):
                        is_HF = True
                        DFTUmat[index] = 2
                        Umat[index] = el[1] ### U TERM
                        Jmat[index] = el[2] ### J0 terM
                ### reset csv reader
                csvfile = open('HF.csv', newline='')
                HF = csv.reader(csvfile, delimiter=',')
            if is_HF:
               file.write("usepawu 1\n")
               file.write("lpawu %s\n" %(' '.join(map(str, DFTUmat))))
               file.write("upawu %s\n" %(' '.join(map(str, Umat))))
               file.write("jpawu %s\n" %(' '.join(map(str, Jmat))))
            #### USE SOC ###
            if self.SOC:
                file.write("pawspnorb 1\n")
                file.write("nspinor 2\n")
                file.write("nsppol  1\n")
                file.write("nspden  1\n")
            file.write("######### END plane wave info ############\n\n\n")
            ######### END plane wave info ############
            #############SCF STEP ##########
            #
            #
            #
            #
            ##################################
            file.write("########### SCF STEP ##########\n")
            file.write("iscf 17\n")
            file.write("tolvrs1  1.00d-12 \n")
            file.write("ngkpt1 8 8 8\n")
            if self.SOC:
                file.write("kptopt1 4\n")
            else:
                file.write("kptopt1 1\n")
            file.write("nshiftk 1 #just one shift is supported by wannier90\n")
            file.write("shiftk 0.00   0.00   0.00  #no shift \n")
            file.write("prtden1 1\n")
            file.write("istwfk1 512*1 #Controls the form of the wavefunctions\n")
            ############ End SCF Step ################
            ############# NSCF STEP  #################
            #
            #
            #
            #############################################
            file.write("########### NSCF STEP ##########\n")
            file.write("prtvol2  1\n")
            file.write("prtden2 1\n")  
            file.write("iscf2 -2 \n")
            file.write("nstep2 0\n")
            file.write("tolwfr2 1.d-12\n")
            file.write("getwfk2 -1 # Get from SCF\n")
            file.write("getden2 -1 # Get den from SCF\n")
            file.write("istwfk2  %i*1  #Controls the form of the wavefunctions\n" %(self.kpts[0]*self.kpts[1]*self.kpts[2]))
            file.write("prtwant2 2   # Call to Wannier90\n")
            file.write("w90iniprj2 2\n")
            file.write("w90prtunk2 0   #Prints UNK files (for plotting the Wannier functions) \n")
            file.write("kptopt2 0\n")
            file.write(kmesh.ABINIT(self.kpts[0],self.kpts[1],self.kpts[2],2))
    ### WRITE the SEEDNAME_FILES file ######
    #def write_files(self):
        with open(self.SEEDDIR + self.seedname + ".files",'w',newline='\n') as file:
            file.write(self.seedname + ".in" + "\n")
            file.write(self.seedname + ".out" + "\n")
            file.write(self.seedname + "_i" + "\n")
            file.write(self.seedname + "_o" + "\n")
            file.write(self.seedname + "_1" + "\n")
                    #### WRITE psuedopotentails to PPDIR
            PPSAVE = []
            znulc = ""
            for x in np.unique(self.struct.species)[::-1]:
                PP = gb.glob(self.pseudodir + str(x) + ".*")
                sendPPdir = self.SEEDDIR + "PP/" + ''.join(map(str,PP))[len(self.pseudodir):] # get directory to send PP
                file.write("PP/" + ''.join(map(str,PP))[len(self.pseudodir):] + " \n")
                copyfile(''.join(map(str,PP)), ''.join(map(str,sendPPdir)))
                #### GET Z for znucl
                PP = gb.glob(self.pseudodir + str(x) + ".*") 
                readpp = minidom.parse(''.join(map(str,PP)))
                items = readpp.getElementsByTagName('atom')
                znulc += str(items[0].attributes['Z'].value) + " "
        return(file)
    #def write_w90(self):
        #########WANNIER90 FILE####################33
        with open(self.SEEDDIR + "w90.win",'w',newline='\n') as f:
            f.write("write_hr = .TRUE.\n")
            f.write("write_xyz = .TRUE.\n")
            f.write("wannier_plot = .TRUE. \n")
            f.write("spinors = %s\n" %(self.SOC))
            f.write("num_wann = %i\n" %(self.wan))
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
            f.write("Hrfile = 'w90_hr.dat'\n")
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
       

    #def write_slurm(self):
        with open(self.SEEDDIR + "abi.slurm",'w',newline= '\n') as f:
            f.write("#!/bin/bash\n")
            #f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --ntasks=1\n" )
            f.write("#SBATCH --time=5-12:00:00\n")
            f.write("#SBATCH --error=job.%J.err\n")
            f.write("#SBATCH --output=job.%J.out\n")
            f.write("#SBATCH --job-name=%s_ABI\n" %(self.seedname))
            f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
            f.write("#SBATCH --mail-type=ALL\n\n\n")
            f.write('#Select how logs get stored\n')
            f.write("mkdir $SLURM_JOB_ID\n")
            f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
            f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
            f.write('#Load Modules\n')
            f.write('ml load abinit\n\n\n')
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
            f.write('abinit <%s.files > log 2> err\n' %(self.seedname))
            f.write('sleep 3\n')
            f.write('grep CONV w90.wout >> Wconv.txt\n')
            f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
            f.write('#End Timestamp\n')
            f.write('date >> $benchmark_logs\n')
            f.write('echo "ulimit -l" >> $benchmark_logs\n')
            f.write('ulimit -l >> $benchmark_logs\n\n\n')
            f.write('#Cleanup\n')
            f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
            f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
            f.write('rm -rf $SLURM_JOB_ID\n')
        ######### WRITE WT SLURM #########
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

    #def WRITE_CLEAN(self):
        with open(self.SEEDDIR + "/cleandft.sh", 'w',newline='\n') as file:
            file.write("rm job.*\n")
            file.write("rm *.1o*\n")
            file.write("rm *.out*\n")
            file.write("rm log\n")
            file.write("rm w90.eig\n")
            file.write("rm w90.mmn\n")
            file.write("rm w90.werr\n")
            file.write("rm err\n")
            file.write("rm w90.nnkp\n")
            file.write("rm w90.wout\n")
            file.write("rm wannier90random.amn\n")
#run = MAT2ABINIT()
# run.write_in()
# run.write_files()
# run.write_w90()
# run.write_wt()
# run.write_slurm()