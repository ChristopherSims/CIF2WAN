#ABI
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
class MAT2ESPRESSO:
    def __init__(self):
        self.key = input.APIKEY # Import API key from input
        self.MaterialID = input.MATERIALID # Set class variable MaterialID from input
        self.seedname = str(input.SEEDNAME) # set class variable SEEDNAME from input file
        self.Ecut = input.ECUT #Class variable; in Ryberg, transfer from input file
        self.kpts = input.KMESH ### KPTS from input file for a generator, class variable
        self.SEEDDIR = "OUTPUT/" + 'QE-' + self.seedname + '/'  #set the seed directory
        self.pseudodir = "Pseudopotentials/Espresso/" # Set the directory for pseudopotentials
        self.EMAIL = input.EMAIL #Class variable; set email for SLURM from input file
        self.NCORE = input.ncore #Class Variable; number of cores, taken from input file
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
        #import structure from PYMATGEN
        with matproj.MPRester(self.key) as m:
            struct = m.get_structure_by_material_id(self.MaterialID,final=selectrelax,conventional_unit_cell=selectconventional) # GET STRUCT FROM MATERIALS PROJECT
            self.struct = m.get_structure_by_material_id(self.MaterialID,final=selectrelax,conventional_unit_cell=selectconventional) # GET STRUCT FROM MATERIALS PROJECT
        C = Composition(str(struct.formula))
        self.cord = struct.frac_coords
        self.ntypat = np.unique(struct.atomic_numbers).size
        self.natom = np.array(struct.atomic_numbers).size
        self.UNIQUE_ATOMS = np.unique(struct.species)
        self.atomnum = np.unique(struct.atomic_numbers)[::1]
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
        else:
             self.WTwan = self.wan
    def PMG2WAN(self):
        with open(self.SEEDDIR + self.seedname + '.scf.in','w',newline='\n') as scf:
            scf.write("&CONTROL\n")
            scf.write("calculation ='scf'\n")
            scf.write("prefix = '" + self.seedname +"'\n")
            scf.write("outdir = './bin'\n")
            scf.write("pseudo_dir = './PP/'\n")
            scf.write("verbosity='high'\n")
            scf.write("/\n\n")
            ################# ATOM PARAM ###############
            scf.write("&system\n")
            scf.write("ibrav = 0\n")
            scf.write("nat = %i\n" %(self.natom))
            scf.write("ntyp = %i\n" %(self.ntypat))
            scf.write("ecutwfc =%i  !Ryberg\n" %(self.Ecut))
            scf.write("ecutrho =%i\n" %(self.Ecut*4))
            scf.write("occupations = 'smearing'\n")
            scf.write("smearing = 'gaussian'\n")
            scf.write("degauss = 0.01\n\n")
            scf.write("!SOC\n")
            scf.write("noncolin = %s\n" %(self.noncolin))
            scf.write("lspinorb = %s\n" %(self.lspinorb))
            for index,atom in enumerate(self.UNIQUE_ATOMS):
                scf.write("starting_magnetization(%i) = %i\n" %(index+1,0))
            scf.write("\n")
            ######## HUBBARD for HF ###########
            is_HF = False
            csvfile = open("HF.csv", newline='')
            HF = csv.reader(csvfile, delimiter=',')
            for index,atom in enumerate(self.UNIQUE_ATOMS):
                for el in HF:
                    if str(atom) in str(el[0]) and len(str(atom)) == len(el[0]):
                        is_HF = True
                        scf.write("Hubbard_U(%i) = %s\n" %(index+1, el[1])) ### U TERM
                        scf.write("Hubbard_J0(%i) = %s\n" %(index+1, el[2])) ### J0 terM
                ### reset csv reader
                csvfile = open('HF.csv', newline='')
                HF = csv.reader(csvfile, delimiter=',')
            if is_HF:
                scf.write("lda_plus_u = .TRUE.\n")
                scf.write("lda_plus_u_kind = 1\n")
            scf.write("\n")
            scf.write("/\n\n")
            #### ELECTRON CARD #####
            scf.write("&ELECTRONS\n")
            scf.write("conv_thr = 1.0d-7\n")
            scf.write("mixing_beta = 0.495\n")
            scf.write("mixing_mode = 'TF'\n")
            scf.write("diagonalization= 'david'\n")
            scf.write("adaptive_thr=.true.\n")
            scf.write("/\n\n")
            scf.write("ATOMIC_SPECIES\n")
            #### ATOMIC SPECIES CARD
            for index,atom in enumerate(self.UNIQUE_ATOMS):
                scf.write("%s " %(str(atom))) #Write element symbol to file
                PP = gb.glob(self.pseudodir + str(atom) + "*") # get pseudopotential
                PP = PP[0] # get first Psuedopotential found, for robustness
                element = mg.Element(str(atom)) # Get element information
                amass = float(element.atomic_mass) # Get atomic mass
                scf.write("%s " %(amass))    #write atomic mass
                sendPPdir = self.SEEDDIR + "PP/" + ''.join(map(str,PP))[len(self.pseudodir):] # get directory to send PP
                scf.write(''.join(map(str,PP))[len(self.pseudodir):] + " \n") # write PP to QE file
                copyfile(''.join(map(str,PP)), sendPPdir)  # send PP to directory of your file
            scf.write("\n\n") # formatting
            scf.write("ATOMIC_POSITIONS (crystal)\n")
            for index,x in enumerate(self.struct.frac_coords):
                scf.write("%s %f %f %f\n" %(str(self.struct.species[index]),x[0], x[1], x[2]))
            scf.write("\n\n")
            ###### CELL PARAMETERS ####
            #str(struct.lattice)
            scf.write("CELL_PARAMETERS (angstrom)\n")
            scf.write(str(self.struct.lattice) + "\n")
            scf.write("\n\n")
            ###KPOINTS####
            scf.write("K_POINTS (automatic)\n")
            scf.write("%i %i %i 0 0 0\n" %(8, 8, 8))
    #def WRITE_NSCF(self):
        with open(self.SEEDDIR + self.seedname + '.nscf.in','w',newline='\n') as nscf:
                nscf.write("&CONTROL\n")
                nscf.write("calculation ='nscf'\n")
                nscf.write("prefix = '" + self.seedname +"'\n")
                nscf.write("outdir = './bin'\n")
                nscf.write("pseudo_dir = './PP/'\n")
                nscf.write("verbosity='high'\n")
                nscf.write("/\n\n")
                ################# ATOM PARAM ###############
                nscf.write("&system\n")
                nscf.write("ibrav = 0\n")
                nscf.write("nat = %i\n" %(self.natom))
                nscf.write("ntyp = %i\n" %(self.ntypat))
                nscf.write("ecutwfc =%i  !Ryberg\n" %(self.Ecut))
                nscf.write("ecutrho =%i\n" %(self.Ecut*4))
                nscf.write("occupations = 'fixed'\n")
                nscf.write("nosym = .true.\n")
                nscf.write("nbnd = %i\n" %(self.NBND))
                nscf.write("!SOC\n")
                nscf.write("noncolin = %s\n" %(self.noncolin))
                nscf.write("lspinorb = %s\n" %(self.lspinorb))
                for index,atom in enumerate(self.UNIQUE_ATOMS):
                    nscf.write("starting_magnetization(%i) = %i\n" %(index+1,0))
                nscf.write("\n")
                ######## HUBBARD for HF ###########
                is_HF = False
                csvfile = open('HF.csv', newline='')
                HF = csv.reader(csvfile, delimiter=',')
                for index,atom in enumerate(self.UNIQUE_ATOMS):
                    for el in HF:
                        if str(atom) in str(el[0]) and len(str(atom)) == len(el[0]):
                            is_HF = True
                            nscf.write("Hubbard_U(%i) = %s\n" %(index+1, el[1])) ### U TERM
                            nscf.write("Hubbard_J0(%i) = %s\n" %(index+1, el[2])) ### J0 terM
                    ### reset csv reader
                    csvfile = open('HF.csv', newline='')
                    HF = csv.reader(csvfile, delimiter=',')
                if is_HF:
                    nscf.write("lda_plus_u = .TRUE.\n")
                    nscf.write("lda_plus_u_kind = 1\n")
                nscf.write("\n")
                nscf.write("/\n\n")
                #### ELECTRON CARD #####
                nscf.write("&ELECTRONS\n")
                nscf.write("conv_thr = 1.0d-7\n")
                nscf.write("mixing_beta = 0.495\n")
                nscf.write("mixing_mode = 'TF'\n")
                nscf.write("diagonalization= 'david'\n")
                nscf.write("adaptive_thr=.true.\n")
                nscf.write("/\n\n")
                nscf.write("ATOMIC_SPECIES\n")
                #### ATOMIC SPECIES CARD
                for index,atom in enumerate(self.UNIQUE_ATOMS):
                    nscf.write("%s " %(str(atom))) #Write element symbol to file
                    PP = gb.glob(self.pseudodir + str(atom) + "*") # get pseudopotential
                    PP = PP[0] # get first Psuedopotential found, for robustness
                    element = mg.Element(str(atom)) # Get element information
                    amass = float(element.atomic_mass) # Get atomic mass
                    nscf.write("%s " %(amass))    #write atomic mass
                    sendPPdir = self.SEEDDIR + "PP/" + ''.join(map(str,PP))[len(self.pseudodir):] # get directory to send PP
                    nscf.write(''.join(map(str,PP))[len(self.pseudodir):] + " \n") # write PP to QE file
                    copyfile(''.join(map(str,PP)), sendPPdir)  # send PP to directory of your file
                nscf.write("\n\n") # formatting
                nscf.write("ATOMIC_POSITIONS (crystal)\n")
                for index,x in enumerate(self.struct.frac_coords):
                    nscf.write("%s %f %f %f\n" %(str(self.struct.species[index]),x[0], x[1], x[2]))
                nscf.write("\n\n")
                ###### CELL PARAMETERS ####
                #str(struct.lattice)
                nscf.write("CELL_PARAMETERS (angstrom)\n")
                nscf.write(str(self.struct.lattice) + "\n")
                nscf.write("\n\n")
                nscf.write("K_POINTS (crystal)\n")
                nscf.write(kmesh.ESPRESSO(self.kpts[0],self.kpts[1],self.kpts[2]))
    #def write_w90(self):
        ########### PW2 WANNIER90 FILE #################
        with open(self.SEEDDIR + self.seedname + '.p2w.in','w',newline='\n') as f:
            f.write("&inputpp \n")
            f.write("outdir         = './bin'\n")
            f.write("prefix         = '%s'\n" %(self.seedname))
            f.write("seedname       = '%s'\n" %(self.seedname))
            f.write("spin_component = 'none'\n")
            f.write("write_mmn      = .true.\n")
            f.write("write_amn      = .true.\n")
            f.write("write_unk      = .false.\n")
            f.write("/\n")
        #########WANNIER90 FILE####################33
        with open(self.SEEDDIR + self.seedname + ".win",'w',newline='\n') as f:
            f.write("write_hr = .TRUE.\n")
            f.write("write_xyz = .TRUE.\n")
            f.write("!wannier_plot = .TRUE. \n")
            f.write("spinors = %s\n" %(self.lspinorb))
            f.write("num_wann = %i\n" %(self.wan))
            f.write("dis_num_iter=1000\n")
            f.write("trial_step=50\n")
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
    #def WRITE_CLEAN(self):
            with open(self.SEEDDIR + "/cleandft.sh", 'w',newline='\n') as file:
                file.write("rm *.amn\n")
                file.write("rm *.mmn\n")
                file.write("rm *.chk\n")
                file.write("rm *.eig\n")
                file.write("rm *.nnkp\n")
                file.write("rm *.nscf.out\n")
                file.write("rm *.scf.out\n")
                file.write("rm *.p2w.out\n")
                file.write("rm *.wout\n")
                file.write("rm *.xyz\n")
                file.write("rm *.dat\n")
                file.write("rm -rf bin\n")
    #def WRITE_SLURM(self):
            ####### SCF SLURM FILE ######################################
            with open(self.SEEDDIR + "scf.slurm",'w',newline= '\n') as f:
                f.write("#!/bin/bash\n")
                f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
                f.write("#SBATCH --time=10:00:00\n")
                f.write("#SBATCH --error=job.%J.err\n")
                f.write("#SBATCH --output=job.%J.out\n")
                f.write("#SBATCH --job-name=%s_SCF\n" %(self.seedname))
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
                f.write('mpirun -np %i pw.x <%s.scf.in> %s.scf.out\n' %(self.NCORE,self.seedname,self.seedname))
                f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
                f.write('#End Timestamp\n')
                f.write('date >> $benchmark_logs\n')
                f.write('echo "ulimit -l" >> $benchmark_logs\n')
                f.write('ulimit -l >> $benchmark_logs\n\n\n')
                f.write('#Cleanup\n')
                f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
                f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
                f.write('rm -rf $SLURM_JOB_ID\n')
                f.write('#sbatch nscf.slurm\n')
            ############# NSCF SLURM FILE #########################
            with open(self.SEEDDIR + "nscf.slurm",'w',newline= '\n') as f:
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
                f.write('mpirun -np %i pw.x <%s.nscf.in> %s.nscf.out\n' %(self.NCORE,self.seedname,self.seedname))
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
                ############ P2W SLURM FILE ###########################
            with open(self.SEEDDIR + "p2w.slurm",'w',newline= '\n') as f:
                f.write("#!/bin/bash\n")
                f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
                f.write("#SBATCH --time=10:00:00\n")
                f.write("#SBATCH --error=job.%J.err\n")
                f.write("#SBATCH --output=job.%J.out\n")
                f.write("#SBATCH --job-name=%s_P2W\n" %(self.seedname))
                f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
                f.write("#SBATCH --mail-type=ALL\n\n\n")
                f.write('#Select how logs get stored\n')
                f.write("mkdir $SLURM_JOB_ID\n")
                f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
                f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
                f.write('#Load Modules\n')
                f.write('ml load wannier90\n')
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
                f.write('wannier90.x -pp %s\n' %(self.seedname))
                f.write('mpirun -np %i pw2wannier90.x <%s.p2w.in> %s.p2w.out\n' %(self.NCORE,self.seedname,self.seedname))
                f.write('wannier90.x -pp %s\n' %(self.seedname))
                f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
                f.write('#End Timestamp\n')
                f.write('date >> $benchmark_logs\n')
                f.write('echo "ulimit -l" >> $benchmark_logs\n')
                f.write('ulimit -l >> $benchmark_logs\n\n\n')
                f.write('#Cleanup\n')
                f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
                f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
                f.write('rm -rf $SLURM_JOB_ID\n')
                f.write('#sbatch wann.slurm\n')
            ################ WANNER90 SLURM
            with open(self.SEEDDIR + "wann.slurm",'w',newline= '\n') as f:
                f.write("#!/bin/bash\n")
                f.write("#SBATCH --ntasks=1\n") ### _HR only works on single core operations
                f.write("#SBATCH --time=10:00:00\n")
                f.write("#SBATCH --error=job.%J.err\n")
                f.write("#SBATCH --output=job.%J.out\n")
                f.write("#SBATCH --job-name=%s_W90\n" %(self.seedname))
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
                f.write('wannier90.x %s\n' %(self.seedname))
                f.write('sleep 3 \n')
                f.write('grep CONV %s.wout >> wconv.txt\n' %(self.seedname))
                f.write('echo "Program is finished with exit code $? at: `date`"\n\n\n')
                f.write('#End Timestamp\n')
                f.write('date >> $benchmark_logs\n')
                f.write('echo "ulimit -l" >> $benchmark_logs\n')
                f.write('ulimit -l >> $benchmark_logs\n\n\n')
                f.write('#Cleanup\n')
                f.write('mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/\n')
                f.write('mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/\n')
                f.write('rm -rf $SLURM_JOB_ID\n')
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

# run = MAT2ESPRESSO()
# run.PMG2WAN()