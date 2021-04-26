#SIESTA
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
class PMG2SIESTA:
    def __init__(self,input):
        #import structure from PYMATGEN
        self.key = input.APIKEY # Import API key from input
        self.MaterialID = input.MATERIALID # Set class variable MaterialID from input
        self.cutoff = input.ECUT
        self.Ecut = input.ECUT #Class variable; in Ryberg, transfer from input file
        self.kpts = input.KMESH ### KPTS from input file for a generator, class variable
        #self.pseudodir = "ABINIT/ABINITPP/ATOMICDATA/" # Set the directory for pseudopotentials
        self.EMAIL = input.Email #Class variable; set email for SLURM from input file
        self.NCORE = input.ncore #Class Variable; number of cores, taken from input file
        self.SOC = input.SOC
        self.nbnd = input.NUMBANDS
        self.lspinorb = '.FALSE.'
        if input.SOC:
            self.noncolin = '.TRUE.'
            self.lspinorb = '.TRUE.'
        selectconventional = input.conventional_cell
        selectrelax = input.Relax
        with matproj.MPRester(self.key) as m:
            struct = m.get_structure_by_material_id(self.MaterialID,final=selectrelax,conventional_unit_cell=selectconventional) # GET STRUCT FROM MATERIALS PROJECT
            self.struct = m.get_structure_by_material_id(self.MaterialID,final=selectrelax,conventional_unit_cell=selectconventional) # GET STRUCT FROM MATERIALS PROJECT
        C = Composition(str(struct.formula))
        self.seedname = str(input.seedname) # set class variable SEEDNAME from input file
        self.SEEDDIR = "OUTPUT/" + 'S-' + self.seedname + '/'  #set the seed directory
        # MAKE DIRECTORIES
        if not os.path.exists(self.SEEDDIR):
            os.makedirs(self.SEEDDIR)
        if not os.path.exists(self.SEEDDIR + "WT"):
            os.makedirs(self.SEEDDIR + "WT")
        ##### Set variables
        self.UNIQUE_ATOMS = np.unique(struct.species)
        self.cord = struct.frac_coords
        self.ntypat = np.unique(struct.atomic_numbers).size
        self.natom = np.array(struct.atomic_numbers).size
        self.numlist = np.array(struct.atomic_numbers)[::1]
        self.atomnum = np.unique(struct.atomic_numbers)[::1]
        self.typat = ""
        for x in struct.atomic_numbers:
            for at in self.atomnum:
                if at == x:
                        self.typat += ''.join(map(str,np.where(self.atomnum == at)
                            [0] + 1)) + " "
        # SOC CASE
        num_bands = input.NUMBANDS
        # for x in struct.species:
        #     PP = gb.glob(self.pseudodir + str(x) + ".*") 
        #     readpp = minidom.parse(''.join(map(str,PP)))
        #     items = readpp.getElementsByTagName('atom')
        #     num_bands += float(items[0].attributes['valence'].value)
        num_wann = num_bands - (num_bands%2) # FORCE EVEN
        self.wan = num_wann
        if self.SOC:
            self.WTwan = self.wan/2
        else:
            self.WTwan = self.wan
    def generate(self):
        with open(self.SEEDDIR + self.seedname + ".fdf",'w',newline='\n') as f:
            f.write("#FDF for SIESTA\n")
            f.write("SystemName%10s\n" %(self.seedname))
            f.write("SystemLabel%10s\n" %(self.seedname))
            f.write("NumberOfSpecies%10i\n"%(self.UNIQUE_ATOMS.size))
            f.write("NumberOfAtoms%10i\n\n\n"%(self.natom))
            f.write("#\n")
            f.write("User-Basis     .false.\n")
            f.write("#\n\n\n")
            f.write("# Output options\n")
            f.write("WriteCoorStep\n")
            f.write("WriteMullikenPop      1\n\n\n")
            f.write("%"+"block ChemicalSpeciesLabel\n")
            # find the number of each unique atoms
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            atomlabel = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.struct.atomic_numbers == self.atomnum[i])).size
                atomlabel[i] = i+1
            for i,x in enumerate(self.UNIQUE_ATOMS):
                f.write("%5i%5i%5s\n" %(i+1,x.Z,x.symbol)) # label, atomic number, symbol
            f.write("%"+"endblock ChemicalSpeciesLabel\n\n")
            f.write("LatticeConstant 1 Ang\n\n")
            f.write("%"+ "block LatticeVectors\n")
            f.write(str(self.struct.lattice) + "\n")
            f.write("%" + "endblock LatticeVectors\n\n\n")
            f.write("AtomicCoordinatesFormat ScaledCartesian\n\n")
            f.write("%" + "block AtomicCoordinatesAndAtomicSpecies\n")
            prtlabel = np.zeros(self.numlist.size)
            for index,x in enumerate(self.struct.frac_coords):
                for i,Z in enumerate(np.unique(self.numlist)):
                    if Z == self.numlist[index]:
                        prtlabel[index] = atomlabel[i]
            for index,x in enumerate(self.struct.frac_coords): 
                f.write("%f %f %f %i\n" %(x[0], x[1], x[2], int(prtlabel[index])))
            f.write("%" + "endblock AtomicCoordinatesAndAtomicSpecies\n")
            f.write("PAO.SplitNorm 0.20\n")
            f.write("PAO.EnergyShift 0.1 eV\n\n\n")
            f.write("\n\n\n#SCF INFO\n")
            #f.write("kgrid_cutoff %5i Ang\n" %(5))
            f.write("MeshCutoff %5i Ry \n" %(self.cutoff))
            f.write("MaxSCFIterations %5i \n" %(100))
            f.write("DM.MixingWeight %5.1f \n" %(0.7))    
            f.write("DM.NumberPulay %5i\n" %(3))
            f.write("DM.Tolerance %5s\n" %("1.d-9"))
            f.write("DM.UseSaveDM %5s .true.\n\n\n")
            f.write("xc.functional %5s GGA %5s# Exchange-correlation functional\n")
            f.write("xc.authors%5sPBE%5s# Exchange-correlation version\n")
            if self.SOC:
                f.write("SpinPolarized%5strue%5s# Logical parameters are: yes or no\n")
            else:
                f.write("SpinPolarized%5sfalse%5s# Logical parameters are: yes or no\n")
            f.write("SolutionMethod%5sdiagon%5s # OrderN or Diagon\n")
            f.write("ElectronicTemperature%5s25 meV%5s# Temp. for Fermi smearing\n")
            f.write("Diag.DivideAndConquer .true.\n")
            f.write("Diag.ParallelOverK .true.\n")
            f.write("WriteEigenvalues .true.\n")
            f.write("WriteKbands .true.\n")
            f.write("WriteBands .true.\n")
            f.write("%" + "block kgrid_Monkhorst_Pack\n")
            f.write("%3i%3i%3i%3i\n" %(self.kpts[0],0,0,0))
            f.write("%3i%3i%3i%3i\n" %(0,self.kpts[1],0,0))
            f.write("%3i%3i%3i%3i\n" %(0,0,self.kpts[2],0))
            f.write("%" + "endblock kgrid_Monkhorst_Pack\n\n\n")

            f.write("\n\n\n")
            f.write("#\n")
            f.write("# Variables related with the Wannierization of the manifolds\n")
            f.write("NumberOfBandManifoldsForWannier 1\n")
            f.write("%"+"block WannierProjections\n")
            f.write("1 # Sequential index of the manifold, from 1 to NumberOfBandManifoldsForWannier\n")
            f.write("%i %3i # Indices of the initial and final band of the manifold\n" %(1,self.nbnd))
            f.write("%i # Number of bands for Wannier transformation\n" %(self.nbnd))
            f.write("num_iter 0 # Number of iterations for the minimization of \Omega\n")
            f.write("write_hr # Write the Hamiltonian in the WF basis\n")
            f.write("%"+"endblock WannierProjections\n\n")
            f.write("%"+"block kMeshforWannier\n")
            f.write("%3i%3i%3i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("%"+"endblock kMeshforWannier\n\n\n")
            ########################
            f.write("###Siesta2wannier inputs###\n")
            f.write("Siesta2Wannier90.WriteMmn .true.\n")
            f.write("Siesta2Wannier90.WriteAmn .true.\n")
            f.write("Siesta2Wannier90.WriteEig .true.\n")
            f.write("Siesta2Wannier90.NumberOfBands%5i\n" %(self.nbnd))
            f.write("\n")
            
    #WRITING SLURM FILE
        with open(self.SEEDDIR + "siesta.slurm",'w',newline= '\n') as f:
            f.write("#!/bin/bash\n")
            #f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --time=5:00:00\n")
            f.write("#SBATCH --error=job.%J.err\n")
            f.write("#SBATCH --output=job.%J.out\n")
            f.write("#SBATCH --job-name=%s_SIESTA\n" %(self.seedname))
            f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
            f.write("#SBATCH --mail-type=NONE\n\n\n")
            f.write('#Select how logs get stored\n')
            f.write("mkdir $SLURM_JOB_ID\n")
            f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
            f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
            f.write('#Load Modules\n')
            f.write('ml load wannier90\n')
            f.write('ml load siesta\n\n\n')
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
            f.write('wannier90.x -pp %s' %(self.seedname))
            f.write('sleep 3\n')
            f.write('mpirun -np %i siesta <%s.fdf> %s.out\n' %(self.NCORE,self.seedname,self.seedname))
            f.write('sleep 3\n')
            f.write('mv %s.eigW %s.eig\n' %(self.seedname,self.seedname))
            f.write('sleep 3\n')
            f.write('wannier90.x %s' %(self.seedname))
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
    #WRITING CLEAN DFT FILE
        with open(self.SEEDDIR + "cleandft.sh",'w',newline='\n') as f:
            f.write("rm INPUT_TMP.*\n")
            f.write("rm *.ion\n")
            f.write("rm *.log\n")
            f.write("rm *ion.xml\n")
            f.write("rm *.XV\n")
            f.write("rm *.FA\n")
            f.write("rm *.DM\n")
            f.write("rm *.ORB_INDX\n")
            f.write("rm OCCS\n")
            f.write("rm MESSAGES\n")
            f.write("rm *.BONDS\n")
            f.write("rm *.BONDS_FINAL\n")
            f.write("rm *.KP\n")
            f.write("rm *.STRUCT_OUT\n")
            f.write("rm *.xml\n")
            f.write("rm *.alloc\n")
            f.write("rm FORCE_STRESS\n")
            f.write("rm NON_TRIMMED_KP_LIST\n")
            f.write("rm CLOCK\n")
            f.write("rm BASIS_ENTHALPY\n")
            f.write("rm BASIS_HARRIS_ENTHALPY\n")
            f.write("rm Wconv.txt\n")
            f.write("rm *.werr\n")
            f.write("rm *.wout\n")
            f.write("rm *.nnkp\n")
            f.write("rm *.EIG\n")
            f.write("rm *.xyz\n")
            f.write("rm *_wsvec.dat\n")
            f.write("rm *.mmn\n")
            f.write("rm *.eig\n")
            f.write("rm *.amn\n")
            f.write("rm *.err\n")
            f.write("rm *_NORMAL_EXIT\n")
            f.write("rm *.DOS \n")
            f.write("rm *.chk \n")
            f.write("rm *.PDOS\n")
            f.write("rm *.out\n")
            f.write("rm *.dat\n")
        #########WANNIER90 FILE####################33
        with open(self.SEEDDIR + self.seedname + ".win",'w',newline='\n') as f:
            f.write("write_hr = .TRUE.\n")
            f.write("write_xyz = .TRUE.\n")
            f.write("!wannier_plot = .TRUE. \n")
            f.write("spinors = %s\n" %(self.lspinorb))
            f.write("num_wann = %i\n" %(self.wan))
            f.write("dis_num_iter=1000\n")
            f.write("!trial_step=50\n")
            f.write("num_iter = 1000\n\n\n\n")
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
    #WRITE Wanniertools file
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
#######################################
#
#
#
#
#
#
#
#
#
#######################################
from ase.io import read, write
from ase.optimize import BFGS
from ase.visualize import view
from ase.io import write as ASEwrite
class CIF2SIESTA:
    def __init__(self,input):
        crystal = read(input.CIF_file)
        self.key = input.APIKEY # Import API key from input
        self.MaterialID = input.MATERIALID # Set class variable MaterialID from input
        self.cutoff = input.ECUT
        self.Ecut = input.ECUT #Class variable; in Ryberg, transfer from input file
        self.kpts = input.KMESH ### KPTS from input file for a generator, class variable
        #self.pseudodir = "ABINIT/ABINITPP/ATOMICDATA/" # Set the directory for pseudopotentials
        self.EMAIL = input.Email #Class variable; set email for SLURM from input file
        self.NCORE = input.ncore #Class Variable; number of cores, taken from input file
        self.SOC = input.SOC
        self.nbnd = input.NUMBANDS
        self.lspinorb = '.FALSE.'
        if input.SOC:
            self.noncolin = '.TRUE.'
            self.lspinorb = '.TRUE.'
        selectconventional = input.conventional_cell
        selectrelax = input.Relax
        self.seedname = str(input.seedname) # set class variable SEEDNAME from input file
        self.SEEDDIR = "OUTPUT/" + 'S-' + self.seedname + '/'  #set the seed directory
        # MAKE DIRECTORIES
        if not os.path.exists(self.SEEDDIR):
            os.makedirs(self.SEEDDIR)
        if not os.path.exists(self.SEEDDIR + "WT"):
            os.makedirs(self.SEEDDIR + "WT")
        ##### Set variables
        self.cell = crystal.get_cell()
        self.cord = crystal.get_scaled_positions(wrap=True)
        self.lattice = np.array(self.cell.lengths())
        lattice_old = pd.DataFrame(np.array(self.cell))
        self.lattice = lattice_old.to_string(header=False, index=False)
        self.species = crystal.get_chemical_symbols()
        self.atomic_numbers = crystal.get_atomic_numbers()
        self.ntypat = np.unique(self.atomic_numbers).size
        self.natom = np.array(self.atomic_numbers).size
        self.Z = np.unique(crystal.get_atomic_numbers())
        self.UNIQUE_ATOMS = np.unique(self.species)
        self.atomnum = np.unique(self.atomic_numbers)[::1]
        self.numlist = np.array(self.atomic_numbers)
        #self.UNIQUE_ATOMS = np.unique(struct.species)
        #self.cord = struct.frac_coords
        #self.ntypat = np.unique(struct.atomic_numbers).size
        #self.natom = np.array(struct.atomic_numbers).size
        #self.numlist = np.array(struct.atomic_numbers)[::1]
        #self.atomnum = np.unique(struct.atomic_numbers)[::1]
        self.typat = ""
        for x in self.atomic_numbers:
            for at in self.atomnum:
                if at == x:
                        self.typat += ''.join(map(str,np.where(self.atomnum == at)
                            [0] + 1)) + " "
        # SOC CASE
        num_bands = input.NUMBANDS
        # for x in struct.species:
        #     PP = gb.glob(self.pseudodir + str(x) + ".*") 
        #     readpp = minidom.parse(''.join(map(str,PP)))
        #     items = readpp.getElementsByTagName('atom')
        #     num_bands += float(items[0].attributes['valence'].value)
        num_wann = num_bands - (num_bands%2) # FORCE EVEN
        self.wan = num_wann
        if self.SOC:
            self.WTwan = self.wan/2
        else:
            self.WTwan = self.wan
    def generate(self):
        with open(self.SEEDDIR + self.seedname + ".fdf",'w',newline='\n') as f:
            f.write("#FDF for SIESTA\n")
            f.write("SystemName%10s\n" %(self.seedname))
            f.write("SystemLabel%10s\n" %(self.seedname))
            f.write("NumberOfSpecies%10i\n"%(self.UNIQUE_ATOMS.size))
            f.write("NumberOfAtoms%10i\n\n\n"%(self.natom))
            f.write("#\n")
            f.write("User-Basis     .false.\n")
            f.write("#\n\n\n")
            f.write("# Output options\n")
            f.write("WriteCoorStep\n")
            f.write("WriteMullikenPop      1\n\n\n")
            f.write("%"+"block ChemicalSpeciesLabel\n")
            # find the number of each unique atoms
            typnum = np.zeros(self.UNIQUE_ATOMS.size)
            atomlabel = np.zeros(self.UNIQUE_ATOMS.size)
            for i in range(0,self.UNIQUE_ATOMS.size):
                typnum[i] = np.asarray(np.where(self.atomic_numbers == self.atomnum[i])).size
                atomlabel[i] = i+1
            for i,x in enumerate(self.UNIQUE_ATOMS):
                f.write("%5i%5i%5s\n" %(i+1,self.Z[i],x)) # label, atomic number, symbol
            f.write("%"+"endblock ChemicalSpeciesLabel\n\n")
            f.write("LatticeConstant 1 Ang\n\n")
            f.write("%"+ "block LatticeVectors\n")
            f.write(str(self.lattice) + "\n")
            f.write("%" + "endblock LatticeVectors\n\n\n")
            f.write("AtomicCoordinatesFormat ScaledCartesian\n\n")
            f.write("%" + "block AtomicCoordinatesAndAtomicSpecies\n")
            prtlabel = np.zeros(self.natom)
            for index,x in enumerate(self.cord):
                for i,Z in enumerate(np.unique(self.numlist)):
                    if Z == self.numlist[index]:
                        prtlabel[index] = atomlabel[i]
            for index,x in enumerate(self.cord): 
                f.write("%f %f %f %i\n" %(x[0], x[1], x[2], int(prtlabel[index])))
            f.write("%" + "endblock AtomicCoordinatesAndAtomicSpecies\n")
            f.write("PAO.SplitNorm 0.20\n")
            f.write("PAO.EnergyShift 0.1 eV\n\n\n")
            f.write("\n\n\n#SCF INFO\n")
            #f.write("kgrid_cutoff %5i Ang\n" %(5))
            f.write("MeshCutoff %5i Ry \n" %(self.cutoff))
            f.write("MaxSCFIterations %5i \n" %(100))
            f.write("DM.MixingWeight %5.1f \n" %(0.7))    
            f.write("DM.NumberPulay %5i\n" %(3))
            f.write("DM.Tolerance %5s\n" %("1.d-9"))
            f.write("DM.UseSaveDM %5s .true.\n\n\n")
            f.write("xc.functional %5s GGA %5s# Exchange-correlation functional\n")
            f.write("xc.authors%5sPBE%5s# Exchange-correlation version\n")
            if self.SOC:
                f.write("SpinPolarized%5strue%5s# Logical parameters are: yes or no\n")
            else:
                f.write("SpinPolarized%5sfalse%5s# Logical parameters are: yes or no\n")
            f.write("SolutionMethod%5sdiagon%5s # OrderN or Diagon\n")
            f.write("ElectronicTemperature%5s25 meV%5s# Temp. for Fermi smearing\n")
            f.write("Diag.DivideAndConquer .true.\n")
            f.write("Diag.ParallelOverK .true.\n")
            f.write("WriteEigenvalues .true.\n")
            f.write("WriteKbands .true.\n")
            f.write("WriteBands .true.\n")
            f.write("%" + "block kgrid_Monkhorst_Pack\n")
            f.write("%3i%3i%3i%3i\n" %(self.kpts[0],0,0,0))
            f.write("%3i%3i%3i%3i\n" %(0,self.kpts[1],0,0))
            f.write("%3i%3i%3i%3i\n" %(0,0,self.kpts[2],0))
            f.write("%" + "endblock kgrid_Monkhorst_Pack\n\n\n")

            f.write("\n\n\n")
            f.write("#\n")
            f.write("# Variables related with the Wannierization of the manifolds\n")
            f.write("NumberOfBandManifoldsForWannier 1\n")
            f.write("%"+"block WannierProjections\n")
            f.write("1 # Sequential index of the manifold, from 1 to NumberOfBandManifoldsForWannier\n")
            f.write("%i %3i # Indices of the initial and final band of the manifold\n" %(1,self.nbnd))
            f.write("%i # Number of bands for Wannier transformation\n" %(self.nbnd))
            f.write("num_iter 0 # Number of iterations for the minimization of \Omega\n")
            f.write("write_hr # Write the Hamiltonian in the WF basis\n")
            f.write("%"+"endblock WannierProjections\n\n")
            f.write("%"+"block kMeshforWannier\n")
            f.write("%3i%3i%3i\n" %(self.kpts[0],self.kpts[1],self.kpts[2]))
            f.write("%"+"endblock kMeshforWannier\n\n\n")
            ########################
            f.write("###Siesta2wannier inputs###\n")
            f.write("Siesta2Wannier90.WriteMmn .true.\n")
            f.write("Siesta2Wannier90.WriteAmn .true.\n")
            f.write("Siesta2Wannier90.WriteEig .true.\n")
            f.write("Siesta2Wannier90.NumberOfBands%5i\n" %(self.nbnd))
            f.write("\n")
            
    #WRITING SLURM FILE
        with open(self.SEEDDIR + "siesta.slurm",'w',newline= '\n') as f:
            f.write("#!/bin/bash\n")
            #f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --ntasks=%i\n" %(self.NCORE))
            f.write("#SBATCH --time=5:00:00\n")
            f.write("#SBATCH --error=job.%J.err\n")
            f.write("#SBATCH --output=job.%J.out\n")
            f.write("#SBATCH --job-name=%s_SIESTA\n" %(self.seedname))
            f.write('#SBATCH --mail-user=%s\n' %(self.EMAIL))
            f.write("#SBATCH --mail-type=NONE\n\n\n")
            f.write('#Select how logs get stored\n')
            f.write("mkdir $SLURM_JOB_ID\n")
            f.write('export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n')
            f.write('export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"\n\n\n')
            f.write('#Load Modules\n')
            f.write('ml load wannier90\n')
            f.write('ml load siesta\n\n\n')
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
            f.write('wannier90.x -pp %s' %(self.seedname))
            f.write('sleep 3\n')
            f.write('mpirun -np %i siesta <%s.fdf> %s.out\n' %(self.NCORE,self.seedname,self.seedname))
            f.write('sleep 3\n')
            f.write('mv %s.eigW %s.eig\n' %(self.seedname,self.seedname))
            f.write('sleep 3\n')
            f.write('wannier90.x %s' %(self.seedname))
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
    #WRITING CLEAN DFT FILE
        with open(self.SEEDDIR + "cleandft.sh",'w',newline='\n') as f:
            f.write("rm INPUT_TMP.*\n")
            f.write("rm *.ion\n")
            f.write("rm *.log\n")
            f.write("rm *ion.xml\n")
            f.write("rm *.XV\n")
            f.write("rm *.FA\n")
            f.write("rm *.DM\n")
            f.write("rm *.ORB_INDX\n")
            f.write("rm OCCS\n")
            f.write("rm MESSAGES\n")
            f.write("rm *.BONDS\n")
            f.write("rm *.BONDS_FINAL\n")
            f.write("rm *.KP\n")
            f.write("rm *.STRUCT_OUT\n")
            f.write("rm *.xml\n")
            f.write("rm *.alloc\n")
            f.write("rm FORCE_STRESS\n")
            f.write("rm NON_TRIMMED_KP_LIST\n")
            f.write("rm CLOCK\n")
            f.write("rm BASIS_ENTHALPY\n")
            f.write("rm BASIS_HARRIS_ENTHALPY\n")
            f.write("rm Wconv.txt\n")
            f.write("rm *.werr\n")
            f.write("rm *.wout\n")
            f.write("rm *.nnkp\n")
            f.write("rm *.EIG\n")
            f.write("rm *.xyz\n")
            f.write("rm *_wsvec.dat\n")
            f.write("rm *.mmn\n")
            f.write("rm *.eig\n")
            f.write("rm *.amn\n")
            f.write("rm *.err\n")
            f.write("rm *_NORMAL_EXIT\n")
            f.write("rm *.DOS \n")
            f.write("rm *.chk \n")
            f.write("rm *.PDOS\n")
            f.write("rm *.out\n")
            f.write("rm *.dat\n")
        #########WANNIER90 FILE####################33
        with open(self.SEEDDIR + self.seedname + ".win",'w',newline='\n') as f:
            f.write("write_hr = .TRUE.\n")
            f.write("write_xyz = .TRUE.\n")
            f.write("!wannier_plot = .TRUE. \n")
            f.write("spinors = %s\n" %(self.lspinorb))
            f.write("num_wann = %i\n" %(self.wan))
            f.write("dis_num_iter=1000\n")
            f.write("!trial_step=50\n")
            f.write("num_iter = 1000\n\n\n\n")
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
        #WRITE Wanniertools file
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
if __name__ == "__main__":
    import sys, os
    sys.path
    parentfolder = os.path.dirname(os.path.dirname( __file__ ))
    sys.path.append(parentfolder)