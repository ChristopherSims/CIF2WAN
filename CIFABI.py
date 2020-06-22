#ABI
import os
import input
class CIF2ABINIT:
    def __init__(self):
        self.seedname = str(input.SEEDNAME)
        self.SEEDDIR = "OUTPUT/" + "abi-" + self.seedname + '/'  
        if not os.path.exists(self.SEEDDIR):
            os.makedirs(self.SEEDDIR)
    def WRITE_IN(self):
        with open(self.SEEDDIR + self.seedname + ".in",'w',newline='\n') as file:
                file.write("ndtset %i \n" % (2)) # for wanniertools
                file.write("eunit 1 #change to EV\n")
                file.write("prtvol 2 #print all kpoints\n")
                file.write("##### " + self.seedname + " stucture ######\n") # header
                file.write("acell" + a + b + c + " angstrom \n") # cell parameters
                file.write("angdeg" + alpha + beta + gamma + "\n") #angle param
                file.write("natom  " + str(natom) + "\n") #number of atoms
                ####
                file.write("xred\n") # Fractional coordinates
                for x in cord:
                   file.write("%12.8f%12.8f%12.8f\n" %(x[0],x[1],x[2]))
                 #####
                file.write("ntypat " + repr(ntypat) + "\n")

                file.write("typat " + typat + "\n")
                file.write("znucl ")
                for nu in np.unique(struct.atomic_numbers):
                    file.write(str(nu) + " ")
                file.write("\n")
                # enforce calculation of forces at each SCF step
                file.write("optforces 0 \n")

                ############ Plane wave  ########
                file.write("########### PLANE WAVE INFO ##########\n")
                file.write("ecut 20  ry \n")
                file.write("pawecutdg 22 ry\n")
                file.write("pawovlp  15\n")
                file.write("nstep 100\n")
                ### SCF STEP
                file.write("########### SCF STEP ##########\n")
                file.write("ngkpt1 4 4 4\n")
                file.write("nshiftk1 1 #just one shift is supported by wannier90\n")
                file.write("shiftk1 0.00   0.00   0.00  #no shift \n")
                file.write("tolvrs1  1.00d-15 \n")
                file.write("prtden1 1\n")
                file.write("kptopt1 1\n")
                file.write("#istwfk1 *1 #Controls the form of the wavefunctions\n")
                ############# NSCF STEP  #################
                #
                #
                #
                #############################################
                file.write("########### NSCF STEP ##########\n")
                file.write("prtvol2  1\n")
                file.write("iscf2 -2 \n")
                file.write("#nstep2 200\n")
                file.write("tolwfr2 1.d-15\n")
                file.write("getwfk2 -1\n")
                file.write("getden2 -1\n")
                file.write("prtden2 1\n")
                file.write("kptopt2 3\n")
                file.write("istwfk2  *1  #Controls the form of the wavefunctions\n")
                #### wannier 90 #####
                #
                #
                #
                #############################
                file.write("########### WANNER 90 ##########\n")
                file.write("#max_ncpus3 1\n")
                file.write("iscf3 -2 #nscf run\n")
                file.write("nstep3 0 #just read the old wave functions\n")
                file.write("tolwfr3 1.e-13 #dummy here\n")
                file.write("getwfk3 -1\n")
                file.write("getden3 -1   # Usual file handling data   \n")
                file.write("istwfk3 *1 #Controls the form of the wavefunctions\n")
                file.write("kptopt3 0   # Option for the automatic generation of k points  \n")
                file.write("nshiftk3 1\n")
                file.write("shiftk3 0.0 0.0 0.0\n")
                file.write("prtwant3 2   # Call to Wannier90\n")
                file.write("w90iniprj3 2\n")
                file.write("w90prtunk2 0   #Prints UNK files (for plotting the Wannier functions) \n")
        return(file)

    def WRITE_FILES(self):
        with open(self.SEEDDIR + self.seedname + ".files",'w',newline='\n') as file:
            file.write(self.seedname + ".in" + "\n")
            file.write(self.seedname + ".out" + "\n")
            file.write(self.seedname + "_i" + "\n")
            file.write(self.seedname + "_o" + "\n")
            file.write(self.seedname + "_1" + "\n")
        return(file)


##### DEBUG  ###################
# run = ABINIT()
# run.WRITE_FILES()
# run.WRITE_IN()