########### GENRATE MKGRID ###############
#
#
#
##########################################
def WANNIER(k1, k2, k3):
    mesh = "mp_grid = %i %i %i \n" % (k1,k2,k3)
    mesh += "begin kpoints\n"
    for x in range(k1):
        for y in range(k2):
            for z in range(k3):
                mesh = mesh +"%12.8f %12.8f %12.8f \n" % (x/k1 , y/k2 , z/k3)
    mesh+= "end kpoints"
    return(mesh)
def ESPRESSO(k1, k2, k3):
    tot = k1*k2*k3
    mesh = "%i\n" %(tot)
    for x in range(k1):
        for y in range(k2):
            for z in range(k3):
                mesh = mesh +"%12.8f %12.8f %12.8f %i\n" % (x/k1 , y/k2 , z/k3, 1)
    return(mesh)
def ABINIT(k1, k2, k3,num):
    tot = k1*k2*k3
    mesh = "nkpt%i %5i\n" %(num,tot)
    mesh += "kpt%i\n" %(num)
    for x in range(k1):
        for y in range(k2):
            for z in range(k3):
                mesh = mesh +"%12.8f %12.8f %12.8f\n" % (x/k1 , y/k2 , z/k3)
    return(mesh)
#msh = KMESH_ABINIT(2,2,2,2)
#print(msh)
