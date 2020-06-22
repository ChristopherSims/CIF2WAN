#MAIN
import input
from ESPRESSO import MAT2ESPRESSO
from ABINIT import MAT2ABINIT
from VASP import MAT2VASP
from SIESTA import MAT2SIESTA

if (input.DFT == 'ESPRESSO'):
    run = MAT2ESPRESSO()
    run.PMG2WAN()
if (input.DFT == 'ABINIT'):
    run = MAT2ABINIT()
    run.PMG2WAN()
if (input.DFT == 'VASP'):
    run = MAT2VASP()
    run.PMG2WAN()
if (input.DFT == 'SIESTA'):
    run = MAT2SIESTA()
    run.PMG2WAN()