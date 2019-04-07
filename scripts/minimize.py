#coding=utf-8
import sys
import os
import re
import pymol
import optimize_plugin as opti

#parameters
path_input = sys.argv[1]
path_output = sys.argv[2]

#call pymol
pymol.pymol_argv = ['pymol','-qc']
pymol.finish_launching()
#load input structure
pymol.cmd.load(path_input, pdb_code)
#local energy minimization
opti.minimize(pdb_code, "MMFF94s", "Conjugate Gradients", 100, 0.001, True, 6.0, 8.0)
#save output structure
pymol.cmd.save(path_output,"(all)",-1)
#quit pymol
pymol.quit()
