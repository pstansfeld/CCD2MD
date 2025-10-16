#######################################
#                                     # 
#    Convert between CCD and CHARMM   # 
#                                     #
#######################################

# Import relevant functions
# -------------------------

import FuncConv
import argparse
import numpy as np
import pandas as pd
from copy import deepcopy
import sys, subprocess, os

# Get command line arguments
# ---------------------------

parser = argparse.ArgumentParser(description='Embed a CHARMM system in a membrane - please note that this requires coarse-graining. Please see keb721/CCD2AT on github for a list of the allowed (CG-mapped) molecules.')

parser.add_argument('inputfile', help='Input file name - .cif, .pdb or .gro.')
parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.')

parser.add_argument('-lc', '--ligchain', help='Output ligands in their own chains - default is off.', action='store_true')

sys_opts = parser.add_argument_group('membrane-embedded system options')

sys_opts.add_argument('-mem', '--membrane', help='After this flag, it is possible to add the arguments for the upper and lower leaflet in the form "-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1" where the codes represent lipids and the numbers represent the ratio between them. The default is two leaflets of pure POPC. Note that membrane embedding may result in minor rearrangements of ligands.', action='store_true')
sys_opts.add_argument('-C',  '--conc',      help='Concentration of NaCl in system - charge balance is maintained, Default = 0.15.', default = 0.15)
sys_opts.add_argument('-mp',  '--mempro',   help='Additional arguments for embedding the protein in the membrane using MemPrO. Add any additional arguments after this flag - default 5 grid points and 15 minimisation operations.', action='store_true')
sys_opts.add_argument('-at', '--cg2at',     help='Additional arguments to pass to CG2AT. Add additional arguments after this flag.', action='store_true')
# sys_opts.add_argument('-sol', '--water',    help='The water model to use - currently only TIP3P') 

info = parser.add_argument_group()
info.add_argument('-V', '--Version', action='version', version='Version '+FuncConv.__version__)


args, unknownargs = parser.parse_known_args()
command_line      = np.array(sys.argv)

# Read input data
# ---------------

if args.inputfile[-3:] == 'cif':
    tmp, title = FuncConv.read_CIF(args.inputfile)
    cryst = []
elif args.inputfile[-3:] == 'gro':
    tmp, title, cryst = FuncConv.read_GRO(args.inputfile)
else:
    if args.inputfile[-3:] != 'pdb':
        print('# WARNING: assuming that input file {} is written in PDB style despite file extension'.format(args.inputfile))
    tmp, title, cryst = FuncConv.read_PDB(args.inputfile)

print('# INFO: Any Hs present will be removed')

input_data   = []
element_name = True
for atom in tmp:
    try:
        if atom['elem'] != 'H' and atom['elem'][0] != 'H':
            input_data.append(atom)
    except (KeyError, IndexError):
        # Either some or all element names missing
        if element_name:
            print('# WARNING: Element names are missing - attempting to infer from atom names. Note that this may cause issues.')
            element_name = False
        atom['elem'] = FuncConv.determine_element(atom['name'])
        if atom['elem'] != 'H':
            input_data.append(atom)                

input_data = pd.DataFrame.from_dict(input_data, orient='columns')

# Find ligands to omit
# ---------------------
    
ligands, atoms, IDs, types, locs = FuncConv.get_residues(input_data, 'CHARMM', [], args.ligchain)

output_data = input_data.to_dict('records')


# Write output
# -------------

if args.outputfile.count('/') != 0:
    output_dir = '/'.join(args.outputfile.split('/')[:-1])
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
            

basename = '.'.join(args.outputfile.split('.')[:-1])

PDBfile = basename+'_nomem.pdb'
FuncConv.write_PDB(PDBfile, output_data, title=title, cryst=cryst, ligand_chains=args.ligchain)

# =================== #
# Embed into membrane #
# =================== # 

# Convert system to CG
# --------------------

prot        = True if len(np.concatenate([np.array(i) for i in IDs])) != len(input_data) else False    
CG_output   = basename+'_CG'
output_data = FuncConv.to_CG(args.inputfile, CG_output+'.pdb', input_data, ['-elastic'], ligands, input_data, types, locs, prot)

CG_output = basename+'_CG_system.pdb'

FuncConv.write_PDB(CG_output, output_data, title=title, cryst=cryst, ligand_chains=False)

# Embed in membrane
# -----------------
    
FuncConv.build_membrane_CG(ligands, CG_output, args.outputfile, command_line, args.mempro, args.membrane, args.conc) 

# Convert back to atomistic
# -------------------------
    
FuncConv.convert_membrane_at(output_data, basename, command_line, args.cg2at)

subprocess.run(['scp', args.outputfile.split('.')[0]+'_CG2AT/FINAL/final_cg2at_aligned.pdb', args.outputfile])


FuncConv.get_topology_atomistic(args.outputfile, True)

