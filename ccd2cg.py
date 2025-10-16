#######################################
#                                     # 
#    Convert between CCD and CHARMM   # 
#                                     #
#######################################

# Import relevant functions
# -------------------------

import FuncConv
import argparse
import pandas as pd
import numpy as np
from copy import deepcopy
import sys, subprocess

# Get command line arguments
# ---------------------------

parser = argparse.ArgumentParser(description='Convert from output from co-folding (specified CCD codes and SMILES strings and userCCD codes) to Martini 3 CG representation. Please see keb721/CCD2MD on github for a list of the allowed SMILES strings and CCD codes.')

parser.add_argument('inputfile',  help='Input file name - .cif, .pdb or .gro.')
parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.')

parser.add_argument('-S',  '--SMILES',      help ='Used SMILES strings, list the order of the name of the ligand used. Note that when multiple of the same ligand are used this can be written either e.g. "POPE POPE" or "POPE 2".', nargs='+', default=[])

parser.add_argument('-nl', '--newlipidome', help='Use updated Martini 3 mappings for lipids as in DOI: 10.1021/acscentsci.5c00755 . Default off.', action='store_true')

CG_opts = parser.add_argument_group('martinize2 options')

# Add user-extensible martinize2 options
protein_network = CG_opts.add_mutually_exclusive_group()
protein_network.add_argument('-E', '--elastic', help='Use the elastic network for Martini when converting an atomisitic protein to coarse-grained. When not present, or only the flag is present only -elastic is passed to martinize2. Include any desired elastic network commands to be passed to martinize2 after this argument (this may but does not need to include -elastic).', action='store_true')
protein_network.add_argument('-G', '--go',      help='Use the Go network for Martini when converting an atomistic protein to coarse-grained. Include any required/desired commands to be passed to martinize2 after this argument.', action='store_true')

CG_opts.add_argument('-M', '--martinize', help='Add additional arguments to martinize2.', action='store_true')


sys_opts = parser.add_argument_group('membrane-embedded system options')

sys_opts.add_argument('-mem', '--membrane', help='Embed the converted system into a membrane. Note that this will lead to minor rearragements of any ligands. After this flag, it is possible to add the arguments for the upper and lower leaflet in the form "-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1" where the codes represent lipids and the numbers represent the ratio between them. The default is two leaflets of pure POPC.', action='store_true')
sys_opts.add_argument('-C',  '--conc',      help='Concentration of NaCl in system - charge balance is maintained, Default = 0.15.', default = 0.15)
sys_opts.add_argument('-mp',  '--mempro',   help='Additional arguments for embedding the protein in the membrane using MemPrO. Add any additional arguments after this flag - default 5 grid points and 15 minimisation operations.', action='store_true')

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

input_data   = []
element_name = True
for atom in tmp:
    try:
        if atom['elem'] != 'H' and atom['elem'][0] != 'H':
            input_data.append(atom)
    except (KeyError, IndexError):
        # Either some or all element names missing
        if element_name:
            print('# WARNING: Element names are missing - attempting to infer from atom names. Note that this may cause issues especially when converting to CG.')
            element_name = False
        atom['elem'] = FuncConv.determine_element(atom['name'])
        if atom['elem'] != 'H':
            input_data.append(atom)               

input_data = pd.DataFrame.from_dict(input_data, orient='columns')

# Split SMILES strings
# ---------------------

if len(args.SMILES) < 2:
    # 0/1 SMILES strings
    SMILES = args.SMILES
else:
    SMILES   = []
    in_range = True ; i = 0
    while in_range:
        try:
            num = int(args.SMILES[i+1])
            SMILES.extend([args.SMILES[i]]*num)
            i += 2
        except (ValueError, IndexError) as e:
            SMILES.append(args.SMILES[i])
            i += 1
        in_range = False if i >= len(args.SMILES) else True
        
print('# INFO: Assuming that chains are labelled sequentially.')

# Find residues to reorder
# ------------------------
    
ligands, atoms, IDs, types, locs = FuncConv.get_residues(input_data, 'CCD', SMILES, False)

# Convert system to CG    
prot = True if len(np.concatenate([np.array(i) for i in IDs])) != len(input_data) else False    

# Generate the martinize2 parameters
mart_params = FuncConv.get_CG_params(np.array(sys.argv), args.martinize, args.elastic, args.go) if prot else []


output_data = FuncConv.to_CG(args.inputfile, args.outputfile, input_data, mart_params, ligands, input_data,
                             types, locs, prot, newlipidome=args.newlipidome)

# Write output
# -------------

basename = '.'.join(args.outputfile.split('.')[:-1])

PDBfile = basename+'_nomem.pdb' if args.membrane else args.outputfile
print('# INFO: Assuming that an entirely new PDB file is being written rather than modified in place')    
FuncConv.write_PDB(PDBfile, output_data, title=title, cryst=cryst, ligand_chains=False)

# Embed in membrane
# -----------------

if args.membrane:
    FuncConv.build_membrane_CG(ligands, PDBfile, args.outputfile, command_line, args.mempro,
                               args.membrane, args.conc, newlipidome=args.newlipidome)
    mempro_output, title, cryst = FuncConv.read_GRO(basename+'_MemPrO/Rank_1/CG_System_rank_1/CG-system.gro')
    
    FuncConv.write_PDB(args.outputfile, mempro_output, title=title, cryst=cryst, ligand_chains=False)

FuncConv.get_topology_CG(args.outputfile, args.membrane, ligands, prot, args.inputfile, newlipidome=args.newlipidome)

