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

parser = argparse.ArgumentParser(description='Convert from output from co-folding (specified CCD codes and SMILES strings, and userCCD codes) to CHARMM. Please see keb721/CCD2MD on github for a list of the allowed SMILES strings and CCD codes.')

parser.add_argument('inputfile',  help='Input file name - .cif, .pdb or .gro.')
parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.')

parser.add_argument('-lc', '--ligchain', help='Output ligands in their own chains - default is off.', action='store_true')
parser.add_argument('-S',  '--SMILES',   help='Used SMILES strings, list the order of the name of the ligand used. Note that when multiple of the same ligand are used this can be written either e.g. "POPE POPE" or "POPE 2".', nargs='+', default=[])

parser.add_argument('-gh', '--pdb2gmx', help='Override ALL defaults of pdb2gmx and optionally pass extra arguments - note that this may require interactivity, and may be necessary for a starting MET. Default is topology in topol.top, OUTPUTNAME_H.pdb, TIP3P water, charmm36-ccd2md forcefield, and charged termini (excepting starting CYST or GLYM which are set to None). Only applicable if NOT embedding the system in a membrane.', action='store_true')

sys_opts = parser.add_argument_group('membrane-embedded system options')

sys_opts.add_argument('-mem', '--membrane', help='Embed the converted system into a membrane. Note that this will lead to minor rearragements of any ligands. After this flag, it is possible to add the arguments for the upper and lower leaflet in the form "-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1" where the codes represent lipids and the numbers represent the ratio between them. The default is two leaflets of pure POPC. Note that membrane embedding may result in minor rearrangements of ligands.', action='store_true')
sys_opts.add_argument('-C',  '--conc',      help='Concentration of NaCl in system - charge balance is maintained, Default = 0.15.', default = 0.15)
sys_opts.add_argument('-mp',  '--mempro',   help='Additional arguments for embedding the protein in the membrane using MemPrO. Add any additional arguments after this flag - default 5 grid points and 15 minimisation operations.', action='store_true')
sys_opts.add_argument('-at', '--cg2at',     help='Additional arguments to pass to CG2AT. Add additional arguments after this flag - note that this may require interactivity.', action='store_true')
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

ligands, atoms, IDs, types, locs = FuncConv.get_residues(input_data, 'CCD', SMILES, args.ligchain)

# Reorder residues
# ----------------
	
output_data = input_data.to_dict('records')

user_CCD_converted = []

for i, lig in enumerate(ligands):
    if types[i] != 'CHARMM':
        output_residues = FuncConv.convert_atomistic(lig, atoms[i], 'CCD', args.ligchain)
        output_data[IDs[i][0]:IDs[i][-1]+1] = output_residues
    else:
        if len(np.where(np.array(user_CCD_converted)==lig)[0])==0:
            nm = lig[:-5] + ' (userCCD)'
            print('# INFO: Gathering database information for molecule name {}'.format(nm))
            user_CCD_converted.append(lig)

# Write output
# -------------

if args.outputfile.count('/') != 0:
    output_dir = '/'.join(args.outputfile.split('/')[:-1])
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

basename = '.'.join(args.outputfile.split('.')[:-1])
        
PDBfile = basename+'_nomem.pdb' if args.membrane else args.outputfile
FuncConv.write_PDB(PDBfile, output_data, title=title, cryst=cryst, ligand_chains=args.ligchain)

# =================== #
# Embed into membrane #
# =================== # 

if args.membrane:
    # Convert system to CG
    # --------------------

    prot        = True if len(np.concatenate([np.array(i) for i in IDs])) != len(input_data) else False    
    output_data = FuncConv.to_CG(args.inputfile, basename+'_CG.pdb', input_data, ['-elastic'], ligands, input_data,
                                 types, locs, prot)

    CG_output = basename+'_CG_system.pdb'

    FuncConv.write_PDB(CG_output, output_data, title=title, cryst=cryst, ligand_chains=args.ligchain)

    # Embed in membrane
    # -----------------
    
    FuncConv.build_membrane_CG(ligands, CG_output, args.outputfile, command_line, args.mempro, args.membrane, args.conc) 

    # Convert back to atomistic
    # -------------------------

    FuncConv.convert_membrane_at(output_data, basename, command_line, args.cg2at)
       
    subprocess.run(['scp', basename+'_CG2AT/FINAL/final_cg2at_aligned.pdb', args.outputfile])

if args.membrane or not args.pdb2gmx:
    FuncConv.get_topology_atomistic(args.outputfile, args.membrane, output_data=pd.DataFrame.from_dict(output_data))
else:
    FuncConv.get_topology_atomistic(args.outputfile, args.membrane, at_command=command_line)




