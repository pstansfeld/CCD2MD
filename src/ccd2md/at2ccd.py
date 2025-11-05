#######################################
#                                     # 
#    Convert between CCD and CHARMM   # 
#                                     #
#######################################

# Import relevant functions
# -------------------------

from ccd2md import FuncConv
import argparse
import pandas as pd
from copy import deepcopy
import sys

def main():
        
    # Get command line arguments
    # ---------------------------
    
    parser = argparse.ArgumentParser(description='Convert from CHARMM to CCD ordering. Please see keb721/CCD2MD on github for a list of the allowed CCD codes.')
    
    parser.add_argument('inputfile', help='Input file name - .cif, .pdb or .gro.')
    parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.')
    
    parser.add_argument('-L', '--ligchain', help='Output ligands in their own chains - default is off.', action='store_true')
    
    info = parser.add_argument_group()
    info.add_argument('-V', '--Version', action='version', version='Version '+FuncConv.__version__)
    
    args, unknownargs = parser.parse_known_args()
    
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
    
    print('# INFO: Assuming that chains are labelled sequentially.')
    
    # Find residues to reorder
    # ------------------------
        
    ligands, atoms, IDs, types, locs = FuncConv.get_residues(input_data, 'CHARMM', [], args.ligchain)
    
    # Reorder residues
    # ----------------
            
    output_data = input_data.to_dict('records')
    
    for i, lig in enumerate(ligands):
        output_residues = FuncConv.convert_atomistic(lig, atoms[i], 'CHARMM', args.ligchain)
        output_data[IDs[i][0]:IDs[i][-1]+1] = output_residues 
    
        
    # Write output
    # -------------
    
    print('# INFO: Assuming that an entirely new PDB file is being written rather than modified in place')    
    FuncConv.write_PDB(args.outputfile, output_data, title=title, cryst=cryst, ligand_chains=args.ligchain)


if __name__ == "__main__":
    main()
