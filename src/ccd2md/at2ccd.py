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
    
    parser = argparse.ArgumentParser(description='Convert from CHARMM to CCD ordering. Please see keb721/CCD2MD on github for a list of the allowed CCD codes.', add_help=False)
    
    parser.add_argument('inputfile', help='Input file name - .cif, .pdb or .gro.')
    parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.')

    parser.add_argument('-T', '--title',     help='Optionally specify the title of the generated pdb file.', default = None)    
    parser.add_argument('-L', '--ligchain', help='Output ligands in their own chains - default is off.', action='store_true')
    

    info = parser.add_argument_group()
    info.add_argument('-V', '--Version', action='version', version='Version '+FuncConv.__version__)
    info.add_argument('-h', '--help',    action='help',    help='Show this message and exit.')
    info.add_argument('-CF', '--configuration', help='Specifiy configuration file. For an example input please see keb721/CCD2MD on GitHub.', default=None)
    
    args, unknownargs = parser.parse_known_args()
    command_line      = np.array(sys.argv)
    
    if args.configuration != None:    
        args, command_line = FuncConv.read_configuration_file(args.configuration, args, command_line)
    
    # Read input data
    # ---------------
    
    input_data, title, cryst = FuncConv.read_in(args.inputfile, args.title)
    
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
