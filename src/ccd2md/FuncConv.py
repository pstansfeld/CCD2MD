##############################################################
#                                                            #
#    Functions to convert between CCD and CHARMM ordering    #
#                                                            #
##############################################################

'''
General file for functions to convert between CCD/CHARMM/Martini files

Last Update: K Blow 30/03/26

Contains:

Input file parsers
-------------------

read_configuration_file(conffile, args, command_line, flag_dict=flag_dict) # Read configuration file and cmpare to command line inputs

read_CIF(INNAME)   # Read CIF file INNAME, return system data (note, only data relevant to PDBs is retained)
read_PDB(INNAME)   # Read PDB file INNAME, return system data and optionally title and crystallographic information
read_GRO(INNAME)   # Read GRO file INNAME, return system data and optionally title and crystallographic information

get_title(title, max_pdb=80)   # Parse command line title and convert to PDB format

determine_element(NAME) # Determine element type from atom name

read_in(INNAME, INTITLE) # Get input dictionary

General conversion information
-------------------------------

get_residues(system_data, data_type, SMILES=[], ligand_chain=False, database=df, degenerate=degenerate_names) 
# Get information about the residues in system_data which need to be converted (data_type gives the current data type. Add SMILES keys and determine if ligands should be separated into chains.
which_degenerate(name, data, Ref_loc=Ref_data)                                            
# Choose between degenerate data for NAME
get_reference(NAME, database=df, Ref_loc=Ref_data)                            
# Get reference information for residue NAME. Prevent displaying information multiple times for the same moleucule, although this may need refinement.

get_command_line_parameters(command_line, flags) # Determine additional command line parameters to pass to various commands.


Atomistic conversion 
---------------------

convert_atomistic(name, inputdict, inputtype, ligand_chains=False, database=df, Ref_loc='Ref_data/') # Change ordering between atomistic orders - CCD (including SMILES strings where accepted) and CHARMM


Coarse-grained conversion 
--------------------------

get_CG_params(command_line, martinize, elastic, go)                                       # Determine additional flags to be passed to martinize2
prot_to_CG(inputfile, basename, input_data, martinizeparams, ligands, types, database=df) # Convert protein to CG using martinize2 on the command line. 
to_CG(inputfile, outputfile, input_data, martinizeparams, ligands, system_data, types, convresi, prot,
          database=df, Ref_loc=Ref_data, mass=masses, newlipidome=False                   # Take in an atomic list of locations, and return a CG representation.


Membrane embedding 
-------------------

build_membrane_CG(ligands, CG_output, outputfile, command_line, mempro_additional, memprod_additional, membrane_comp, ion_conc, database=df, newlipidome=False, num_CPUs=1) # Builds a CG membrane around the system using MemPrO
convert_membrane_at(system_data, basename, command_line, CG2AT_additional, CG_name, executable, newlipidome=False)
                                                      # Converts membrane using CG2AT-lite

Output file writing
--------------------

write_PDB_atom_line(OUTFILE, counter, data)                            # For a given OUTFILE, write data corresponding to a given atom
write_PDB(OUTFILE, ordered_dict, title=[], cryst=[], ligchain = False) # For a given OUTFILE, write all of the data for the system (ordered_dict)

check_residue_number(ordered_dict) # Renumber residues where this exceeds max PDB can handle
convert_vectors(box_vecs)          # Convert box information from GRO format to PDB format


Globular protein preparation
-----------------------------

make_glob(system, command_line, executable, simtype, editconf, solvate, ions, conc, topol='topol.top')
# Create box, solvate and add ions for non-membrane protein


MD generation
--------------

make_gmx(system, command_line, tpr_type, executable, ndx=None, run=False, topol='topol.top') # Write and optionally run tpr file for system


Steering
---------

steer_reference(ref_data, st_ligs, steerpath, executable) # Get reference positions for ligands
create_ref_files(ref_pos, lig, steerpath, executable)     # Create all necessary reference files
get_ligand(output_data, ligand, loc)                      # Get position of predicted ligands
steer_atomistic(predicted_file, ligand, executable)       # Steer reference file to prediction
steer_ligand(ligand, steer_files, executable)             # Perform steering for a ligand


Topology generation
--------------------

get_topology_CG(outputfile, membrane, ligands, prot, inputfile, newlipidome=False, topoldir=None, database=df) # Write topology for CG system
get_topology_atomistic(outputfile, membrane, executable, at_command=None, output_data=None)                    # Write topology for atomistic system


'''

__version__ = "1.1.0"

import numpy as np
import pandas as pd
import warnings
import subprocess, os, shutil, sys, platform
import re, glob
import ast

# Locations of packages
# ---------------------

import ccd2md

CCD2MD_dir = os.path.dirname(ccd2md.__file__)+'/'  # Location of CCD2MD
# CCD2MD_dir = os.path.dirname(os.path.abspath(__file__))+'/'  # Location of CCD2MD
Ref_data   = CCD2MD_dir + 'Ref_data/'
CHARMMPath = CCD2MD_dir + 'charmm36-ccd2md.ff/'
oldmartini = CCD2MD_dir + 'martini_v3.itp'
newmartini = CCD2MD_dir + 'martini_new_lipidome.itp'
oldinsane  = CCD2MD_dir + 'MemPrO/Insane4MemPrO.py'
newinsane  = CCD2MD_dir + 'MemPrO/Insane4MemPrO_new_lipidome.py'
mdpPath    = CCD2MD_dir + 'mdp_files/'

# General information
# -------------------

# Create an empty dictionary linking CIF labels to PDB output - note worls with CHAI and CCD/PDB output but due to
# nomenclature reduncancies not guaranteed to work in all cases

CIF_keywords = { # Keywords given by Chai
                '_atom_site.group_PDB'             : 'entry', # PDB entry type
                '_atom_site.label_atom_id'         : 'name',  # Unique name for the atom
                '_atom_site.label_comp_id'         : 'resnm', # Residue name
                '_atom_site.auth_asym_id'          : 'chain', # Chain identifier
                '_atom_site.auth_seq_id'           : 'resi',  # Residue number 
                '_atom_site.Cartn_x'               : 'x',     # x coordinate
                '_atom_site.Cartn_y'               : 'y',     # y coordinate
                '_atom_site.Cartn_z'               : 'z',     # z coordinate
                '_atom_site.occupancy'             : 'occ',   # Occupancy
                '_atom_site.B_iso_or_equiv'        : 'B',     # B factor
                '_atom_site.type_symbol'           : 'elem',  # Element symbol
                # Keywords from CCD - note some are missing, duplicates are the other possible name
                '_chem_comp_atom.pdbx_component_atom_id'   : 'name',  # Unique name for the atom
                '_chem_comp_atom.pdbx_component_comp_id'   : 'resnm', # Residue name
                '_chem_comp_atom.pdbx_model_Cartn_x_ideal' : 'x',     # x coordinate
                '_chem_comp_atom.pdbx_model_Cartn_y_ideal' : 'y',     # y coordinate
                '_chem_comp_atom.pdbx_model_Cartn_z_ideal' : 'z',     # z coordinate
                '_chem_comp_atom.type_symbol'              : 'elem',  # Element symbol
                # Keywords for userCCD cif file - note some are missing, duplicates are the other possible name
                '_chem_comp_atom.atom_id'   : 'name',  # Unique name for the atom
                '_chem_comp_atom.comp_id'   : 'resnm'  # Residue name
                }

PDB_keywords = {'entry' : [0, 6],
                'atomi' : [6, 12],
                'name'  : [12, 16],
                'resnm' : [17, 21],
                'chain' : [21, 22],
                'resi'  : [22, 26],
                'x'     : [30, 38],
                'y'     : [38, 46],
                'z'     : [46, 54],
                'occ'   : [54, 60],
                'B'     : [60, 66],
                'elem'  : [76, 78]}

GRO_keywords = {'name'  : [10, 15],
                'resnm' : [5, 10],
                'resi'  : [0, 5],
                'x'     : [20, 28],
                'y'     : [28, 36],
                'z'     : [36, 44]}

masses = {'C'  : 12.01100,
          'O'  : 15.99940,
          'N'  : 14.00700,
          'P'  : 30.97400,
          'S'  : 32.06000}

float_keys = ['x', 'y', 'z', 'occ', 'B']
chars      = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
              'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

# All flags read by parts of CCD2MD

flag_dict = {'CF'    : 'configuration', # NEED TO ADD TO FILES
             'T'     : 'title',         # NEED TO ADD TO FILES
             'S'     : 'SMILES',
             'L'     : 'ligchain',
             'gh'    : 'pdb2gmx',
             'nl'    : 'newlipidome',
             'M'     : 'martinize',
             'G'     : 'go',
             'E'     : 'elastic',
             'mem'   : 'membrane',
             'mp'    : 'mempro',
             'mdef'  : 'memprod',
             'C'     : 'conc',
             'ncpu'  : 'num_cpus',
             'at'    : 'cg2at',
             'CGEM'  : 'CG_energy_minimise',
             'rCGEM' : 'run_CG_energy_minimise',
             'AAEM'  : 'AA_energy_minimise',
             'rAAEM' : 'run_AA_energy_minimise',
             'CGeq'  : 'CG_equil',
             'rCGeq' : 'run_CG_equil',
             'AAeq'  : 'AA_equil',
             'rAAeq' : 'run_AA_equil',
             'CGprd' : 'CG_prod',
             'AAprd' : 'AA_prod',
             'ndx'   : 'make_ndx',
             'St'    : 'steer',
             'Stnm'  : 'steer_name',
             'Stref' : 'steer_ref',
             'g'     : 'gmx',
             'B'     : 'box',
             'SV'    : 'solvate',
             'I'     : 'ions'}    


short_flag = lambda x: '-' + x
long_flag  = lambda x: '--'  + flag_dict[x]

all_flags = [flag(key) for key in flag_dict.keys() for flag in (short_flag, long_flag)]


df = pd.read_csv(Ref_data+"database.csv", index_col='Name', skipinitialspace=True)

base_ptms     = ['CYST', 'CYSD', 'CYSP', 'CYSG', 'CYSF', 'GLYM']
PTMs          = set(base_ptms + [ptm + '_user' for ptm in base_ptms])
terminal_PTMs = ['CYST', 'GLYM']

degenerate_names = ['ATP'] # CCD code and CHARMM code the same 

# Global information
# ------------------

lig_names = []
converted = []

element_names = True
steer_warn    = True

# MD information
# ---------------

mdp = {'CGEM'       : 'em.mdp',
       'AAEM'       : 'em.mdp',
       'CGeq_prot'  : 'cg-equil.mdp',
       'CGeq_lip'   : 'cg-equil-membrane.mdp',
       'AAeq_prot'  : '10ns-pr.mdp',
       'AAeq_lip'   : '10ns-pr-membrane.mdp',
       'CGprd_prot' : 'cg-10us.mdp',
       'CGprd_lip'  : 'cg-10us-membrane.mdp',
       'AAprd_prot' : '500ns.mdp',
       'AAprd_lip'  : '500ns-membrane.mdp'}


# ======================== #
# ------------------------ #
#    INPUT FILE PARSERS    #
# ------------------------ #
# ======================== #


# Configuration file parser
# -------------------------

def read_configuration_file(configuration, args, command_line, flag_dict=flag_dict):
    ''' Read information from configuration file and command line to generate input data. '''

    # Read line, get the corresponding flag and check against args
    # Change args and command_line if needed

    # ======================================= #
    # Get information from configuration file #
    # ======================================= #
    
    # Read configuration file and strip newlines
    data = open(configuration, 'r').read().split('\n')
    data = [line for line in data if len(line)!=0 and line.strip()[0]!='#']    

    config_info = {}
    
    for line in data:
        # Split into flag and value
        l    = line.split('=')
        flag = l[0].strip()
        flag = flag if flag_dict.get(flag) == None else flag_dict[flag]
        val  = l[1].split('#')[0].strip() # Remove trailing comments

        if config_info.get(flag) != None:
            print('ERROR: Information for the {} flag has been set multiple times within the configuration file.'.format(flag))
            exit()

        if flag == 'num_cpus':
            val = int(val)
        elif flag == 'SMILES':
            val = val.split()
        elif flag == 'steer_name':
            val = val.split()
        elif val.lower() == 'true':
            val = True
        elif val.lower() == 'false':
            val = False
            
        config_info[flag] = val

    # Convert to long flags and check every argument for agreement 

    test = list(config_info.keys())
    test = np.array(test)

    # ======================= #
    # Compare to command line #
    # ======================= #

    # Generate putative values
    # ------------------------
    
    for key in args.__dict__:
        if key[0] == '_':
            continue
        if key == 'inputfile' or key == 'outputfile':
            # Positional argument
            # Can also have input/output 
            flag = [f for f in test if f.startswith(key[:5])]
            if len(flag) == 0:
                continue
            elif len(flag) != 1:
                print('ERROR: The {} has been set multiple times within the configuration file.'.format(key))
            else:
                arg_test = config_info.pop(flag[0])

        elif len(np.where(test == key)[0]) == 0:
            if key == 'elastic':
                # Test for secondary structure information
                # Note, not testing go as mutually exclusive groups

                if len(np.where(test == 'secondary structure')[0]) != 0:
                    secondary_structure =  config_info.pop('secondary structure')
                    if secondary_structure[0] != '-':
                        # Written as elastic or go
                        secondary_command = '--'+secondary_structure
                        if len(np.where(test == 'secondary structure additional')[0]) != 0:
                            secondary_command += config_info.pop('secondary structure additional')
                    command_line = np.append(command_line, secondary_command)
                    
                    secondary_info = secondary_structure.split()[0]
                    
                    if secondary_info == 'elastic' or secondary_info == '--elastic' or secondary_info == '-E':
                        args.elastic = True
                    elif secondary_info == 'go' or secondary_info == '--go' or secondary_info == '-G':
                        args.go = True
                    else:
                        print('ERROR: Unknown secondary structure method')
                        exit()
            
            continue
        else:
            # Optional argument
            arg_test = config_info.pop(key)

        # Update command line input if not set
        # -------------------------------------
        
        # Have the value from configuration file, now need to compare to the command_line
        # Note changes to e.g. new lipidome and concentration to take defaults into account
        
        if args.__dict__[key] == None:
            # No value set, and not a flag for additional command line information
            # Can replace without issue
            args.__dict__[key] = arg_test

        elif args.__dict__[key] == False:
            
            # Not set, but represents a command line presence flag 
            # Could be membrane composition which also has additional fields in configuration file

            # Also need to test for some MD parameters in case run is added and not the initial one

            if arg_test == False:
                continue
            
            args.__dict__[key] = True

            command_line = np.append(command_line, ['--{}'.format(key)])
            
            # May be membrane = true
            if arg_test != True:
                command_line = np.append(command_line, arg_test.split())

            if key == 'membrane':
                # Test for presence of composition flags
                leaflet = {'lower': '-l', 'upper': '-u'}

                # Get composition information
                for leaf in leaflet.keys(): 
                    if len(np.where(test == '{} leaflet composition'.format(leaf))[0]) != 0:
                        lflet = str(config_info.pop('{} leaflet composition'.format(leaf)))
                        lflet = lflet.split()
                        lflet = [entry.split(':') for entry in lflet]
                        lflet = [val for lipid in lflet for val in lipid]
                        if len(lflet)%2 != 0 and len(lflet) != 1:
                            print('ERROR: {} leaflet composition is in an incorrect format.'.format(leaf))
                            exit()
                        
                        composition = lflet if len(lflet)==1 else [lflet[2*i]+':'+lflet[2*i+1] for i in range(int(0.5*len(lflet)))]
                        composition = [entry for lipid in composition for entry in ('{}'.format(leaflet[leaf]), lipid)]
                        command_line = np.append(command_line, composition)
            
        else:
            # This has been set on the command line
            print('ERROR: The input flag --{} has been set on the command line and within the configuration file.'.format(key))
            exit()

        
    # Check for entries in configuration file which don't correspond to args
    
    if len(config_info) != 0:
        print('# WARNING: some parameters set in {} are not applicable to the programme chosen. Please ensure that you have specified the correct conversion.')
        print('# WARNING: these parameters are {}'.format(', '.join(config_info.keys())))

        
    return args, command_line



# Read CIF files
# ---------------

def read_CIF(name):
    ''' Read CIF file into molecule dictionary. '''

    fle = open(name, 'r').read()
    print('# INFO: Reading CIF file {}. Note that this assumes all fields are consistent in CIF ordering'.format(name))

    title      = ['TITLE     Converted from {} using CCD2MD.'.format(name)]
    CCD_direct = False
    resimiss   = False
    
    models  = fle.split('\ndata_')
    mod_out = [] 
    
    # Check for comments only
    models = [model for model in models if '\nloop' in model]
    
    
    for model in models:
        CIF = model.split('loop_')
        
        out_info = [{}]
    
        for block in CIF:
            # Identify what is in each block of model information (separated by data labels)
            data      = block.split('\n')
            data_labs = np.array([i for i in data if (len(i)!=0 and i[0]=='_')])
            is_key    = [CIF_keywords[i.strip()] if (CIF_keywords.get(i.strip()) != None) else '_'  for i in data_labs]

            if (is_key.count('_') ==  len(is_key)):
                continue

            if not CCD_direct:
                keys      = [i.strip() if (CIF_keywords.get(i.strip()) != None) else '_' for i in data_labs]
                for k in keys:
                    if k[:5]=='_chem':
                        CCD_direct = True
                        print('# WARNING: This file may have possible alternative labels and missing data, which could cause issues.')
                        break
                
            data = np.array([i for i in data if (len(i)!=0 and i[0]!='_' and i[0]!='#')])
        
            if len(out_info) == 1:
                out_info = [{} for i in range(len(data))]
            elif len(out_info)!=len(data):
                if list(out_info[0].keys()) == ['elem']:
                    # Only element saved - for some CCD files there is an entry per residue
                    out_info = [{} for i in range(len(data))]
                else:
                    print('ERROR: Number of atoms inconsistent in CIF file.')
                    sys.exit()

            line_start  = True
            delete_dict = []
        
            for i, line in enumerate(data):
                # Read in each line of data 
                line = line.split()
                # Assuming consistent ordering between blocks
                if len(line) != len(is_key):
                    # Chai now appears to be implementing a maximum line length. Test for this
                    if line_start:
                        line.extend(data[i+1].split())
                        line_start = False
                    else:
                        line_start = True
                        delete_dict.append(i)
                        continue
                
                for e, key in enumerate(is_key):
                    if key!='_':
                        if key == 'resi':
                            out_info[i][key] = int(line[e])
                        elif key in float_keys:
                            out_info[i][key] = float(line[e])
                        elif key == 'name':
                            # May have a dash, important but difficult to deal with!
                            out_info[i][key] = line[e][1:-1] if line[e][0] == '"' else line[e]
                        else:
                            out_info[i][key] = line[e]

            for i in range(len(delete_dict)-1, -1, -1):
                _ = out_info.pop(delete_dict[i])
        
        # Deal with missing data in CCD files
        if CCD_direct:
            # Check for the presence/absence of information
            missing = {'occ': 1.0, 'B': 0.0, 'entry': 'ATOM', 'chain': 'A'}
            for key in missing.keys():
                if list(out_info[0].keys()).count(key) == 0:
                    out_info = [dict(atom, **{key: missing[key]}) for atom in out_info]

            if list(out_info[0].keys()).count('resi') == 0:
                if not resimiss:
                    print('# WARNING: Inferring residue IDs - these will be incorrect if there are two sequential residues of the same type')
                    resimiss = True
                resi               = 1
                out_info[0]['resi'] = resi
                for i in range(1, len(out_info)):
                    resi = resi+1 if out_info[i]['resnm'] != out_info[i-1]['resnm'] else resi
                    out_info[i]['resi'] = resi

        else:
            print('# INFO: Co-folding predictions usually have accuracy encoded in B factor.')

        mod_out.append(out_info)

    return mod_out, title



# Read PDB files
# ---------------

def read_PDB(name):
    ''' Read in relevant information from PDB file. '''

    print('# INFO: Reading PDB file {}.'.format(name))
    
    PDBfile = open(name)
    PDB     = PDBfile.read().split('\nEND')[0]
    PDB    += '\n'

    PDB     = PDB.split('\n')[:-1]

    out_info = []

    title = []
    cryst = []
    
    for line in PDB:
        # Only look for certain molecules
        if line.count('TITLE') == 1:
            title.append(line)
        elif line.count('CRYST') == 1:
            cryst.append(line)
        elif line.count('REMARK') == 1:
            title.append(line)
        elif line.count('HETATM') == 0 and line.count('ATOM') == 0:
            continue
        else:
            # Append relevant information
            out_info.append({})
            for i, key in enumerate(PDB_keywords.keys()):
                sect = line[PDB_keywords[key][0]:PDB_keywords[key][1]]
                if key == 'resi':
                    out_info[-1][key] = int(sect.strip(' '))
                elif key in float_keys:
                    out_info[-1][key] = float(sect.strip(' '))
                elif key == 'chain':
                    out_info[-1][key] = sect.strip(' ') if len(sect.strip(' '))!= 0 else ''
                else:
                    out_info[-1][key] = sect.strip(' ')
                    
                    
    if len(title) == 0:
        title = ['TITLE     Converted from {} using CCD2MD.'.format(name)]
                    
    return out_info, title, cryst



# Read gro files
# ---------------

def read_GRO(name):
    ''' Read in relevant information from GRO file. '''

    print('# INFO: Reading GRO file {}.'.format(name))
    print('# WARNING: Inferring chain information from residue IDs.')
    
    GROfile = open(name)
    GRO     = GROfile.read().split('\n')[:-1]

    out_info = []

    title = [GRO[0]]
    atoms = int(GRO[1])

    prev_id = 0 
    char_id = 0
    
    for i in range(atoms):
        # Get information for each line

        line = GRO[2+i]
        
        # Append missing information
        out_info.append({'entry': 'ATOM', 'B': 0, 'occ': 1, 'elem': ''})
        
        for i, key in enumerate(GRO_keywords.keys()):
            sect = line[GRO_keywords[key][0]:GRO_keywords[key][1]]
            if key == 'resi':
                out_info[-1][key] = int(sect.strip(' '))
            elif key in float_keys:
                out_info[-1][key] = float(sect.strip(' ')) * 10  # nm -> A
            else:
                out_info[-1][key] = sect.strip(' ')

        if out_info[-1]['resi'] < prev_id:
            char_id += 1

        out_info[-1]['chain'] = chars[char_id]
        prev_id = out_info[-1]['resi']

    cryst  = convert_vectors(GRO[atoms+2])
        
    return out_info, title, cryst



# Convert title to  PDB format
# -----------------------------

def get_title(command_line_title, max_pdb=80):
    ''' Convert a title given via the command line to PDB format.'''
    
    # Check for PDB format
    if command_line_title[:8] == 'TITLE   ' and len(command_line_title) <= max_pdb:
        title = [command_line_title]
    else:
        title = []
        # Either missing starting title or too long
        if len(command_line_title) <= (max_pdb - 10):
            title = ['TITLE     ' + command_line_title]
        else:
            if command_line_title[:8] != 'TITLE   ':
                command_line_title = 'TITLE     ' + command_line_title
            
            curr_line   = []
            split_title = command_line_title.split(' ')
            curr_line.append(split_title.pop(0))
            len_curr    = len(curr_line[-1])
            num_line    = 1

            while len(split_title)!=0:
                test = len_curr + 1 + len(split_title[0])
                if test > max_pdb:
                    # Append to title, create new line
                    title.append(' '.join(curr_line))
                    num_line += 1
                        
                    curr_line = 'TITLE   {:>2}'.format(num_line)
                                                            
                    curr_line = curr_line.split(' ')
                    curr_line.append(split_title.pop(0))
                    len_curr  = 11 + len(curr_line[-1])
                    
                else:
                    # Add to current line
                    len_curr = test
                    curr_line.append(split_title.pop(0))
                
            # Add to title array when finished parsing string
            title.append(' '.join(curr_line))        
    
    return title        
                        


# Get element types
# ------------------

def determine_element(name):
    ''' Determine element type from atom name - note this may be problematic. '''

    # Prioritise element types - assume all H names will be hydrogen
    if name.count('H')==1:
        elem = 'H'
    elif name.count('O')!=0:
        elem = 'O'
    elif name.count('N')!=0:
        elem = 'N'
    elif name.count('C')!=0:
        elem = 'C'
    elif name.count('P')!=0:
        elem = 'P'
    elif name.count('S')!=0:
        elem = 'S'
    else:
        elem = ''
             
    return elem



# Get input dictionary
# ---------------------

def read_in(inputfile, in_title):
    ''' Generate input dictionary from an input file. '''
    
    if inputfile[-3:] == 'cif':
        tmp, title = read_CIF(inputfile)
        if len(tmp)!=1:
            print('ERROR: Multiple models in input structure.')
            sys.exit()
        tmp   = tmp[0]
        cryst = []
    elif inputfile[-3:] == 'gro':
        tmp, title, cryst = read_GRO(inputfile)
    else:
        print('pdb')
        if inputfile[-3:] != 'pdb':
            print('# WARNING: Assuming that input file {} is written in PDB style despite file extension'.format(inputfile))
        tmp, title, cryst = read_PDB(inputfile)
        
    print('# INFO: Any Hs present will be removed')

    if in_title != None:
        title = get_title(in_title)
        
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
            atom['elem'] = determine_element(atom['name'])
            if atom['elem'] != 'H':
                input_data.append(atom)        
    
    input_data = pd.DataFrame.from_dict(input_data, orient='columns')
    
    return input_data, title, cryst
    

   
# ==================================== #
# ------------------------------------ #
#    GENERAL CONVERSION INFORMATION    #
# ------------------------------------ #
# ==================================== #


# Determine residues to convert
# ------------------------------

def get_residues(system_data, data_type, SMILES=[], ligand_chain=False, database=df, degenerate=degenerate_names):
    ''' Get the names of the residues for which there is an associated database entry, and the atoms 
        corresponding to these residues. Where the input is SMILES, change the residue name from LIG. '''
    
    # Get list of residue names
    # --------------------------
    
    chains = (set(system_data['chain']))
    
    full_to_convert = []
    full_convert_ID = []
    convres         = []
    convert         = []
    convresi        = []
    types           = []
    lig_IDs         = []
    
    prev_max   = 0
    prev_chain = ''
    max_chnID  = len(set(system_data['chain'])) # 1 offset applied due to 0 indexing

    first_smiles = True
    first_skip   = True

    sort_chn = sorted(chains)
    if sort_chn[0] == '':
        # Put empty chains last
        sort_chn.append(sort_chn.pop(0))
    
    for j, chain in enumerate(sort_chn):
        try:
            chain_data    = system_data.loc[system_data['chain']==chain]
            first_residue = min(chain_data['resi'])
            last_residue  = max(chain_data['resi'])
        except ValueError:
            chain_data    = system_data.loc[system_data['chain']=='']
            first_residue = min(chain_data['resi'])
            last_residue  = max(chain_data['resi'])

        residues  = []
        get_kind  = []
        
        num_residues = last_residue+1-first_residue
        
        for i in range(first_residue, last_residue+1):
            try:
                residues.append(list(chain_data.loc[chain_data['resi']==i, 'resnm'])[0])
                if residues[-1] in degenerate:
                    # Store dengerate residue information to determine if user code, or CCD code
                    get_kind.append(which_degenerate(residues[-1], chain_data.loc[chain_data['resi']==i, 'name']))
                    
            except IndexError:
                # Allow for skips, but assume rare
                if first_skip:
                    print('# WARNING: Missing residues may cause issues')
                    first_skip = False
                
        # Find list of those which are in the database for reordering
        # ------------------------------------------------------------
        N_smiles = 0
        if len(SMILES)!=0:
            if first_smiles:
                print('# WARNING: Ligands represented as SMILES strings must have the same ordering in input file and as specified through the `-S`/`--SMILES` flag.')
                print('# WARNING: Attempting to internally assess and convert ligand naming - this may cause issues.')
                new_suffix = re.compile(r'\d_1') # Check for \d_1 to determine if new r old Chai naming -
                                                # Prevents error with C_1 etc.
                first_smiles = False

            # SMILES are their own chain - test residue names for different co-folding programmes
            smiles_loc = np.where(np.array([r[:3] for r in residues])=='LIG')[0]
            if len(smiles_loc) == 0:
                # Test different name
                smiles_name = re.compile(r'l\d\d')
                smiles_loc = np.where([smiles_name.search(r)!=None for r in residues])[0]
                    
            for res in smiles_loc:
                # This chain corresponds to a ligand taken from a SMILES string
                N_smiles += 1; convres.append(residues[res])
                global lig_names
                lig_names.append(residues[res])
                
                new_res       = SMILES.pop(0)
                res_data      = chain_data.loc[chain_data['resnm']==residues[res]]
                residues[res] = new_res

                for index, atom_data in res_data.iterrows():
                    curr          = system_data.iloc[index].to_dict()
                    curr['resnm'] = new_res
                    # Determine if renaming is necessary - default naming is e.g. C1
                    if '_' in curr['name']:
                        # AF3 naming omits the '_' and appends LIG with the chain
                        # Either C1_1 or C_1
                        if new_suffix.search(curr['name'])!=None:  
                            # Rename SMILES atoms to account for new Chai naming system - assuming all
                            # names have '_1' and are missing intial '_' => strip final '_2'
                            curr['name'] = curr['name'][:-2]
                        else:
                            # Old Chai name e.g. C_1
                            curr['name'] = ''.join(curr['name'].split('_'))
                    system_data.iloc[index] = curr

        cres  = [res for res in residues if len(database[database[data_type+'Name']==res])!=0]
        cert  = [database[database[data_type+'Name']==res].index[0] for res in cres]
        
        cert = [res if res not in degenerate else get_kind.pop() for res in cert]
        
                
        for residue in cert:
            # Rename the user CCD codes to function as CHARMM if running through e.g. ccd2at
            if residue[-5:]=='_user':
                types.append('CHARMM')
            elif residue[-7:]=='_SMILES':
                types.append('SMILES')
            else:
                types.append(data_type)
        
        convres.extend(cres)
        convert.extend(cert)
        
        if len(cres)==0: 
            # Protein only

            prev_max   = last_residue
            prev_chain = chain

        elif len([res for res in cres if res not in PTMs])==0: 
            # Modified protein only

            prev_max   = last_residue
            prev_chain = chain
            
            convresi.extend([[chain, i+first_residue] for i in range(len(residues)) if (residues[i] in cres)])
                        
        elif (not (ligand_chain)) and (len(cres) == num_residues):
            # Chain of just ligands
            # No PTMs possible            
            # Ligands should not have their own chain, here add to previous protein

            # For multi-protein multi-ligand need to ensure that there is no overlap of ligand ID for empy chain
            
            for i in range(num_residues):
                res_data = chain_data.loc[chain_data['resi']==i+first_residue]
                curr_resi = prev_max + i
                for offset in range(100):
                    curr_resi += 1
                    if len(np.where(np.array(lig_IDs) == curr_resi)[0])==0:
                        break
                    
                for index, atom_data in res_data.iterrows():
                    curr          = system_data.iloc[index].to_dict()
                    curr['chain'] = ''
                    curr['resi']  = curr_resi

                    system_data.iloc[index] = curr
                    
                convresi.extend([['', curr_resi]])
                lig_IDs.append(curr_resi)
            prev_max = prev_max + num_residues

        elif (len(cres) != num_residues):
            # Mixed chain
            # Ligands should have their own chain, here need to seperate from previous protein 
            # Also need to check for PTMs

            resids = []
            for res in set(cres):
                resis = set(chain_data['resi'].loc[chain_data['resnm']==res])
                resids.extend(resis)
            resids = sorted(resids)
            for resi in resids:
                if resi in PTMs:
                    convresi.extend([[chain, resi]])
                else:
                    curr_chain = chars[max_chnID] if ligand_chain else ''
                    # Assume that all in correct order => keep resi if not ligand chain
                    curr_resi  = 1                if ligand_chain else resi
                    res_data = chain_data.loc[chain_data['resi']==resi]

                    for offset in range(100):
                        if len(np.where(np.array(lig_IDs) == curr_resi)[0])==0:
                            break
                        curr_resi += 1

                    for index, atom_data in res_data.iterrows():
                        curr          = system_data.iloc[index].to_dict()
                        curr['chain'] = curr_chain
                        curr['resi']  = curr_resi

                        system_data.iloc[index] = curr
                    convresi.extend([[curr_chain, curr_resi]])
                    if ligand_chain:
                        max_chnID += 1
                    lig_IDs.append(curr_resi)
            prev_max   = 1                  if ligand_chain else curr_resi
            prev_chain = chars[max_chnID-1] if ligand_chain else chain
            
        elif (chain != 'A' and prev_chain == ''):
            # Ligands first?
            
            for i in range(num_residues):
                res_data  = chain_data.loc[chain_data['resi']==i+first_residue]
                new_chain = chars[max_chnID] if ligand_chain else ''
                new_resi  = 1                if ligand_chain else prev_max + i + first_residue + 1

                for offset in range(100):
                    if len(np.where(np.array(lig_IDs) == new_resi)[0])==0:
                        break
                    new_resi += 1

                    
                for index, atom_data in res_data.iterrows():
                    curr          = system_data.iloc[index].to_dict()
                    curr['chain'] = new_chain
                    curr['resi']  = new_resi
                
                    system_data.iloc[index] = curr
                convresi.extend([[new_chain, int(i)+1+prev_max] for i in range(len(residues)) if (residues[i] in cres)])
                if ligand_chain:
                    max_chnID += 1
                lig_IDs.append(new_resi)
                    
            prev_max   = 1                if ligand_chain else new_resi
            prev_chain = chars[max_chnID] if ligand_chain else ''
            
        else:
            # Ligands have their own chain, keep this

            prev_max   = num_residues
            prev_chain = chain
            
            resids     = [i+first_residue for i in range(num_residues)]
            convresi.extend([[chain, resid] for resid in resids])

    # For each occurence of a resiude to reorder, get the data and locations within system_data
    # -----------------------------------------------------------------------------------------
        
    for [chain, ID] in convresi:
        full_to_convert.append(system_data.loc[(system_data['chain']==chain) & (system_data['resi']==ID)])
        full_convert_ID.append(system_data.index[(system_data['chain']==chain) & (system_data['resi']==ID)].tolist())

    return convert, full_to_convert, full_convert_ID, types, convresi



# Choose between degenerate names
# --------------------------------

def which_degenerate(name, data, Ref_loc=Ref_data):
    ''' Determine if data with degenerate name corresponds to user CCD code or standard CCD code (i.e. CCD or CHARMM).'''

    CHARMM_order = np.genfromtxt(Ref_loc+name+'_CHARMM.txt', dtype=str)
    CCD_order    = np.genfromtxt(Ref_loc+name+'_CCD.txt', dtype=str)
    
    yes_CCD    = list(data)==list(CCD_order)
    yes_CHARMM = list(data)==list(CHARMM_order)

    assert (yes_CCD or yes_CHARMM), ' ERROR: does not match CCD or CHARMM'
    
    if yes_CCD and not yes_CHARMM:
        return name
    elif yes_CHARMM and not yes_CCD:
        return name+'_user'

    

# Extract database information
# ----------------------------

def get_reference(name, database=df, Ref_loc=Ref_data):
    ''' Get reference information from database for the relevant molecule. '''

    name = name.strip()  # Strip whitespace once

    with warnings.catch_warnings():
        # Suppress warning about elementwise comparison
        warnings.filterwarnings('ignore', 'elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison')
        global converted

        if len(np.where(np.array(converted)==name)[0])==0:
            if name[-5:] == '_user':
                nm = name[:-5] + ' (userCCD)'
            else:
                nm = name
            print('# INFO: Gathering database information for molecule name {}'.format(nm))
            converted.append(name)
        
    CHARMM_order = np.genfromtxt(Ref_loc+database.at[name.strip(), 'CHARMMName']+'_CHARMM.txt', dtype=str)
    CCD_order    = np.genfromtxt(Ref_loc+database.at[name.strip(), 'CCDName']+'_CCD.txt', dtype=str)
    with warnings.catch_warnings():
        # Suppress warning when CHARMM and CCD have the same atom names
        warnings.filterwarnings('ignore', 'genfromtxt: Empty input file: "{}"'.format(Ref_loc+name+'.txt'))
        rename = np.genfromtxt(Ref_loc + name.strip() +'.txt', dtype=str) # CHARMM name -> CCD name where different
    rename = rename.reshape(-1, 2)
    
    return CHARMM_order, CCD_order, rename



# Get additional command line parameters
# ---------------------------------------

def get_command_line_parameters(command_line, flag):
    ''' Determine additional command line parameters to pass to various commands. '''
        
    flags = ['-'+flag, '--'+flag_dict[flag]]
    
    CCD2MD_Flags = [f for f in all_flags if f not in flags]
    
    # User options - need to test, otherwise write defaults
    init = np.where(command_line==flags[0])[0]
    init = init[0] if len(init) != 0 else np.where(command_line==flags[1])[0][0]
    
    possargs = command_line[init+1:]

    # Possargs now contains possible arguments for the command.
    # Need to strip out any possibilities from different CCD2MD flags

    for f in CCD2MD_Flags:
        if len(np.where(possargs==f)[0])!=0:
            # This flag is present, remove this
            possargs = possargs[:np.where(possargs==f)[0][0]]
            
    # Should now only have the relevant parameters
    build_sys = [str(arg) for arg in possargs]
    
    return np.array(build_sys)



# ========================== #
# -------------------------- #
#    ATOMISTIC CONVERSION    #
# -------------------------- #
# ========================== #


def convert_atomistic(name, inputdict, inputtype, ligand_chains=False, database=df, Ref_loc=Ref_data):
    ''' Convert between atomstic representations (CCD and CHARMM). '''
    
    CHARMM_order, CCD_order, rename = get_reference(name, database, Ref_loc)
    
    output_order = []

    if inputtype=='CCD':
        outputorder = CHARMM_order
        inputnames  = rename[:, 1]
        outputnames = rename[:, 0]
        outputtype  = 'CHARMMName'
    else:
        outputorder = CCD_order
        inputnames  = rename[:, 0]
        outputnames = rename[:, 1]
        outputtype  = 'CCDName'
        
    for i, atomname in enumerate(outputorder):
        currname = atomname if len(np.where(outputnames == atomname)[0])==0 else \
                   inputnames[np.where(outputnames == atomname)[0][0]]

        output_order.extend(inputdict.loc[inputdict['name'] == currname].to_dict('records'))

        output_order[-1]['resnm'] = database.at[name, outputtype]
        output_order[-1]['name']  = atomname
        if not ligand_chains:
        # Preserve original chain ID from inputdict
            original_chain = inputdict.loc[inputdict['name'] == currname, 'chain'].values
            if len(original_chain) > 0:
                output_order[-1]['chain'] = original_chain[0]
            else:
                output_order[-1]['chain'] = ''


    return output_order



# ============================== #
# ------------------------------ #
#    COARSE-GRAINED CONVERSION   #
# ------------------------------ #
# ============================== #


# Martinize2 parameters
# ---------------------

def get_CG_params(command_line, martinize, elastic, go):
    ''' Determine additional command line parameters to pass to martinize2. '''

    martini = []
    # Note telling martinize2 to ignore Hs due to possible presence in input pdb file

    order = {'mart' : ['-M', '--martinize', martinize, 0],
             'elas' : ['-E', '--elastic',   elastic,   0],
             'go'   : ['-G', '--go',        go,        0]}
    
    for flag in order.keys():
        if order[flag][2]:
            # This flag is present, check for generic parameters
            # Command line options added after this flag                                        
            init = np.where(command_line==order[flag][0])[0]
            init = init[0] if len(init) != 0 else np.where(command_line==order[flag][1])[0][0]

            possargs = command_line[init+1:]

            # Possargs now contains possible arguments for the martinize2 command.
            # Need to strip out any possibilities from different flags

            for f in all_flags:
                if len(np.where(possargs==f)[0])!=0:
                    # This flag is present, remove this
                    possargs = possargs[:np.where(possargs==f)[0][0]]

            # Should now only have the relevant parameters
            order[flag][3] = len(possargs)
            if flag == 'elas':
                if len(np.where(possargs=='-elastic')[0])==0:
                    martini.append('-elastic')            
            elif flag == 'go':
                if len(possargs)==0:
                    print('# WARNING: A go network was specified but no additional commands will be passed to martinize2 - to generate a go network please include the relevant commands.')
                    martini.append('-go')
                if len(np.where(possargs=='-go')[0])==0:
                    martini.append('-go')            

            martini.extend(possargs)

    if (not order['elas'][2]) and (not order['go'][2]):
        print('# INFO: No network information was provided. Defaulting to elastic network with parameters as used in MemProtMD_Insane on pstansfeld GitHub.')
        martini.extend(['-elastic', '-ef', '500', '-eu', '1.0', '-el', '0.5', '-ea', '0', '-ep', '0'])


    SetFlags = {'-f'      : [1, 'input file'], '-x' : [1, 'output file'], '-o' : [1, 'output file'],
                '-ignore' : [1, 'additional ligands']}
    
    martini = np.array(martini)
    
    # Remove any pre-set flags
    for flag in SetFlags.keys():
        if len(np.where(martini==flag)[0])!=0:
            print('# WARNING: You have tried to overwrite the {} passed to Martinize2. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))
            martini = np.hstack((martini[:np.where(martini == flag)[0][0]], martini[np.where(martini == flag)[0][0]+1+SetFlags[flag][0]:]))

    if not martinize:
        # Add defaults
        martini = np.append(martini, ['-dssp'])
            
    return martini



# Convert protein
# ----------------

def prot_to_CG(inputfile, basename, input_data, martinizeparams, ligands, types, database=df):
    ''' Convert protein to CG using martinize2 on the command line. '''

    martfile = None
    
    if inputfile[-3:] == 'pdb' or inputfile[-3:] == 'gro':
        # Can directly utilise in martinize2
        martfile = inputfile
    else:
        # Cif file - need to test vermouth version and presence of PyCifRW
        version = subprocess.check_output(['martinize2', '-V'], universal_newlines=True)
        version = version.split()[-1].split('.')
        if int(version[0]) > 0 or int(version[1]) >= 14:
            # Martinize2 can handle vermouth, test PyCifRW presence
            try:
                import CifFile
                martfile = inputfile
            except ModuleNotFoundError:
                martfile = None
                
    if martfile == None:
        # Need to write an intermediate pdb file for use of martinize2
        # SMILES have been renamed in generating this
        print('# INFO: Writing intermediate pdb file for martinize2')
        martfile = basename + '_convert.pdb'
        prot_data = input_data.to_dict('records')
        write_PDB(martfile, prot_data)
        ligs = [database.at[lig, types[i]+'Name'] if types[i] != 'SMILES'
                else database.at[lig, 'CCDName'] for i, lig in enumerate(ligands)]
    else:
        # Initial input - need to consider LIG/LIG1 etc in residue names
        ligs = [database.at[lig, types[i]+'Name'] if types[i] != 'SMILES'
                else 'LIG' for i, lig in enumerate(ligands)]
        global lig_names
        ligs.extend(lig_names)
        ligs = list(set(ligs))

    outputpdb = basename + '_proteinCG.pdb'
    outputtop = basename + '_proteinCG.top'

    martini = [
        'martinize2',
        '-f', martfile,
        '-o', outputtop,
        '-x', outputpdb,
        '-ignore', ','.join(ligs),
        '-ignh'
    ]


    # # Check for mapping directory
    # if len(np.where(martinizeparams=='-map-dir')[0])!=0:
    #     print('# WARNING: You have tried to overwrite the mapping directory passed to Martinize2. This may cause an error.')
    # else:
    #     martini.extend('-map-dir', CCD2MD_dir+'martini')
        
    
    martini.extend(martinizeparams)    
    result = subprocess.run(martini)

    # Test output

    assert result.returncode==0, "ERROR: Failed to run martinize2, please check for errors in your input"
    
    return None



# Return full CG representation
# ------------------------------

def to_CG(inputfile, outputfile, input_data, martinizeparams, ligands, system_data, types, convresi, prot,
          database=df, Ref_loc=Ref_data, mass=masses, newlipidome=False):
    ''' Take in an atomic list of locations, and return a CG representation. '''
    
    filtered_ligands = [lig for lig in ligands if lig not in PTMs]
    basename         = '.'.join(outputfile.split('.')[:-1])
    
    output_dict = []
    if outputfile.count('/') != 0:
        output_dir = '/'.join(outputfile.split('/')[:-1])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    if prot:
        PTM_check = False
        if len(filtered_ligands) != len(ligands):
            # PTMs present - may be issues with merged chains or renumbering input chains
            test_params = np.array(martinizeparams)
            renumber = np.where(test_params == '-resid')[0]
            if len(renumber) != 0:
                if test_params[renumber + 1] == 'mol':
                    print('# WARNING: Martinize2 will renumber residues from 1. This may cause PTMs to be written with a 3 letter code instead of the correct 4 letter code.')
                    PTM_check = True # Double check properties before converting PTM
            else:
                martinizeparams = list(martinizeparams)
                martinizeparams.extend(['-resid', 'input'])
                
        print('# INFO: Running martinize2 on protein.')       
        prot_to_CG(inputfile, basename, input_data, martinizeparams, filtered_ligands, types, df)
        output_dict, _, _ = read_PDB(basename+'_proteinCG.pdb')

    complete = []

    print('# INFO: Bead positions are being calculated without H data.')

    lipidome = ' (updated martini 3 lipidome)' if newlipidome else ''
    
    for i, lig in enumerate(ligands):
        if lig in PTMs:
            # Need to rename, but not convert (handled by martinize2)

            redo = pd.DataFrame.from_dict(output_dict, orient='columns')
            new_name = lig if lig in base_ptms else lig[:-5] # Trim _user if present

            if PTM_check:
                # Tests to perform
                # 1. Check residue - assume that first 3 letters are unmodified AA
                # 2. Check number of beads and their names
                # If it passes, convert but print warning
                # If it fails, do not convert and print warning

                # Test residue name
                # -----------------
                
                curr_AA = list(redo.loc[(redo['chain'] == convresi[i][0]) & (redo['resi'] == convresi[i][1]), 'resnm']) # Should all be the same

                if len(curr_AA) == 0:
                    # No matching chain/residue combination
                    print('# WARNING: Residue {} on chain {} does not exist. The PTM {} will be present somewhere with the 3 letter code {}.'.format(convresi[i][1], convresi[i][0], new_name, new_name[:3]))
                    continue
                
                if curr_AA[0] != new_name[:3]:
                    print('# WARNING: Residue {} on chain {} is {} which does not correspond to the expected amino acid type for {}. This PTM will be present somewhere with the 3 letter code {}.'.format(convresi[i][1], convresi[i][0],  curr_AA[0], new_name, new_name[:3]))
                    continue

                # Test number of beads
                CG = eval(open(Ref_loc + database.at[lig, 'CGName'] + '_CG.txt').read())
                if len(CG) != len(curr_AA):
                    print('# WARNING: Residue {} on chain {} is does not have the expected number of beads for {}. This PTM will be present somewhere with the 3 letter code {} and {} beads.'.format(convresi[i][1], convresi[i][0], new_name, new_name[:3], len(CG)))
                    continue

                # Test bead names
                curr_AA = list(redo.loc[(redo['chain'] == convresi[i][0]) & (redo['resi'] == convresi[i][1]), 'name'])
                mismatch = [name for i, name in curr_AA if CG[i] != name]
                if len(mismatch) == 0 :
                    print('# WARNING: Residue {} on chain {} has the correct AA type and bead names as PTM {}. This will be renamed to reflect this, but this may be incorrect.'.format(convresi[i][1], convresi[i][0], new_name))
                    redo.loc[(redo['chain'] == convresi[i][0]) & (redo['resi'] == convresi[i][1]), 'resnm'] = new_name
                    output_dict = redo.to_dict('records')
                else:
                    print('# WARNING: Residue {} on chain {} does not have the correct bead names for as PTM {}. This PTM will be present somewhere with the 3 letter code {} and {} beads - check the bead types.'.format(convresi[i][1], convresi[i][0], new_name, new_name[:3], len(CG)))
            else:
                redo.loc[(redo['chain'] == convresi[i][0]) & (redo['resi'] == convresi[i][1]), 'resnm'] = new_name
                output_dict = redo.to_dict('records')
                
            continue
        
        if len(complete) == 0 or len(np.where(np.array(complete) == lig)[0]) == 0:
            if len(np.where(np.array(complete)==lig)[0])==0:
                if lig[-5:] == '_user':
                    nm = lig[:-5] + ' (userCCD)'
                else:
                    nm = lig
                print('# INFO: Converting {} from {} atomistic representation to coarse-grained{}.'.format(nm, types[i], lipidome))
                complete.append(lig)

        resi = convresi[i][1]
        residue = system_data.loc[(system_data['chain'] == convresi[i][0]) & (system_data['resi'] == convresi[i][1])]

        # Input file is dictionary
        if newlipidome:
            try:
                CG = eval(open(Ref_loc + 'newlipidome_'+database.at[lig, 'CGName'] + '_CG.txt').read())
            except FileNotFoundError:
                print('ERROR: A new lipidome mapping is not available for some of the ligands in your input script.')
                sys.exit()
        else:
            CG = eval(open(Ref_loc + database.at[lig, 'CGName'] + '_CG.txt').read())
        for bead_name in CG.keys():
            bead_dict = {
                'entry': 'ATOM',
                'resnm': database.at[lig, 'CGName'],
                'resi': resi,
                'chain': '',
                'occ': 1.0,
                'elem': ''
            }
            bead_dict['name'] = bead_name

            B = 0
            pos = []
            weights = []

            for atom_name in CG[bead_name][types[i]]:
                atom = residue.loc[residue['name'] == atom_name]
                if not atom.empty:
                    B += atom['B'].iloc[0]
                    pos.append([atom['x'].iloc[0], atom['y'].iloc[0], atom['z'].iloc[0]])
                    weights.append(mass[atom['elem'].iloc[0]])
                else:
                    print(f"# WARNING: Atom '{atom_name}' not found in residue {resi}. Using default values.")
                    B += 0.0
                    pos.append([0.0, 0.0, 0.0])
                    weights.append(1.0)

            bead_dict['B'] = B / len(CG[bead_name][types[i]])
            bead_pos = np.average(pos, weights=weights, axis=0)
            bead_dict['x'] = bead_pos[0]
            bead_dict['y'] = bead_pos[1]
            bead_dict['z'] = bead_pos[2]
            output_dict.append(bead_dict)

    return output_dict



# ======================== #
# ------------------------ #
#    MEMBRANE EMBEDDING    #
# ------------------------ #
# ======================== #


# Build CG membrane
# ------------------

def build_membrane_CG(ligands, CG_output, outputfile, command_line, mempro_additional, memprod_additional, membrane_comp, ion_conc, database=df, newlipidome = False, num_CPUs = 1):
    ''' Build membrane using MemPrO, optionally MemPrOD, and Insane4MemPrO. '''
    
    os.environ['NUM_CPU'] = str(num_CPUs) # This may be best set by the user
    
    CG_ligands = [database.at[name, 'CGName'] for name in ligands]

    MemPrO_dir = '.'.join(outputfile.split('.')[:-1])+'_MemPrO'
    
    MemPrO = ['MemPrO', '-f', CG_output, '-res', ','.join(set(CG_ligands)), '-o', MemPrO_dir]

    os.environ['PATH_TO_MARTINI'] = newmartini if newlipidome else oldmartini
    os.environ['PATH_TO_INSANE']  = oldinsane # Not used, but needs to be set.

    try:
        insane_params = get_command_line_parameters(command_line, 'mem')
    except IndexError:
        # at2mem doesn't require the same flag
        insane_params = np.array([])
    
    extra_bd_args = ''  # bd_args all need to be passed together

    dual_membrane = False

    if mempro_additional:
        build_sys = get_command_line_parameters(command_line, 'mp')
        SetFlags = {'-f' : [1, 'input file'], '--file_name' : [1, 'input file'],
                    '-res' : [1, 'additional CG ligands'], '--additional_residues' : [1, 'additional CG ligands'],
                    '-o' : [1, 'output directory'], '--output' : [1, 'output directory'],
                    '-bd_args' : [1], '--build_arguments': [1]}
        # Remove any pre-set flags - and separate out bd_args
        for flag in SetFlags.keys():
            if len(np.where(build_sys==flag)[0])!=0:
                if flag == '-bd_args' or '--build_arguments':
                    extra_bd_args += build_sys[np.where(build_sys == flag)[0][0]+1]
                else:
                    print('# WARNING: You have tried to overwrite the {} passed to MemPrO. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1])) 
                build_sys = np.hstack((build_sys[:np.where(build_sys == flag)[0][0]], build_sys[np.where(build_sys == flag)[0][0]+1+SetFlags[flag][0]:]))

        if len(np.where(build_sys=='-dm')[0])!=0 or len(np.where(build_sys=='--dual_membrane')[0])!=0:
            # Need to separately pass to Insane4MemPrO
            dual_membrane = True
                
        if len(np.where(build_sys=='-bd')[0])!=0 or len(np.where(build_sys=='--build_system')[0])!=0:
            flag = '-bd' if len(np.where(build_sys=='-bd')[0])!=0 else '--build_system'
            # Making the reasonable assumption of just one flag
            print('# WARNING: You have tried to generate a membrane embedded system directly via MemPrO. By default the top-ranked output will be used to create the output directly via Insane4MemPrO. This command has therefore been ignored.') 
            build_sys = np.hstack((build_sys[:np.where(build_sys == flag)[0][0]], build_sys[np.where(build_sys == flag)[0][0]+1+SetFlags[flag][0]:]))
                
        MemPrO.extend(list(build_sys))
                    
    extra_flags = {'-ni' : '15', '-ng' : '5', '-res_itp': os.environ['PATH_TO_MARTINI']}

    for flag in extra_flags.keys():
        if MemPrO.count(flag)==0:
            MemPrO.extend([flag, extra_flags[flag]])    
            
    result = subprocess.run(MemPrO)

    # Test output
    assert result.returncode==0, "ERROR: Failed to run MemPrO, please check for errors in your input"

    # Remove dummy membrane - for MemPrOD and Insane4MemPrO
    oriented = MemPrO_dir+'/Rank_1/'+'.'.join(outputfile.split('/')[-1].split('.')[:-1])+'_oriented.pdb'
    subprocess.run(['sed', '/DUM/d', MemPrO_dir+'/Rank_1/oriented_rank_1.pdb'], stdout=open(oriented, 'w'))
    
    all_insane = ['python', newinsane] if newlipidome else ['python', oldinsane]
    all_insane.extend(['-f', oriented, '-p', MemPrO_dir+'/Rank_1/CG_System_rank_1/topol.top', '-o', MemPrO_dir+'/Rank_1/CG_System_rank_1/CG-system.gro'])

    # Now optionally run MemPrOD    
    if memprod_additional:
        
        MemPrOD = ['MemPrOD', '-f', oriented, '-res', ','.join(set(CG_ligands)), '-o', MemPrO_dir+'/Rank_1/Deformations/']
        
        deform_sys = get_command_line_parameters(command_line, 'mdef')
        SetFlags = {'-f' : [1, 'input file'], '--file_name' : [1, 'input file'],
                    '-res' : [1, 'additional CG ligands'], '--additional_residues' : [1, 'additional CG ligands'],
                    '-o' : [1, 'output directory'], '--output' : [1, 'output directory']}
        # Remove any pre-set flags
        for flag in SetFlags.keys():
            if len(np.where(deform_sys==flag)[0])!=0:
                print('# WARNING: You have tried to overwrite the {} passed to MemPrOD. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1])) 
                build_sys = np.hstack((deform_sys[:np.where(deform_sys == flag)[0][0]], deform_sys[np.where(deform_sys == flag)[0][0]+1+SetFlags[flag][0]:]))                
                
                
        extra_flags = {'-res_itp': os.environ['PATH_TO_MARTINI']} # Setting grid size and iterations can lead to errors

        for flag in extra_flags.keys():
            if MemPrOD.count(flag)==0:
                MemPrOD.extend([flag, extra_flags[flag]])    
                
        result = subprocess.run(MemPrOD)
        
        # Test output
        assert result.returncode==0, "ERROR: Failed to run MemPrOD, please check for errors in your input"

        all_insane.extend(['-def', MemPrO_dir+'/Rank_1/Deformations/Membrane_Data/'])
        
    # Now run Insane4MemPrO
    # insane_params contains arguments passed via the -mem flag
    # Compare with bd_args, insane_params
    
    extra_bd_args = extra_bd_args.split()
    if len(extra_bd_args)!=0:
        extra_bd_args = np.array(extra_bd_args)
        flags = np.array([flag for flag in extra_bd_args if flag[0]=='-']) # Preserves order

        for f in range(len(flags)-1):
            # Insane4MemPrO has short flags only => don't need to consider duplication
            next_flag_loc = np.where(extra_bd_args==flags[f+1])[0][0]
            if len(np.where(insane_params == flags[f])[0]) != 0:
                print('# ERROR: The flag {} has been specified via -mem and -mp. Taking value from -mem.'.format(flag))
            else:
                # No duplicate specification
                insane_params = np.append(insane_params, extra_bd_args[:next_flag_loc])
            # Delete this flag from consideration

            extra_bd_args = extra_bd_args[next_flag_loc:]
            
        # Deal with last flag
        if len(np.where(insane_params == flags[-1])[0]) != 0:
            print('# ERROR: The flag {} has been specified via -mem and -mp. Taking value from -mem.'.format(flags[-1]))
        else:
            # No duplicate specification
            insane_params = np.append(insane_params, extra_bd_args)
        
            
    SetFlags = {'-f' : [1, 'input file'], '-p' : [1, 'topology file'], '-o' : [1, 'output file']}
    # Remove any pre-set flags
    for flag in SetFlags.keys():
        if len(np.where(insane_params==flag)[0])!=0:
            print('# WARNING: You have tried to overwrite the {} passed to Insane4MemPrO. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1])) 
            insane_params = np.hstack((insane_params[:np.where(insane_params == flag)[0][0]], insane_params[np.where(insane_params == flag)[0][0]+1+SetFlags[flag][0]:]))                

    all_insane = np.append(np.array(all_insane), insane_params)

    # Test for basic system set-up parameters
    extra_flags = {'-sol' : 'W', '-negi_c0': 'CL', '-posi_c0': 'NA', '-l': 'POPC', '-ion_conc': ','.join([str(ion_conc)]*3)}
            
    for flag in extra_flags.keys():
        if flag=='-ion_conc' and len(np.where(all_insane==flag)[0])==1:
                print('# WARNING: Overwriting any ion concentration specified via -C/--conc with alternative specification passed via -mem/--membrane.')
        if len(np.where(all_insane==flag)[0])==0:
            all_insane = np.append(all_insane, [flag, extra_flags[flag]])

    # Test for extent of system
    where = {}
    read_in = False

    for dim in ['x', 'y', 'z']:
        where[dim] = np.where(all_insane=='-'+dim)[0]
        if len(where[dim])==0:
            # Manually determine extent of system in unspecified dimensions
            # Based on MemPrO determination
            
            if not read_in:
                get_pdb, _, _ = read_PDB(oriented)
                oriented_data = pd.DataFrame.from_dict(get_pdb, orient='columns')
                read_in       = True

            max_dim = oriented_data[dim].max() ; min_dim = oriented_data[dim].min()
            dim_len = 8+0.1*(max_dim-min_dim)
            if dim == 'z':
                dim_len += 2
            
            all_insane = np.append(all_insane, ['-'+dim, str(dim_len)])

    if dual_membrane:
        # Manually pass the inter-membrane distance
        # Check for clashes

        rank1_info = open(MemPrO_dir+'/Rank_1/info_rank_1.txt').read().split('\n')
        rank1_info = [line for line in rank1_info if 'distance' in line][0]
        IM_dist    = float(rank1_info.split(':')[1].split()[0])

        IM_Insane  = 0.05 * IM_dist # Half distance, nm not A

        ps_included = np.where(insane_params=='-ps')[0]
        
        if len(ps_included)!=0:
            print('# WARNING: You have overwritten the inter-membrane distance passed to Insane4MemPrO. The value calculated by MemPrO is {} A, and the value which has been passed is {} A. Insane4MemPrO will proceed using your specified value.'.format(IM_Insane, 20*float(insane_params[ps_included+1])))
        else:
            all_insane = np.append(all_insane, ['-ps', IM_Insane])

            
    if not os.path.exists(MemPrO_dir+'/Rank_1/CG_System_rank_1/'):
        os.makedirs(MemPrO_dir+'/Rank_1/CG_System_rank_1/')
            
    result = subprocess.run(all_insane)
    
    # Test output
    assert result.returncode==0, "ERROR: Failed to run Insane4MemPrO, please check for errors in your input"
    
    return None



# Convert CG membrane to atomsitic
# ---------------------------------

def convert_membrane_at(system_data, basename, command_line, CG2AT_additional, CG_name, executable, newlipidome=False):
    ''' Convert membrane from CG to atomistic via CG2AT-lite. '''

    cg2at_args = ['cg2at_lite', '-a', basename+'_nomem.pdb', '-c', basename+'_MemPrO/Rank_1/CG_System_rank_1/{}'.format(CG_name), '-loc', basename+'_CG2AT']

    if os.path.isdir(basename+'_CG2AT/INPUT'):
        print('# WARNING: It appears that there is already a CG2AT directory of name {}. This will be deleted to allow for current conversion.'.format(basename+'_CG2AT/'))
        shutil.rmtree(basename+'_CG2AT/')
    
    if CG2AT_additional:
        cg2at_extra = get_command_line_parameters(command_line, 'at')
        SetFlags = {'-a'   : [1, 'atomistic input file'], '-c' : [1, 'CG input file'],
                    '-loc' : [1, 'output directory']}
        # Remove any pre-set flags
        for flag in SetFlags.keys():
            if len(np.where(cg2at_extra==flag)[0])!=0:
                cg2at_extra = np.hstack((cg2at_extra[:np.where(cg2at_extra == flag)[0][0]], cg2at_extra[np.where(cg2at_extra == flag)[0][0]+1+SetFlags[flag][0]:]))
                print('# WARNING: You have tried to overwrite the {} passed to CG2AT-lite. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1])) 

        if len(np.where(cg2at_extra == '-fg')[0]) !=0:
            if cg2at_extra[np.where(cg2at_extra == '-fg')[0][0]+1] != 'martini_3-0_charmm36':
                print('# WARNING: There may be issues when not using martini 3')
        else:
            print('# WARNING: There may be issues when not using martini 3')

    else:
        if newlipidome:
            cg2at_extra = ['-fg', 'martini_3-1_new_lipidome_charmm36']
        else:
            cg2at_extra = ['-fg', 'martini_3-0_charmm36']
        cg2at_extra.extend(['-gmx', executable])
            
    cg2at_args.extend(list(cg2at_extra))
        
    result = subprocess.run(cg2at_args)
    assert result.returncode==0, "ERROR: Failed to run CG2AT-lite, please check for errors in your input"

    return None
    


# ========================= #
# ------------------------- #
#    OUTPUT FILE WRITING    #
# ------------------------- #
# ========================= #


# Write line of PDB file
# -----------------------

def write_PDB_atom_line(f, counter, data):
    ''' Write a single line of a PDB file. '''
    
    f.write('{:<6}{:>5} '.format(data['entry'], counter))
    if len(data['name'])==4:
        f.write('{:<4} '.format(data['name']))
    else:
        f.write(' {:<3} '.format(data['name']))

    if len(data['resnm'])<=3:
        f.write('{:>3} '.format(data['resnm']))

    else:
        f.write('{:>4}'.format(data['resnm']))

    if len(data['chain'])==0:
        f.write(' {:>4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}'.format(data['resi'], data['x'],   data['y'],
                                                                            data['z'],    data['occ'], data['B']))
    else:
        f.write('{}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}'.format(data['chain'], data['resi'], data['x'],
                                                                             data['y'], data['z'], data['occ'], data['B']))
    if len(data['elem'])!=0:
        f.write('          {:>2}'.format(data['elem']))
    f.write('\n')
        
    return None



# Write full PDB file
# --------------------

def write_PDB(name, ordered_dict, title=[], cryst=[], ligand_chains=False):
    ''' Write full PDB file. '''

    f = open(name, 'w')

    ordered_dict = check_residue_number(ordered_dict)
    
    header = title
    header.extend(cryst)
    
    # Write title and UC information, if applicable
    if header != None:
        for line in header:
            f.write(line+'\n')

    j = 1
    for i, line in enumerate(ordered_dict):
        write_PDB_atom_line(f, j, line)
        j += 1
        if ligand_chains:
            if ((i+1) == len(ordered_dict) or (line['chain'] != ordered_dict[i+1]['chain'])):
                f.write('TER   {:>5}      '.format(j))
                if len(line['resnm'])<=3:
                    f.write('{:>3} '.format(line['resnm']))
                else:
                    f.write('{:>4}'.format(line['resnm']))

                f.write('{}{:>4}\n'.format(line['chain'], line['resi']))
                j += 1 
    f.write('END')
    f.close()

    return None



# Ensure residue numbers PDB-compatible
# --------------------------------------

def check_residue_number(ordered_dict):
    ''' Renumber residues where this exceeds max PDB can handle. '''
    
    dataframe = pd.DataFrame.from_dict(ordered_dict, orient='columns')
    
    if max(dataframe['resi']) > 9999:
        print('# WARNING: Number of residues exceeds 9999, which is a limitation of the PDB format. Renumbering all residues with high residue IDs from 0 (will repeat as often as needed).')
        
        while max(dataframe['resi']) > 9999:
            dataframe['resi'] = dataframe['resi'].apply(lambda x: (x if x < 10000 else x - 10000))
                
    return dataframe.to_dict('records')



# Get PDB box from GRO file
# --------------------------

def convert_vectors(box_vecs):
    ''' Convert from GRO box format to PDB box format. '''

    print('# INFO: Converting unit cell vectors into a, b, c, alpha, beta, gamma format.')
    
    box_vecs = box_vecs.split()
    box_vecs = np.array([float(vec) if abs(float(vec)) > 1e-8 else 0 for vec in box_vecs])

    # Test for orthogonal axes
    if len(box_vecs) == 3:
        # Only first 3 values given
        a = box_vecs[0] ; b = box_vecs[1] ; c = box_vecs[2]
        alpha = 90.0    ; beta = 90.0     ; gamma = 90.0

    elif sum(abs(box_vecs[3:])) == 0:
        # Last 6 values all 0
        a = box_vecs[0] ; b = box_vecs[1] ; c = box_vecs[2]
        alpha = 90.0    ; beta = 90.0     ; gamma = 90.0

    else:
        # Not orthogonal unit cell

        # Determine length of box vectors
        a = np.sqrt(box_vecs[0]*box_vecs[0] + box_vecs[3]*box_vecs[3] + box_vecs[4]*box_vecs[4])
        b = np.sqrt(box_vecs[1]*box_vecs[1] + box_vecs[5]*box_vecs[5] + box_vecs[6]*box_vecs[6])
        c = np.sqrt(box_vecs[2]*box_vecs[2] + box_vecs[7]*box_vecs[7] + box_vecs[8]*box_vecs[8])

        # Create unit vectors
        A = np.array([box_vecs[0], box_vecs[3], box_vecs[4]])/a
        B = np.array([box_vecs[5], box_vecs[1], box_vecs[6]])/b
        C = np.array([box_vecs[7], box_vecs[8], box_vecs[2]])/c
        
        # Determine angles
        convert = 180/np.pi
        
        alpha = np.arccos(np.dot(B, C)) * convert
        beta  = np.arccos(np.dot(A, C)) * convert
        gamma = np.arccos(np.dot(A, B)) * convert


    # nm -> A 
    cryst = 'CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P 1           1'.format(a*10, b*10, c*10, alpha, beta, gamma)
    
    return [cryst]


# ================================== #
# ---------------------------------- #
#    GLOBULAR PROTEIN PREPARATION    #
# ---------------------------------- #
# ================================== #


# Solvate globular protein environment
# -------------------------------------

def make_glob(system, command_line, executable, simtype, conc, topol='topol.top'):
    ''' Place a globular protein into a box, and solvate and add ions to allow for subsequent simulation. '''

    # CG following https://github.com/marrink-lab/martini-workshop/blob/main/02_protein_basics/tutorial.md
    
    outputdir  = '/'.join(system.split('/')[:-1]) 
    sysname    = system.split('/')[-1]
    workingdir = os.getcwd()

    # Place protein into a box - for AA and CG
    # -----------------------------------------
    
    SetFlags    = {'-f' : [1, 'input file']} # Don't allow overwriting file name
    
    gro_editconf = [executable, 'editconf', '-f', sysname]
    
    # Generate defaults
    # -----------------
    box = '.'.join(sysname.split('.')[:-1])+'_box.pdb'

    # Generate editconf parameters from command line
    editconf_extra = get_command_line_parameters(command_line, 'B')
    
    if len(editconf_extra)==0:
        gro_editconf.extend(['-o', box, '-c', '-d', '2'])
    else:
        # Remove any pre-set flags
        for flag in SetFlags.keys():
            if len(np.where(editconf_extra==flag)[0])!=0:
                editconf_extra = np.hstack((editconf_extra[:np.where(editconf_extra == flag)[0][0]], editconf_extra[np.where(editconf_extra == flag)[0][0]+1+SetFlags[flag][0]:]))
                print('# WARNING: You have tried to overwrite the {} passed to editconf. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))

        # Ensure name is set
        if len(np.where(editconf_extra=='-o')[0])==0:
            # This flag has not been set, add to editconf_extra
            editconf_extra = np.hstack((editconf_extra, np.array(['-o', box])))
        else:
            box = editconf_extra[np.where(editconf_extra=='-o')[0][0]+1]
            
        gro_editconf.extend(list(editconf_extra))
            
    if outputdir != '':
        os.chdir(outputdir)
        
    gmx = subprocess.run(gro_editconf)
    assert gmx.returncode ==0, ' ERROR: GROMACS could not create a box.'
            

    # Solvate box - for AA and CG
    # ----------------------------
    
    SetFlags = {'-cp' : [1, 'input file'], '-p' : [1, 'topology file']} # Don't allow overwriting file name or topology

    # Generate defaults
    # -----------------
    
    gro_solvate = [executable, 'solvate', '-cp', box, '-p', topol]
    solv        = '.'.join(sysname.split('.')[:-1])+'_solv.pdb'
    
    if simtype == 'AA':
        solvate_extra = ['-cs', '-o', solv]
        defaults      = {'-o': solv, '-cs': ''}    
    elif simtype == 'CG':
        solvate_extra = ['-cs', 'martini_water.gro', '-radius', '0.21', '-o', solv]
        defaults      = {'-o': solv, '-cs': 'martini_water.gro'}    

    solv_extra = get_command_line_parameters(command_line, 'SV')
        
    if len(solv_extra)==0:
        gro_solvate.extend(solvate_extra)
        if simtype == 'CG':
            subprocess.run(['scp', '-r', mdpPath+'martini_water.gro', '.'])
        
    else:
        # Remove any pre-set flags
        for flag in SetFlags.keys():
            if len(np.where(solv_extra==flag)[0])!=0:
                solv_extra = np.hstack((solv_extra[:np.where(solv_extra == flag)[0][0]], solv_extra[np.where(solv_extra == flag)[0][0]+1+SetFlags[flag][0]:]))
                print('# WARNING: You have tried to overwrite the {} passed to solvate. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))
                
        # Ensure name is set
        if len(np.where(solv_extra=='-o')[0])==0:
            # This flag has not been set, add to solv_extra
            solv_extra = np.hstack((solv_extra, [flag, defaults[flag]]))
        else:
            solv = solv_extra[np.where(solv_extra=='-o')[0][0]+1]

        # Check for '-cs' flag
        if len(np.where(solv_extra=='-cs')[0])==0:
            # This flag has not been set, add to solv_extra
            solv_extra = np.hstack((solv_extra, [flag, defaults[flag]]))
            if simtype == 'CG':
                subprocess.run(['scp', '-r', mdpPath+'martini_water.gro', '.'])
        gro_solvate.extend(list(solv_extra))
        
    gmx = subprocess.run(gro_solvate)
    assert gmx.returncode ==0, ' ERROR: GROMACS could not solvate the system. '
            
            
    # Add ions
    # ---------

    subprocess.run(['scp', '-r', mdpPath+'genion.mdp', '.'])
    print('# WARNING: Ignoring warnings for creating a tpr to add ions')
    ion_tpr = [executable, 'grompp',  '-f', 'genion.mdp', '-c', solv, '-p', topol, '-o', 'ion.tpr', '-maxwarn', '50']
    gmx = subprocess.run(ion_tpr)
    assert gmx.returncode ==0, ' ERROR: GROMACS could not create a tpr to add ions'
        
    SetFlags = {'-s' : [1, 'tpr file'], '-p' : [1, 'topology file']} # Don't allow overwriting tpr or topology

    # Generate defaults
    # -----------------

    outname  = '.'.join(sysname.split('.')[:-1])+'_ions_{}.pdb'
    gmx_ions = [executable, 'genion', '-s', 'ion.tpr', '-p', topol]

    genion_extra = get_command_line_parameters(command_line, 'I')
    
    if len(genion_extra)==0:
        ions_extra = ['-neutral', '-conc', str(conc), '-o', outname.format(conc)]
        gmx_ions.extend(ions_extra)
        output = outname.format(conc)
    else:
        extra_flags = {'-conc': conc, '-o': outname}
        
        if len(np.where(genion_extra=='-conc')[0])==1:
            print('# WARNING: Overwriting any ion concentration specified via -C/--conc with alternative specification passed via -I/--ions.')
            conc = genion_extra[np.where(genion_extra=='-conc')[0][0]+1]
        else:
            genion_extra = np.hstack((genion_extra, np.array(['-conc', str(conc)])))

                
        if len(np.where(genion_extra=='-o')[0])==1:
            output = genion_extra[np.where(genion_extra=='-o')[0][0]+1]
        else:
            output       = outname.format(conc)
            genion_extra = np.hstack((genion_extra, np.array(['-o', output])))
            
        gmx_ions.extend(list(genion_extra))


    gmx = subprocess.Popen(gmx_ions, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out = gmx.communicate(input='SOL')
    print(out[1])
    assert out[0] ==None, ' ERROR: grompp could not add ions.'
    
    print('# INFO: Final output is {}.'.format(output))
    
    return output
                  

    
# =================== #
# ------------------- #
#    MD GENERATION    #
# ------------------- #
# =================== #


# Make and optionally run tpr files
# ----------------------------------

def make_gmx(system, command_line, tpr_type, executable, ndx=None, run=False, topol='topol.top'):
    ''' Generate (and optionally run) one of a variety of gromacs tpr files. '''

    outputdir  = '/'.join(system.split('/')[:-1]) 
    sysname    = system.split('/')[-1]
    workingdir = os.getcwd()

    gro_grompp = [executable, 'grompp', '-c', sysname, '-p', topol]
        
    SetFlags = {'-c' : [1, 'input file'], '-p' : [1, 'topology file']} # Don't allow overwriting file name
    
    path_flags = ['-c', '-r', '-rb', '-t', '-e', '-qmi']
    
    
    if ndx==None:
        path_flags.append('-n')
    else:
        # Utilise created index file
        gro_grompp.extend(['-n', ndx])
        SetFlags['-n'] = [1, 'index file']
        
    # Generate defaults
    # -----------------
    outputname = '.'.join(system.split('/')[-1].split('.')[:-1])
    defaults = {'-f' : mdp[tpr_type], '-o': outputname+'_'+tpr_type+'.tpr'}

    # Add default restrain file for AAeq
    if tpr_type.split('_')[0] == 'AAeq':
        defaults['-r'] = sysname

    # Generate grompp parameters from command line
    grompp_extra = get_command_line_parameters(command_line, tpr_type.split('_')[0])
    
    # Remove any pre-set flags
    for flag in SetFlags.keys():
        if len(np.where(grompp_extra==flag)[0])!=0:
            grompp_extra = np.hstack((grompp_extra[:np.where(grompp_extra == flag)[0][0]], grompp_extra[np.where(grompp_extra == flag)[0][0]+1+SetFlags[flag][0]:]))
            print('# WARNING: You have tried to overwrite the {} passed to grompp. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))

    # Check flags which are path-dependent
    for flag in path_flags:
        flag_loc = np.where(grompp_extra==flag)[0]
        if len(flag_loc)!=0:
            grompp_extra[flag_loc[0]+1] = workingdir+'/'+grompp_extra[flag_loc[0]+1]
            
    actual     = {'-f': False, '-o': False}
    # Ensure name and mdp file are set
    for flag in defaults.keys():
        if len(np.where(grompp_extra==flag)[0])==0:
            # This flag has not been set, add to grompp_extra
            grompp_extra = np.hstack((grompp_extra, [flag, defaults[flag]]))
            actual[flag] = True
            
    tpr_name = grompp_extra[np.where(grompp_extra=='-o')[0][0]+1]
    gro_name = '.'.join(tpr_name.split('.')[:-1])+'.gro'

    if outputdir != '':
        os.chdir(outputdir)
    
    if actual['-f']:
        # The default tpr has been used and should be copied, need flag to prevent
        # overwriting other file of same name
        subprocess.run(['scp', '-r', mdpPath+mdp[tpr_type], '.']) # All other files should exist
        
    gro_grompp.extend(grompp_extra)

    
    gmx = subprocess.run(gro_grompp)
    assert gmx.returncode ==0, ' ERROR: grompp could not generate a tpr file.'
    
    # Run tpr if applicable

    if run:
        gro_mdrun = [executable, 'mdrun', '-deffnm', '.'.join(tpr_name.split('.')[:-1])]

        # Generate command line parameters
        mdrun_extra = get_command_line_parameters(command_line, 'r'+tpr_type.split('_')[0])
        
        SetFlags = {'-deffnm' : [1, 'input file']} # Don't allow overwriting file name
        # Remove any pre-set flags
        for flag in SetFlags.keys():
            if len(np.where(mdrun_extra==flag)[0])!=0:
                mdrun_extra = np.hstack((mdrun_extra[:np.where(mdrun_extra == flag)[0][0]], mdrun_extra[np.where(mdrun_extra == flag)[0][0]+1+SetFlags[flag][0]:]))
            print('# WARNING: You have tried to overwrite the {} passed to mdrun. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))

        if len(np.where(mdrun_extra == '-c')[0])!=0:
            # Changing the name of the output file
            gro_name = mdrun_extra[np.where(mdrun_extra == '-c')[0][0]+1]
            
        gro_mdrun.extend(mdrun_extra)
            
        gmx = subprocess.run(gro_mdrun)
        assert gmx.returncode ==0, ' ERROR: tpr file would not run.'
    
    os.chdir(workingdir)
                    
    return gro_name



# ============== #
# -------------- #
#    STEERING    #
# -------------- #
# ============== #


# Extract reference positions for ligands to steer
# --------------------------------------------------

def steer_reference(reffiles, st_ligs, steerpath, executable):
    ''' Extract the reference ligand positions. '''

    reflig_locs = dict.fromkeys(st_ligs)
    findref     = st_ligs
    
    # ==================== #
    # Read reference files # 
    # ==================== #
    
    if reffiles != None:
        for reffile in reffiles:
            # Read files
            # ----------
            
            if reffile[-3:] == 'gro':
                tmp, _, _ = read_GRO(args.inputfile)
            elif reffile[-3:] == 'pdb':
                tmp, _, _ = read_PDB(args.inputfile)
            elif reffile[-3:] == 'cif':
                tmp, _ = read_CIF(args.inputfile)
                tmp = [atom for model in tmp for atom in model]
            else:
                print('# WARNING: Ligand reference file must be PDB/GRO/CIF but file {} is neither. This will be ignored.'.format(reffile))
    
            ref_data = pd.DataFrame.from_dict(tmp, orient='columns')            
            refligs  = ref_data['resnm'].unique()

            # Consider each ligand in file
            # ----------------------------
            
            for reflig in refligs:
                reflig = reflig.upper()
                if len(np.where(st_ligs == reflig)[0])==0:
                    # Not a ligand to steer
                    continue

                # Test for multiply defined ligands
                # ----------------------------------
                
                if reflig_locs[reflig] != None:
                    print('ERROR: Reference structure for ligand {} is mutliply defined - in {} and {}.'.format(reflig, reflig_locs[reflig], reffile))
                    sys.exit()
                    
                # Check that there is only one instance of the ligand - assume on chain/resid
                liglocs = ref_data.loc[data['resnm']==reflig]
                if len(liglocs.groupby(['resi', 'chain']).size()) > 1:
                    print('ERROR: Reference structure for ligand {} is mutliply defined in {}.'.format(reflig, reffile))
                    sys.exit()
                
                # Create appropriate referernce files
                # -----------------------------------
                    
                reflig_locs[reflig] = reffile
                findref.remove(reflig)
                print('# INFO: Creating steering reference file for {} from {}.'.format(reflig, reffile))
                create_ref_files(liglocs, reflig, steerpath, executable)
                
    # =========== #
    # Read CCD2MD # 
    # =========== #

    # First want to check if there are any unassigned ligands
    if len(findref)==0:
        return None

    ccdcif, _ = read_CIF(CCD2MD_dir+'CCD2MD.cif')
    ccdcif = [atom for ligand in ccdcif for atom in ligand]
    ccdcif = pd.DataFrame.from_dict(ccdcif, orient='columns')

    ccdligs  = ccdcif['resnm'].unique() # Guaranteed to be only one copy
        
    # Test each remaining ligand
    for lig in findref:
        if len(np.where(ccdligs == lig)[0])==0:
            print('ERROR: Attempting to steer {}, but no reference structure was found in input files or CCD2MD.cif'.format(lig))
            exit()

        ligloc = ccdcif.loc[ccdcif['resnm']==lig]
        print('# INFO: Creating steering reference file for {} from CCD2MD.cif.'.format(lig))
        create_ref_files(ligloc, lig, steerpath, executable)

    return None
    
    

# Create reference files 
# ------------------------

def create_ref_files(ref_pos, lig, steerpath, executable):
    ''' Create all necessary reference files for ligand.'''

    # Create directory and copy reference PDB
    # ----------------------------------------

    refpath = steerpath+'/steer_{}/'.format(lig)
    if not os.path.exists(refpath):
        os.makedirs(refpath)
    else:
        print('# WARNING: Steering output directory {} exists. Any relevant files will be overwritten.'.format(refpath))

    subprocess.run(['scp', '-r', 'charmm36-ccd2md.ff', refpath])
        
    write_PDB(refpath+'{}_ref.pdb'.format(lig), ref_pos)

    # Generate topology
    # -----------------
    
    subprocess.run([executable, 'pdb2gmx', '-f', refpath+'{}_ref.pdb'.format(lig), '-o',
                    refpath+'{}_ref.pdb'.format(lig), '-p', refpath+'{}_ref.top'.format(lig), '-ff',
                    'charmm36-ccd2md', '-water', 'tip3p'])

    subprocess.run([executable, 'editconf', '-f', refpath+'{}_ref.pdb'.format(lig), '-box', '10', '10', '10',
                    '-c', '-o', refpath+'{}_ref.pdb'.format(lig)])


    # Write index file
    # -----------------
    
    # Note, need to read from the pdb2gmx output in case Hs missing
    
    refH, _, _ = read_PDB(refpath+'{}_ref.pdb'.format(lig))
    
    index = open(refpath+'{}.ndx'.format(lig), 'w')
    index.write('[ Heavy ]')
    index.write('\n')
    
    for i, atom in enumerate(refH):
        try:
            if atom['elem'] != 'H' and atom['elem'][0] != 'H':
                index.write(str(i+1)+'\t')
        except (KeyError, IndexError):
            # Either some or all element names missing
            global element_name
            if element_name:
                print('# WARNING: Element names are missing - attempting to infer from atom names. Note that this may cause issues.')
                element_name = False
            elem = FuncConv.determine_element(atom['name'])
            if elem != 'H':
                index.write(str(i+1)+'\t')
                
    index.close()
    
    return None
    


# Extract predicted positions for ligands
# ----------------------------------------

def get_ligand(output_data, ligand, loc):
    ''' Extract the ligands to steer from the simulation output. '''
    
    # Get current positions of ligands to steer                                                                   
    all_lig = [atom for atom in output_data if atom['resnm'] == ligand]
    
    if len(all_lig) == 0:
        return None, None
    
    test_lig  = all_lig
    lig_IDs   = {}
    to_steer  = []
    
    while len(test_lig)!=0:
        curr_chain = test_lig[0]['chain']
        curr_resi  = test_lig[0]['resi']
        
        # Determine if multiple different instances or not                                                    
        curr_lig = [line for line in test_lig if (line['chain'] == curr_chain) and (line['resi'] == curr_resi)]
        test_lig = [line for line in test_lig if (line['chain'] != curr_chain) or (line['resi'] != curr_resi)]

        
        curr_chain = curr_chain if len(curr_chain) == 0 else curr_chain+'_'

        lig_IDs['steer_{}_{}{}.pdb'.format(ligand, curr_chain, curr_resi)] = [curr_lig[0]['atomi'], curr_lig[-1]['atomi']]
                
        # Need to renumber atoms to correspond to 1-N and keep a record of this
                                                      
        for i, atom in enumerate(curr_lig):
            atom['atomi'] = i
        
        write_PDB(loc+'steer_{}_{}{}.pdb'.format(ligand, curr_chain, curr_resi), curr_lig)
        to_steer.append('steer_{}_{}{}.pdb'.format(ligand, curr_chain, curr_resi))
        
    return to_steer, lig_IDs



# Perform steering on single file
# --------------------------------

def steer_atomistic(predicted_file, ligand, executable):
    ''' Steer reference file towards prediction. '''
        
    # Add prediction to box
    # ----------------------

    gmx = subprocess.run([executable, 'editconf', '-f', predicted_file, '-box', '10', '10', '10', '-o', predicted_file])
    assert gmx.returncode ==0, ' ERROR: GROMACS could not create a box for {}.'.format(predicted_file)

    # Get restraints
    # --------------
    
    gmx = subprocess.run([executable, 'genrestr', '-f', predicted_file, '-o', 'posre.itp', '-n', '{}.ndx'.format(ligand)])
    # Only one entry in index file => no interactivity needed
    assert gmx.returncode ==0, ' ERROR: GROMACS could not create restraints for {}.'.format(predicted_file)

    # Perform steering
    # -----------------

    global steer_warn
    if steer_warn:
        print('# WARNING: Ignoring GROMACS warnings for steering. This may obscure errors in input files.')
        steer_warn = False
        
    steer_name = predicted_file.split('.')[0]
    
    gmx = subprocess.run([executable, 'grompp', '-f', mdpPath+'steer.mdp', '-c', predicted_file, '-r',
                          ligand+'_ref.pdb', '-p', ligand+'_ref.top', '-o', steer_name, '-maxwarn', '10'])
    assert gmx.returncode ==0, ' ERROR: GROMACS could not create tpr file for {}.'.format(predicted_file)    
    gmx = subprocess.run([executable, 'mdrun', '-deffnm', steer_name])
    assert gmx.returncode ==0, ' ERROR: GROMACS could not steer {}.'.format(predicted_file)
                        
    return steer_name+'.gro'




# Perform steering on all files for a ligand
# -------------------------------------------

def steer_ligand(ligand, steer_loc, steer_files, steerIDs, output, executable):
    ''' Perform steering for all files for a ligand. '''
    
    workingdir = os.getcwd()
    # Change to the output directory
    os.chdir(steer_loc)
    
    # Add prediction to box
    # ----------------------

    for fle in steer_files:
        steernm = steer_atomistic(fle, ligand, executable)
        steered, _, _ = read_PDB(steernm)

        for atomid in range(len(steered)):
            output[atomid+ligand_IDs[fle][0]-1]['x'] = steered['x']
            output[atomid+ligand_IDs[fle][0]-1]['y'] = steered['y']
            output[atomid+ligand_IDs[fle][0]-1]['z'] = steered['z']

    os.chdir(workingdir)
                    
    return output



# ========================= #
# ------------------------- #
#    TOPOLOGY GENERATION    #
# ------------------------- #
# ========================= #


# Generate topology for CG systems
# ---------------------------------

def get_topology_CG(outputfile, membrane, ligands, prot, inputfile, newlipidome=False, topoldir=None, database=df):
    ''' Write topology files for CG systems. '''

    if prot:
        martinize       = open('.'.join(outputfile.split('.')[:-1])+'_proteinCG.top', 'r').read()
        martinize_lines = martinize.split('\n')

        includes = [line for line in martinize_lines if len(line)!=0 and line[0] == '#']
        includes = [line for line in includes if line[:17] != '#include "martini']        
        files    = [line.split('"')[1] for line in includes]
        
    outputdir = '/'.join(outputfile.split('/')[:-1]) if topoldir is None else topoldir
    if len(outputdir)!=0:
        # Copy molecule files into correct directory
        outputdir += '/'
        for f in files:
            subprocess.run(['scp', f, outputdir])
    else:
        outputdir = '.'

    os.environ['PATH_TO_MARTINI'] = newmartini if newlipidome else oldmartini
    
    copy_martini = ['scp', os.environ['PATH_TO_MARTINI'], outputdir]
    subprocess.run(copy_martini)

    loc = outputdir+'topol.top' if outputdir != '.' else 'topol.top'
    
    if membrane:
        # Generate topology from MemPrO output
        # ------------------------------------

        # Should have a protein for membrane association
        
        # Replace MemPrO protein-cg.itp with martinize inputs
        copy     = ['scp', '.'.join(outputfile.split('.')[:-1])+'_MemPrO/Rank_1/CG_System_rank_1/topol.top', outputdir]
    
        subprocess.run(copy)

        includes = '\\n'.join(includes)
        # Note: assuming that MemPrO will NOT have the same martini file(s) as here

        # Get correct number of molcules
        mols = martinize.split('[ molecules ]')[-1].split('\n')
        mols = [mol for mol in mols if len(mol)!=0]

        # MemPrO as default does not include additional ligands - add here
        for lig in set([l for l in ligands if l not in PTMs]):
            num = len([l for l in ligands if l==lig])
            CG_name = database.at[lig.strip(), 'CGName']
            mols.extend(['{} \t{}'.format(CG_name, num)])
        
        # Note: checking for Mac/Linux as the sed behaviour is different

        curr_sys = platform.system()
        if curr_sys == 'Darwin' or curr_sys == 'darwin':
            sedflag     = ['sed', '-i', '']
            martiniflag = [r'/"protein-cg.itp"/i\ \n#include "{}"\n'.format(copy_martini[1].split('/')[-1]), loc]
            proteinflag = ['s/Protein *1/'+'\\n'.join(mols)+'/g', loc]
        else:
            sedflag     = ['sed', '-i']
            martiniflag = ['/"protein-cg.itp"/i #include "{}"'.format(copy_martini[1].split('/')[-1]), loc]
            proteinflag = [r's/Protein\s*1/'+'\\n'.join(mols)+r'/g', loc]
            
        subprocess.run(sedflag + ['/"martini.*itp"/d', loc])
        subprocess.run(sedflag + martiniflag)
        subprocess.run(sedflag + ['s/#include "protein-cg.itp"/'+includes+'/g', loc])
        subprocess.run(sedflag + proteinflag) # Changes protein name and adds ligands
        
        print('# INFO: Topology file has been created based on martinize2 and Insane4MemPrO outputs.')

    else:
        # Create topology from martinize2 output, if exists
        # -------------------------------------------------
        
        topol = open(loc, 'w')
        
        if prot:
            include_martini = False
            topol_gen       = ' based on martinize2 output'
            
            for i, line in enumerate(martinize_lines):
                if line[:17] != '#include "martini':
                    # Insert correct version of martini itp file(s)
                    if not (i==len(martinize_lines)-1 and len(line.strip())==0):
                        topol.write(line+'\n')
                elif not include_martini:                
                    topol.write('#include "{}"\n'.format(copy_martini[1].split('/')[-1]))
                    include_martini = True

            
        else:
            # No martinize information => write file
            topol.write('\n#include "{}"\n'.format(copy_martini[1].split('/')[-1]))
            topol.write('\n[ system ]\n')
            topol.write('CG system created by CCD2MD from {}\n'.format(inputfile))
            topol.write('\n[ molecules ]\n')
            topol_gen = ''
            
        # Add topology for additional ligands
        # -----------------------------------
        
        for lig in set([l for l in ligands if l not in PTMs]):
            num = len([l for l in ligands if l==lig])
            CG_name = database.at[lig.strip(), 'CGName']
            topol.write('{} \t{}\n'.format(CG_name, num))

        topol.close()
        print('# INFO: Topology file has been created{}.'.format(topol_gen))
        
    return None



# Generate topology for atomistic systems
# ----------------------------------------

def get_topology_atomistic(outputfile, membrane, executable, at_command=None, output_data=None):
    ''' Write topology for atomistic systems. '''
    
    outputdir = '/'.join(outputfile.split('/')[:-1])+'/' if len(outputfile.split('/')[:-1])!=0 else './'
    
    if membrane:
        # Copy topology and itp files from output of CG2AT
        # -------------------------------------------------
        
        itpfiles = glob.glob('.'.join(outputfile.split('.')[:-1])+'_CG2AT/FINAL/*itp')
        subprocess.run(['scp', '.'.join(outputfile.split('.')[:-1])+'_CG2AT/FINAL/topol_final.top', outputdir+'topol.top'])
        
        for itp in itpfiles:
            subprocess.run(['scp', itp, outputdir])
        print('# INFO: Topology file generated by CG2AT-lite')

        curr_sys = platform.system()
        if curr_sys == 'Darwin' or curr_sys == 'darwin':
            subprocess.run(['sed', '-i', '', r's|../FINAL/charmm36-cg2at.ff|charmm36-cg2at.ff|g', 'topol.top'])
        else:
            subprocess.run(['sed', '-i', r's|../FINAL/charmm36-cg2at.ff|charmm36-cg2at.ff|g', 'topol.top'])
        
        final = outputfile
        topol = 'topol.top'
        
    else:
        # Generate topology and add hydrogens using pdb2gmx
        # -------------------------------------------------
        
        subprocess.run(['scp', '-r', CHARMMPath, outputdir])        
        os.chdir(outputdir)
        outputname = outputfile.split('/')[-1]
        
        if at_command is None:
            # No additional parameters for pdb2gmx passed - use defaults
            
            gromacs = [executable,  'pdb2gmx',  '-f', outputname, '-o',
                       '.'.join(outputname.split('.')[:-1])+'_H.pdb', '-water', 'tip3p', '-ff', 'charmm36-ccd2md']

            # Check for PTMs
            # --------------
            
            if len(output_data[output_data['resnm'].isin(terminal_PTMs)]) != 0:
                
                subprocess.run(['scp', 'charmm36-ccd2md.ff/residuetypes.dat', '.']) # Needed for PTMs                

                # Need to check for terminal PTMs as this changes the termini
                # Loop through chains and first and last residue in these

                termini = []
                
                chains = output_data.chain.unique()
                for chain in chains:
                    chain_data    = output_data.loc[output_data['chain']==chain]

                    # Assuming no C terminal modifications and 6 possible starting termini for terminal_PTMs
                    resnm = list(chain_data.loc[chain_data['resi']==min(chain_data['resi']), 'resnm'])[0]
                    # Assuming that there are 6 starting termini and 5 ending termini
                    if resnm in terminal_PTMs:
                        termini.extend(['6', '0'])
                    else:
                        termini.extend(['0', '0'])
                        
                    # Note gromacs doesn't seem to use stdout and some usage of stderr (but not complete)

                ter = '\n'.join(termini)+'\n'

                gromacs = [executable,  'pdb2gmx',  '-f', outputname, '-o',
                           '.'.join(outputname.split('.')[:-1])+'_H.pdb', '-water', 'tip3p',
                           '-ff', 'charmm36-ccd2md', '-ter']
            
            
                gmx = subprocess.Popen(gromacs, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                out = gmx.communicate(input=ter)
                print(out[1])
                assert out[0] == None, 'GROMACS could not run pdb2gmx'
                
            else:
                gmx = subprocess.run(gromacs)

            final = '.'.join(outputname.split('.')[:-1])+'_H.pdb'
            topol = 'topol.top'                
        else:
            # Generate command line input and pass to pdb2gmx
            # Command line options added after -gh or --pdb2gmx
            pdb2gmx_args = get_command_line_parameters(at_command, 'gh')

            SetFlags = {'-f' : [1, 'input file']} # Don't allow overwriting file name
            # Remove any pre-set flags
            for flag in SetFlags.keys():
                if len(np.where(pdb2gmx_args==flag)[0])!=0:
                    pdb2gmx_args = np.hstack((pdb2gmx_args[:np.where(pdb2gmx_args == flag)[0][0]], pdb2gmx_args[np.where(pdb2gmx_args == flag)[0][0]+1+SetFlags[flag][0]:]))
                    print('# WARNING: You have tried to overwrite the {} passed to pdb2gmx. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))
                        
            gromacs = [executable, 'pdb2gmx', '-f', outputname]
            gromacs.extend(pdb2gmx_args)
            
            gmx = subprocess.run(gromacs)

            final = 'conf.gro'
            topol = 'topol.top'
            
            if len(np.where(pdb2gmx_args=='-o')[0])!=0:
                final = pdb2gmx_args[np.where(pdb2gmx_args == '-o')[0][0]+1]
            if len(np.where(pdb2gmx_args=='-p')[0])!=0:
                topol = pdb2gmx_args[np.where(pdb2gmx_args == '-p')[0][0]+1]

            
        if gmx.returncode==0:
            print('# INFO: Topology file generated by pdb2gmx')
        else:
            print('# ERROR: pdb2gmx could not generate a topology file.')
    
    return final, topol

