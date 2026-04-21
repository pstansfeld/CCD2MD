#######################################
#                                     #
#   FUNCTIONS FOR CCD2MD - POS2CIF    #
#                                     #
#######################################

'''
General file for functions to write userCCD fields and AF3 inputs

Last Update: K Blow 21/04/26

Contains:

Input file parsers
-------------------

read_configuration_file(conffile, args, flag_dict=flag_dict) # Read configuration file and cmpare to command line inputs
read_fasta(fastafile)                                        # Determine protein information from FASTA file

get_elements(NAME)        # Determine element type from atom name
get_components(COMPONENT) # Split list of proteins/ligands to determine consituent parts
 

User CCD information generation
-------------------------------

get_charge(charge, charge_cutoff)       # Determine integer charge of atoms
get_bonds_dist(atom_data, bond_cutoff)  # Determine if atoms are bonded based on distance cutoff 

get_files(all_files)  # Split files into position and bonding files


Position/bonding file parsers
-----------------------------

read_mol2(mol2_data, name,  keywords=mol2_keywords)     # Get position and bonding information from mol2 file

read_position(data, name)            # Get position from non-mol2 file
read_rtp_bonding(data, name, charge) # Get bonding information from rtp file 
read_rtf_bonding(data, name, charge) # Get bonding information from rtf file
read_itp_bonding(data, name, charge) # Get bonding information from itp file


General conversion information
-------------------------------

CCD_from_name(name, position_files, bonding_files, bond_cutoff, charge, include_H=False) # Generate userCCD information for NAME

get_position(name, position_files)     # Scan all files for position information
get_bonds(name, bonding_files, charge) # Scan all files for bonding information

read_covalent_mol2(molfile, keywords=mol2_keywords) # Get covalently bonded molecule position and bonds from mol2 file
get_covalent_position(molfile, filetype='pdb')      # Get covalently bonded molecule position from pdb file
def get_covalent_bonding_itp(bondfl, atomic)        # Get covalently bonded molecule bonds from itp file


Output file writing
--------------------

cif_information(nname, descript, posdata)                                                # Write userCCD cif data
covalent_cif(ligname, atomic, bonds, components)                                         # Write userCCD cif data for covalently bonded ligands
extract_covalent(currlig, lig, json, inCIF, CIF, title, lig_chains)                      # Extract and write covalent ligand information
extract_ligand(currlig, lig_check, allchains, inCIF, CIF, json, title, lignums, ligfile) # Extract and write non-covalent ligand information

'''

__version__ = "1.1.1"


import sys
import numpy as np
import pandas as pd

# General parameters
# -------------------


flag_dict = [{'e'   : 'charge', 
              'b'   : 'bond', 
              'H'   : 'Hydrogen',
              'j'   : 'json',
              't'   : 'title',
              'A'   : 'afvers',
              'd'   : 'dialect',
              'CF'  : 'configuration',
              'nC'  : 'no_CCD2MD',
              'n'   : 'names', 
              'r'   : 'rename', 
              'f'   : 'files',
              'c'   : 'covalent',
              'l'   : 'ligand',
              'u'   : 'userCCDPath',
              's'   : 'seeds',
              'p'   : 'protein',
              'ncl' : 'no_command_line',
              'ptm' : 'post_trans_mod'},

             {'charge'          : 1,
              'bond'            : 1,
              'Hydrogen'        : 1,
              'json'            : 1,
              'title'           : 1,
              'afvers'          : 1,
              'dialect'         : 1,
              'configuration'   : 1,
              'no_CCD2MD'       : 1,
              'names'           : 0,
              'rename'          : 0,
              'files'           : 0,
              'covalent'        : 0,
              'ligand'          : 0, 
              'userCCDPath'     : 0,
              'seeds'           : 0,
              'protein'         : 0,
              'post_trans_mod'  : 0,
              'no_command_line' : 1}]



# Position file keywords
# -----------------------

mol2_keywords = {'ID'     : 0,
                 'name'   : 1,
                 'x'      : 2,
                 'y'      : 3,
                 'z'      : 4,
                 'elem'   : 5,
                 'resi'   : 6,
                 'resnm'  : 7}
            
keyword_dict = {'pdb' : {'ID'    : [6, 11],
                         'name'  : [12, 16],
                         'resnm' : [17, 21],
                         'resi'  : [22, 26],
                         'x'     : [30, 38],
                         'y'     : [38, 46],
                         'z'     : [46, 54],
                         'elem'  : [76, 78]},
            
                'crd' : {'resnm' : [22,   30],
                         'name'  : [32,   40],
                         'x'     : [41,   60],
                         'y'     : [61,   80],
                         'z'     : [81,  100],
                         'elem'  : [101, 102]}, # Empty element identifier to force self-identification


                'gro' : {'resnm' : [5,  10],
                         'name'  : [10, 15],
                         'x'     : [20, 28],
                         'y'     : [28, 36],
                         'z'     : [36, 44],
                         'elem'  : [14, 15]}} # Empty element identifier to force self-identification - assume 4 letter code only

                                
# resnm_dict = {'pdb' : [17, 21], 'crd' : [22, 30], 'gro' : [5, 10]}


mol2bond_map       = {'1' : ['SING', 'N'], '2': ['DOUB', 'N'], '3': ['TRIP', 'N'], 'ar': ['AROM', 'Y']}
first_unknown_bond = True
first_general_bond = True
# first_covalent     = True
first_unknown      = True
first_bond         = True
first_component    = True


# ================== #
# ------------------ #
# INPUT FILE PARSERS #
# ------------------ #
# ================== #


# Read configuration file
# ------------------------

def read_configuration_file(configuration, args, flag_dict=flag_dict):
    ''' Read information from configuration file and command line to generate input data. '''

    # Read line, get the corresponding flag and check against args
    # Change args if needed

    no_command_line_extra = False if args.no_command_line is None else args.no_command_line
    
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
        flag = flag if flag_dict[0].get(flag) == None else flag_dict[0][flag]
        val  = l[1].split('#')[0].strip() # Remove trailing comments

        if flag_dict[1][flag] == 0:
            # Multiple arguments 
            val = val.split()
        else:
            if   val.lower() == 'true':
                val = True
            elif val.lower() == 'false':
                val = False
            
        
        if config_info.get(flag) != None:
            if flag_dict[1][flag] == 1:
                # Single argument 
                print('ERROR: Information for the {} flag has been set multiple times within the configuration file.'.format(flag))
                exit()
            else:
                config_info[flag].extend(val)

        else:
            config_info[flag] = val
        
    # Convert to long flags and check every argument for agreement 

    if args.no_command_line == None: 
        # no_command_line not set on the command line, need to check if set in configuration file
        if config_info.get('no_command_line') != None: 
            no_command_line_extra = config_info.get('no_command_line')
                    
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
        if len(np.where(test == key)[0]) == 0:
            # No entry
            continue


        # Optional argument
        arg_test = config_info.pop(key)

        # Update command line input if not set
        # -------------------------------------
        
        # Have the value from configuration file, now need to compare to the command line
        if flag_dict[1][key] == 0:
            # Check if additional command line arguments can be added
            if no_command_line_extra and len(args.__dict__[key])!=0:
                print('# WARNING: parameters for --{} have been set via the command line, but passing additional arguments via the command line has been disabled.'.format(key))
                args.__dict__[key] = arg_test
            else:
                # Append to the command line argument
                args.__dict__[key].extend(arg_test)
                
        elif args.__dict__[key] == None:
            # No value set, and not a flag for additional command line information
            # Can replace without issue
            args.__dict__[key] = arg_test

        else:
            # Single argument has been set on the command line and in the configuration file
            print('ERROR: The input flag --{} has been set on the command line and within the configuration file.'.format(key))
            exit()
                    
        
    # Check for entries in configuration file which don't correspond to args
    
    if len(config_info) != 0:
        print('# WARNING: some parameters set in {} are not applicable to the programme chosen. Please ensure that you have specified the correct conversion.')
        print('# WARNING: these parameters are {}'.format(', '.join(config_info.keys())))
    
    return args

    


# Fasta file parsing
# -------------------

def read_fasta(fastafile):
    ''' Read in and split a fasta file.'''
    
    fasta = open(fastafile, 'r').read().split('>')   # Split into different seqeunces                                 
    fasta = [line for line in fasta if len(line)!=0] # Trim newlines                                                  
    fasta = [entry.split('\n') for entry in fasta]
    fasta = [[line for line in sequence if len(line)!=0] for sequence in fasta] # Trim newlines                       
    protein = [''.join(sequence[1:]) for sequence in fasta]

    return protein

    

# Determine elements from names
# -----------------------------
    
def get_elements(name):
    ''' Use atom name to determine element. '''
    if name.count('H')!=0:
        return 'H'
    elif name.count('O')!=0:
        return 'O'
    elif name.count('N')!=0:
        return 'N'
    elif name.count('C')!=0:
        return 'C'
    elif name.count('S')!=0:
        return 'S'
    elif name.count('P')!=0:
        return 'P'



# Split list into components
# --------------------------

def get_components(component):
    ''' Split a list into unique consituent parts. '''
    
    components = []
    in_range = True ; i = 0
    
    while in_range:
        try:
            num = int(component[i+1])
            if '.' in component[i]: # Should be irrelevant for ligands => assume FASTA wlog
                comp_i = read_fasta(component[i])
            else:
                comp_i = [component[i]]
            components.extend(list(comp_i)*num)
            i += 2
        except (ValueError, IndexError) as e:
            if '.' in component[i]:
                comp_i = read_fasta(component[i])
            else:
                comp_i = [component[i]]
                components.extend(list(comp_i))
            i += 1
        in_range = False if i >= len(component) else True        

    return components


    
# =============================== #
# ------------------------------- #
# USER CCD INFORMATION GENERATION #
# ------------------------------- #
# =============================== #


# Charge function
# ----------------

def get_charge(charge, charge_cutoff):
    ''' Convert float charge to integer charge based on cutoff. '''
    if abs(charge) < charge_cutoff:
        return 0
    # Some charge - want charge of 2 to require > 1 + charge_cutoff
    int_charge = np.floor(abs(charge))
    rem_charge = abs(charge) - int_charge
    tot_charge = int_charge + 1 if rem_charge > charge_cutoff else int_charge
    tot_charge = tot_charge if charge > 0 else -1* tot_charge
    return tot_charge


    
# Bond function
# --------------

def get_bonds_dist(atom_data, bond_cutoff):
    ''' Determine if two atoms are bonded based on distance cutoff. '''
    bonds = []
    for i, atomi in enumerate(atom_data[:-1]):
        posx = np.array([atomi['x'], atomi['y'], atomi['z']])
        for j, atomj in enumerate(atom_data[i+1:]):
            posy = np.array([atomj['x'], atomj['y'], atomj['z']])
            dist = np.sqrt(np.vdot((posx-posy), (posx-posy)))
            if dist < bond_cutoff:
                bonds.append([atomi['name'], atomj['name'], 'SING', 'N'])
    return bonds



# Split files based on containing position or bonding information
# ----------------------------------------------------------------

def get_files(all_files):
    ''' Split files into containing position and bonding information. '''

    position_files = [] ;  bonding_files  = []
    
    first_rtp = True ; first_rtf = True ; first_mol = True ; first_itp = True ; first_gro = True

    for fl in all_files:
        if fl[-3:] == 'pdb':
            position_files.append(fl)
        elif fl[-3:] == 'crd':
            position_files.append(fl)
        elif fl[-3:] == 'gro':
            position_files.append(fl)
            if first_gro:
                print('# WARNING: gro files cannot be used for userCCD codes longer than 4 characters.')
                first_gro = False
        elif fl[-4:] == 'mol2':
            position_files.append(fl)
            if first_mol:
                print('# WARNING: if mol2 files do not contain NAME within @<TRIPOS>ATOM this must be given as the molecule name in @<TRIPOS>MOLECULE')
                print('# WARNING: there may be issues with determining bonding if ligands to be converted from mol2 files also appear in pdb/crd/gro files')
        elif fl[-3:] == 'rtp':
            bonding_files.append(fl)
            if first_rtp:
                print('# WARNING: rtp file entries must start with [ NAME ], then [ atoms ] and [ bonds ]')
                first_rtp = False
        elif fl[-3:] == 'rtf':
            bonding_files.append(fl)
            if first_rtf:
                print('# WARNING: rtf files must contain only a single ligand')
                first_rtf = False
        elif fl[-3:] == 'itp':
            bonding_files.append(fl)
            if first_itp:
                print('# WARNING: itp files must contain only a single ligand')
                first_itp = False
        else:
            print('# WARNING: {} was input but is not a .crd/.gro/.itp/.mol2/.pdb/.rtf/.rtp file'.format(fl))
            
    return position_files, bonding_files



# ============================= #
# ----------------------------- #
# POSITION/BONDING FILE PARSERS #
# ----------------------------- #
# ============================= #


# Get position information from mol2 file
# ----------------------------------------

def read_mol2(mol2_data, name, keywords=mol2_keywords):
    ''' Get information for a single ligand from a mol2 file. '''
    
    # Read mol2 file - consider both full-length and short atom information
    # ----------------------------------------------------------------------
    
    # Test if atomic or molecular data contains residue name - avoid similar names
    atomic = mol2_data.split('@<TRIPOS>ATOM')[1].split('@')[0].split('\n')
    molec  = mol2_data.split('@<TRIPOS>MOLECULE')[1].split('@')[0].split('\n')

    atomic = [[p for p in line.split(' ') if len(p)!=0] for line in atomic if len(line)!=0]
    molec  = [dat for dat in molec  if len(dat)!=0]

    mbonds = mol2_data.split('@<TRIPOS>BOND')[1].split('@')[0].split('\n')
    mbonds = [[b for b in line.split(' ') if len(b)!=0] for line in mbonds if len(line)!=0]

    if len(atomic[0]) >= 8:
        # Optional keyword of name present - scan through molecules
        atoms = np.array(atomic)
        ligs  = ' '.join(atoms[:, keywords['resnm']])
        if ligs.count(' {} '.format(name))==0 and ligs.count(' {} '.format(name.lower()))==0:
            # Reject similar names
            return None, False, None
        else:
            descript = molec[0].strip().strip(';').strip() if molec[0].strip() != '[ atoms ]' else '?'
            all_ats  = False # Check descriptor of atom before adding
            if len(atomic[0]) >= 9:
                keywords['charge'] = 8 # Optional keyword
    else:
        # Check molecular descriptor
        if molec[0].strip()==name:
            if len(molec) == 6:
                descript = descript if len(molec[5].strip())==0 else molec[5].strip()
            all_ats  = True # Consider all atoms
        else:
            return None, False, None
            
    molmap   = {}               # Map the name of the atom to its number for bonding information
    pos_info = []
    
    # Loop through atoms looking for constituents of molecule
    # --------------------------------------------------------

    for atom in atomic:
        if not all_ats:
            if atom[keywords['resnm']].strip() != name:
                # Reject similar names
                continue
        # Append relevant information
        pos_info.append({})
        pos_info[-1]['charge'] = 0
        for i, key in enumerate(keywords.keys()):
            sect = atom[keywords[key]]
            if key == 'name':
                pos_info[-1][key] = sect.strip()
                molmap[atom[0]] = pos_info[-1][key]
            elif key == 'elem':
                pos_info[-1][key] = sect.strip().split('.')[0]
            elif key != 'resnm':
                pos_info[-1][key] = float(sect.strip())            

    bonds = []
    
    # Read mol2 file for bonding
    # ---------------------------
                
    for line in mbonds:
        # Check for correct atoms in bond
        try:
            currbond = [molmap[line[1]], molmap[line[2]]]
        except KeyError:
            continue

        try:
            currbond.extend(mol2bond_map[line[3]])
        except KeyError:
            global first_unknown_bond
            if first_unknown_bond:
                print('# WARNING: Allocating unknown bonds as single bonds - this should not affect output.')
                first_unknown_bond = False
            currbond.extend(mol2bond_map['1'])
        bonds.append(currbond)

    return pos_info, bonds, descript



# Read positions if not mol2
# ---------------------------

def read_position(data, name, keywords):
    ''' Get positions from pdb/crd/gro file. '''
    
    data   = data.split('END')[0]
    # Consider if a similar name has been used - formatting depends on name length and file type
    
    pos           = data.split('\n')[:-1]
    element_names = True
    
    pos_info      = []

    for line in pos:
        # Only look for named molecules
        if line[0:3] == 'TER':
            continue
        if line[keywords['resnm'][0]:keywords['resnm'][1]].strip() == name:
            # Append relevant information
            pos_info.append({})
            for i, key in enumerate(keywords.keys()):
                sect = line[keywords[key][0]:keywords[key][1]]
                if key == 'name' or key=='resnm':
                    pos_info[-1][key] = sect.strip()
                elif key == 'elem':
                    pos_info[-1][key] = sect.strip() if len(sect) != 0 else sect
                    if len(pos_info[-1]['elem']) == 0:
                        element_names = False
                        pos_info[-1]['elem'] = get_elements(pos_info[-1]['name'])
                else:
                    pos_info[-1][key] = float(sect.strip())

    return pos_info, element_names



# Read bonds from rtp
# ---------------------

def read_rtp_bonding(data, name, charge):
    ''' Get descriptor, bonds and charges from rtp file. '''
    
    # Open and read rtp file
    # ----------------------

    if data.count('[ {} ]'.format(name)) == 0:
        # Check for similar names
        return None
                    
                   
    RTP      = data.split('[ {} ]'.format(name))[1]      # Strip previous molecules - consider presence of similar names
    RTP      = RTP.split(']')[:3]                        # Select atoms and bonds
    descript = RTP[0].split('\n')[1].strip(';').strip()  # Name may be below the CHARMM code
    
    # Get charges
    atom_data = RTP[1].split('\n')[1:-1]

    charges = {}
    bonds   = []
    
    for line in atom_data:
        if len(line)==0:
            continue
        info = [split for split in line.strip().split(' ') if len(split)!=0]
        if len(info)==0 or info[0][0] == ';':
            # Skip comments
            continue
        charges[info[0].strip()] = get_charge(float(info[2].strip()), charge)
        
    # Get bonds
    bond_data = RTP[2].split('\n')[1:-1]
    for line in bond_data:
        if len(line)==0:
            continue
        info = [split for split in line.strip().split(' ') if len(split)!=0]
        if len(info)==0 or info[0][0] == ';':
            # Skip comments
            continue
        bonds.append([info[0], info[1], 'SING', 'N'])

    return descript, bonds, charges
        


# Read bonds from rtf
# --------------------

def read_rtf_bonding(data, name, charge):
    ''' Get bonds and charges from rtf file. '''
    
    # Open and read rtf file
    # ----------------------
    if data.count(' {} '.format(name)) == 0 and data.count(' {} '.format(name.lower())) == 0:
        # Check for similar names
        return None
    
    RTF = data.split('\n')
    
    charges = {}
    bonds   = []
    
    for line in RTF:
        if line[:4] == 'ATOM':
            # Get charges
            info = [split for split in line.strip().split(' ') if len(split)!=0]
            charges[info[0]] = get_charge(float(info[3]), charge)
            
        elif line[:4] == 'BOND':
            info = [split for split in line.strip().split(' ') if len(split)!=0]
            bonds.append([info[1], info[2], 'SING', 'N'])
            
    return '', bonds, charges



# Read bonds from itp
# ---------------------

def read_itp_bonding(data, name, charge):
    ''' Get bonds and charges from itp file. '''

    # Open and read itp file
    # ----------------------
    if data.count(' {} '.format(name)) == 0 and data.count(' {} '.format(name.lower())) == 0 :
        # Check for similar names
        return None
    
    # Generate mapping
    # ----------------
    molmap  = {}
    charges = {}
    bonds   = []
    itpmap  = data.split('[ atoms ]')[1].split('[')[0].split('\n')
    
    atoms = [[b for b in line.split(' ') if len(b)!=0] for line in itpmap if len(line)!=0]
    for atom in atoms:
        # Get charges and mapping
        if atom[0][0] == ';':
            # Skip comments
            continue
        molmap[atom[0]] = atom[4]
        charges[atom[4]] = get_charge(float(atom[6]), charge)
    
    # Get bonds
    # ---------
    ITP = data.split('[ bonds ]')[1].split('[')[0].split('\n')
    ITP = [[b for b in line.split(' ') if len(b)!=0] for line in ITP if len(line)!=0]
    
    for line in ITP:
        if line[0][0] == ';':
            # Skip comments
            continue
        bonds.append([molmap[line[0]], molmap[line[1]], 'SING', 'N'])

    return '', bonds, charges



# Get full information about covalently bonded molecule from mol2file
# --------------------------------------------------------------------

def read_covalent_mol2(molfile, keywords=mol2_keywords):
    ''' Get position and bonding information about covalently bonded molecule from mol2file. '''
    
    full_bonds = []
    
    data    = open(molfile).read()                                                      # Open and read mol2 file - must have substructure information 
    ligname = data.split('@<TRIPOS>MOLECULE')[1].split('@')[0].split('\n')[1].strip()   # Get full name of ligand from @<TRIPOS>MOLECULE

    data = data.replace('\t', '    ')
    
    # Get atomic positions and unique component names within molecule
    # ----------------------------------------------------------------
            
    atomic = data.split('@<TRIPOS>ATOM')[1].split('@')[0].split('\n')
    atomic = [[p for p in line.split(' ') if len(p)!=0] for line in atomic if len(line)!=0]
    if np.shape(atomic)[1] == 9:
        keywords[8] = 'charge'
    elif np.shape(atomic)[1] == 10:
        keywords[8] = 'charge' ; keywords[9] = 'status'
    atomic = pd.DataFrame(atomic, columns = keywords.keys()) ; atomic = atomic.set_index('ID')
    atomic['elem'] = atomic['elem'].apply(lambda x: x.split('.')[0])      
    
    components = atomic['resnm'].unique()
    
    # Get full bonding information - intra- and inter-component
    # ----------------------------------------------------------
    mbonds = data.split('@<TRIPOS>BOND')[1].split('@')[0].split('\n')
    mbonds = [[b for b in line.split(' ') if len(b)!=0] for line in mbonds if len(line)!=0]
    
    for line in mbonds:
        # Need to take into account both name and residue number 
        currbond = [atomic.at[line[1], 'name'], atomic.at[line[1], 'resi'],
                    atomic.at[line[2], 'name'], atomic.at[line[2], 'resi']]
        try:
            currbond.extend(mol2bond_map[line[3]])
        except KeyError:
            global first_unknown
            if first_unknown:
                print('# WARNING: Allocating unknown bonds as single bonds - this should not affect output.')
                first_unknown = False
            currbond.extend(mol2bond_map['1'])
        full_bonds.append(currbond)

    return atomic, full_bonds, components


# Get covalent position information
# ----------------------------------

def get_covalent_position(molfile, keywords=keyword_dict['pdb']):
    ''' Get position information for covalently bonded molecule - NOTE: currently only pdb file. '''

    posdata  = open(molfile).read()
    posdata  = posdata.split('END')[0]
    posdata  = posdata.split('\n')[:-1]
    pos_info = []
    
    element_names = True
            
    for line in posdata:
        if line.count('TER') == 1:
            # Termination of residue chain, skip
            continue
        if line.count('ATOM')==0 and line.count('HETATM')==0:
            # Not residue information
            continue
        pos_info.append({})
        for i, key in enumerate(keywords.keys()):
            sect = line[keywords[key][0]:keywords[key][1]]
            if key =='x' or key == 'y' or key == 'z':
                pos_info[-1][key] = float(sect.strip())
            elif key == 'elem':
                pos_info[-1][key] = sect.strip()
                if len(pos_info[-1]['elem']) == 0:
                    if element_names:
                        print('# WARNING: Element names are missing - attempting to infer from atom names. Note that this may cause issues.')
                        element_names = False
                    pos_info[-1]['elem'] = get_elements(pos_info[-1]['name'])
                else:
                    pos_info[-1][key] = sect.strip()
    
    atomic     = pd.DataFrame.from_dict(pos_info, orient='columns') ; atomic = atomic.set_index('ID')
    components = atomic['resnm'].unique()

    return atomic, components



# Get covalent bonding information
# ---------------------------------

def get_covalent_bonding_itp(bondfl, atomic):
    ''' Get bonding information for covalently bonded molecule - NOTE: currently only itp file. '''

    bondfile = open(bondfl).read()
    bondfile = bondfile.replace('\t', '    ')
            
    # Get full name of ligand from [ moleculetype ]  
    # ---------------------------------------------
    ligname = bondfile.split('[ moleculetype ]')[1].split('[')[0].split('\n')
    ligname = [l for l in ligname if len(l)!=0 and l[0:4].count(';')==0]
    ligname = ligname[0].split(' ')[0]
            
    # Generate mapping
    # ----------------
    molmap = {}
    itpmap = bondfile.split('[ atoms ]')[1].split('[')[0].split('\n')
            
    atoms = [[b for b in line.split(' ') if len(b)!=0] for line in itpmap if len(line)!=0]
    for atom in atoms:
        # Get charges and mapping
        if atom[0][0] == ';':
            # Skip comments
            continue
        molmap[atom[0]] = atom[4]
        atomic.loc[atomic['name']==atom[4], 'charge'] = get_charge(float(atom[6]), args.charge)
    
    # Get bonds
    # ---------
    ITP = bondfile.split('[ bonds ]')[1].split('[')[0].split('\n')
    ITP = [[b for b in line.split(' ') if len(b)!=0] for line in ITP if len(line)!=0]
            
    for line in ITP:
        if line[0][0] == ';':
            # Skip comments
            continue
        currbond = [atomic.at[line[0], 'name'], atomic.at[line[0], 'resi'],
                    atomic.at[line[1], 'name'], atomic.at[line[1], 'resi'],
                    'SING', 'N']
        full_bonds.append(currbond)

        return atomic, full_bonds
        
    
    
# ============================== #
# ------------------------------ #
# GENERAL CONVERSION INFORMATION #
# ------------------------------ #
# ============================== #


# Get CCD code from position_file
# --------------------------------

def CCD_from_name(name, position_files, bonding_files, bond_cutoff, charge, include_H=False):
    ''' Generate userCCD information for NAME using position and optionally bonding information. '''
        
    # Get position information
    # -------------------------

    position, bonds, descript = get_position(name, position_files)
    position_data             = pd.DataFrame.from_dict(position, orient='columns')

    # Attmept to get bonding information from files (excluding mol2)
    # ---------------------------------------------------------------

    if bonds==None:
        descript, bonds, charges = get_bonds(name, bonding_files, charge)
        if charges!=None and len(charges)!=0:
            for atom in charges.keys():
                position_data.loc[position_data['name']==atom, 'charge'] = charges[atom]

    
    if bonds == None:
    
        # Generate bonding information from position if needed
        # ------------------------------------------------------
    
        # Specify bonds from pos file - either no RTP/RTF file or wrong.
    
        # Note that AF README states that the bond order and aromacity don't matter
        # Currently not considering but that may change in future
    
        print('# WARNING: Inferring bonding from proximity, ignoring Hs.')
        print('# WARNING: Giving all atoms a charge of 0.')

        global first_bond
        if first_bond:
            print('# WARNING: Allocating all bonds as single bonds - this should not affect output.')
            first_bond = False
    
        position_data.loc[:, 'charge'] = 0
        noHdata = position_data.loc[position_data['elem'] != 'H'].to_dict(orient='records')
    
        bonds = get_bonds_dist(noHdata, bond_cutoff)
    
    if include_H:
        posdata = position_data.to_dict(orient='records')
        Hs      = []
    else:
        posdata = position_data.loc[position_data['elem'] != 'H'].to_dict(orient='records')
        Hs      = position_data.loc[position_data['elem'] == 'H', 'name'].to_list()

    return descript, posdata, bonds, Hs
    



# Get position information
# -------------------------

def get_position(name, position_files):
    ''' Scan through position files to get position information. '''
    
    # Scan through position files
    # ----------------------------

    bonds    = None
    descript = '?'
    
    for posfl in position_files:
        pos = open(posfl).read()
        
        if pos.count(name) == 0 and pos.count(name.lower())==0:
            # No reference to the ligand, move on to next position file
            continue

        if posfl[-3:] != 'pdb':
            # Don't replace in pdb where column index matters
            pos = pos.replace('\t', '    ')
    
        if posfl[-4:] == 'mol2':
            positions, bonds, descript = read_mol2(pos, name)
            if positions == None:
                continue

            mol2info ='# INFO: Gathering position '
            if len(bonds)!=0:
                mol2info += 'and bonding '
                
            print(mol2info + 'information from {}.'.format(posfl))

            return positions, bonds, descript
        
            
        else:
            keywords                 = keyword_dict[posfl[-3:]]
            positions, element_names = read_position(pos, name, keywords)

            if positions == None:
                continue

            print('# INFO: Gathering position information from {}. Note it is assumed ligand information is only present as a single ligand in one position file.'.format(posfl))
            if not element_names:
                print('# WARNING: Element names are missing - attempting to infer from atom names. Note that this may cause issues.')

            return positions, bonds, descript

    # No reference to ligand has been found
    print('ERROR: {} has not been found in {}'.format(name, position_files))
    exit()      




# Get bonding information
# ------------------------
    
def get_bonds(name, bonding_files, charge):
    ''' Scan through bonding files to get bonding information. '''
    
    # Scan through RTP/RTF/ITP files
    # -------------------------------

    # Ignore mol2files as these should be picked up by position generation
    
    for bndfl in bonding_files:

        bonds = None

        global first_general_bond
        if first_general_bond:
            print('# WARNING: Allocating all bonds as single bonds - this should not affect output.')
            first_general_bond = False

        data = open(bndfl).read()
        if data.count(name) == 0 and data.count(name.lower()) == 0:
            # Note crd file from CHARMMGUI may convert to lowercase
            continue

        data = data.replace('\t', '    ')
        
        if bndfl[-3:]== 'rtp':
            descript, bonds, charges = read_rtp_bonding(data, name, charge)
    
        elif bndfl[-3:] == 'rtf':
            descript, bonds, charges = read_rtf_bonding(data, name, charge)

        elif bndfl[-3:] == 'itp':
            descript, bonds, charges = read_itp_bonding(data, name, charge)

        if bonds != None:
            # Do returns etc.
            
             print('# INFO: Gathering bonding and charge information from {}.'.format(bndfl))

             return descript, bonds, charges


    # No bonding information found
    # ----------------------------

    return None, None, None

         
# =================== #
# ------------------- #
# OUTPUT FILE WRITERS #
# ------------------- #
# =================== #


# Write cif information
# ----------------------

def cif_information(nname, descript, posdata):
    ''' Write mmcif field information for new molecule (used in userCCD). '''
    
    cif_content = []
    cif_content.append("data_"+nname+"\n#")
    cif_content.append("_chem_comp.id "+nname)
    cif_content.append("_chem_comp.name '{}'".format(descript))
    cif_content.append("_chem_comp.type lipid")
    cif_content.append("_chem_comp.formula ?")
    cif_content.append("_chem_comp.mon_nstd_parent_comp_id ?")
    cif_content.append("_chem_comp.pdbx_synonyms ?")
    cif_content.append("_chem_comp.formula_weight ?")
    cif_content.append("#")
    cif_content.append("loop_")
    cif_content.append("_chem_comp_atom.comp_id")
    cif_content.append("_chem_comp_atom.atom_id")
    cif_content.append("_chem_comp_atom.type_symbol")
    cif_content.append("_chem_comp_atom.charge")
    cif_content.append("_chem_comp_atom.pdbx_leaving_atom_flag")
    cif_content.append("_chem_comp_atom.pdbx_model_Cartn_x_ideal")
    cif_content.append("_chem_comp_atom.pdbx_model_Cartn_y_ideal")
    cif_content.append("_chem_comp_atom.pdbx_model_Cartn_z_ideal")
    
    for atom in posdata:
        cif_content.append(f"{nname} {atom['name']} {atom['elem']} {int(atom['charge'])} N {atom['x']:.3f} {atom['y']:.3f} {atom['z']:.3f}")        
        
    cif_content.append("#")
    cif_content.append("loop_")
    cif_content.append("_chem_comp_bond.atom_id_1")
    cif_content.append("_chem_comp_bond.atom_id_2")
    cif_content.append("_chem_comp_bond.value_order")
    cif_content.append("_chem_comp_bond.pdbx_aromatic_flag")
    
    return cif_content



# Write cif information for covalently bonded molecule
# -----------------------------------------------------

def covalent_cif(ligname, atomic, bonds, components, Hydrogens):
    ''' Write mmcif field information for new covalently-bonded molecule (used in userCCD). '''

    chains = []
    
    molecule_file = open(ligname+'_output.json', 'w')
    molecule_file.write('{\n\t"userCCD" : "')    
    
    resimapping = {}
    for counter, resi in enumerate(list(atomic['resi'].unique())):
        chains.append([counter+1, 'A', atomic.loc[atomic['resi']==resi, 'resnm'][0]])
        # AF3 indexes from 0, add residues to the same chain
        resimapping[resi] = str(counter+1)
    
    for j, component in enumerate(components):
        # Determine if multiple components are the same
        resis = list(atomic.loc[atomic['resnm']==component, 'resi'].unique())
        global first_component
        if len(resis) > 1 and first_component:
            print('# INFO: Where multiple of the same component are present the first component present in the file is used to construct userCCD')
            first_component = False
            
        component_data        = atomic.loc[atomic['resi']==resis[0]]
        intra_component_bonds = bonds.loc[bonds['resi1']==resis[0]].loc[bonds['resi2']==resis[0]]
                    
        # Write userCCD data for component only 
        # --------------------------------------
                    
        if Hydrogens:
            posdata = component_data.to_dict(orient='records')
            Hs      = []
        else:
            posdata = component_data.loc[component_data['elem'] != 'H'].to_dict(orient='records')
            Hs      = component_data.loc[component_data['elem'] == 'H', 'name'].to_list()
    
        cif_content = cif_information(component.upper(), '?',  posdata)
                
        for row, bond in intra_component_bonds.iterrows():
            if Hs.count(bond['atom1'])==0 and Hs.count(bond['atom2'])==0:
                cif_content.append(f"{bond['atom1']} {bond['atom2']} {bond['type']} {bond['ar']}")
        cif_content.append('#')
            
        # Write to CIF file
        molecule_file.write("\n".join(cif_content))

        # =================== #
        # Write modifications # 
        # =================== #
    
        # Close userCCD section
        
        molecule_file.write('"\n  "sequences": [ \n')
    
        chains = pd.DataFrame(chains, columns=['resi', 'chain', 'resnm']) ; chains = chains.set_index('resi')
        
        # Add correct number of ligands - want to keep a record of this for JSON
        # ----------------------------------------------------------------------
        lig = '   {"ligand": {"id": ["A"], "ccdCodes": ["'+'", "'.join(chains['resnm'])+'"]}}'
    
        molecule_file.write(lig+', \n')
            
        # Add modifications - currently won't work for proteins but could extend in future
        # --------------------------------------------------------------------------------
    
        modification_info = []
        molecule_file.write('\n    ], \n   "bondedAtomPairs": [\n')
        for row, bond in flbnds.loc[flbnds['resi1']!=flbnds['resi2']].iterrows():
            # Assume all ligands at residue 1
            mod  = '    [["A", '+resimapping[bond['resi1']]+', "'+bond['atom1']+'"],'
            mod +=      '["A", '+resimapping[bond['resi2']]+', "'+bond['atom2']+'"]]'
            modification_info.append(mod)

        molecule_file.write(',\n'.join(modification_info))
        molecule_file.write('],')
            
        molecule_file.close()

        return None


    
# Extract cif information from covalently bonded molecule
# --------------------------------------------------------

def extract_covalent(currlig, lig, json, inCIF, CIF, title, lig_chains, this_session=False):
    ''' Generate information for json file from cif file containing covalent ligand information. '''

    # Covalent ligand - find components
    # ---------------------------------
                
    CCDs       = currlig.split('"userCCD" : "')[1].split('"')[0]
    components = CCDs.split('data_')[1:]
    compnames  = [comp.split('\n')[0] for comp in components]
    
    # Check against other defined ligands
    # ------------------------------------
    
    for i, compname in enumerate(compnames):
        if inCIF.get(compname) == None:
            inCIF[compname] = '{}_output.cif'.format(lig)
            inCIF[compname] = inCIF[compname] + ' (this session)' if this_session else inCIF[compname]
            CIF.write('data_'+components[i])
        else:
            print('# WARNING: {} is defined in {}_output.cif but {}_CCD.cif contains a version from {}'.format(compname, lig, title, inCIF[compname])) 
                        
    # Find covalent information
    # --------------------------

    lig_info  = currlig.split('"id": ["A"],')[1].split('}}')[0]
    json.write('{"ligand": {"id": ["'+', '.join(lig_chains)+'"],'+lig_info+'}}')    
                
    bonding_info = currlig.split('"bondedAtomPairs": [')[1]
    bonding_info = bonding_info[:-2]

    bonded = []
    for chain in lig_chains:
        # Rename the chains from "A" to what is appropriate
        bonded.append(bonding_info)
        bonded[-1].replace('["A"', '["{}"'.format(lig_chain))

    return json, CIF, inCIF, bonded



# Extract ligand information
# ---------------------------

def extract_ligand(currlig, lig_check, allchains, inCIF, CIF, json, title, lignums, ligfile):
    ''' Generate information for json file from cif file containing ligand information. '''

    # Find ligands in file
    # ----------------------
                
    ligands  = currlig.split('data_')[1:]
    lignames = [lgs.split('\n')[0] for lgs in ligands]
    
    # Check against other defined ligands
    # ------------------------------------

    for i, ligname in enumerate(lignames):
        if lig_check.get(ligname) != None:
            if not lig_check[ligname] or inCIF.get(ligname) != None:
                print('# WARNING: {} is defined in {} but {}_CCD.cif contains a version from {}'.format(ligname, ligfile, title, inCIF[ligname]))
            # Want to add this ligand
            elif inCIF.get(ligname) == None:
                inCIF[ligname] = ligfile
                CIF.write('data_'+ligands[i])

            if lig_check[ligname]:
                
                numligs = lignums.pop(ligname)
                
                lig_chains = allchains[:numligs]
                allchains  = allchains[numligs:]
            
                lig_check[ligname] = False
                if numligs!=0:
                    json.write('{"ligand": {"id": ["'+', '.join(lig_chains)+'"], "ccdCodes": ["{}"]}}'.format(ligname))
                if len(lignums) != 0:
                    json.write(',\n\t\t')                

    return json, lig_check, lignums, allchains, CIF, inCIF
