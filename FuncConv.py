#############################################################
#                                                            #
#    Functions to convert between CCD and CHARMM ordering    #
#                                                            #
##############################################################

# User-dependent variables
# --------------------------

CCD2MD_dir = '/home/user/Packages/CCD2MD/'  # Location of CCD2MD - to be set by the user
num_CPUs   = 1 # Number of CPUs to use for membrane insertion. Integer only.

'''
General file for functions to convert between CCD/CHARMM/Martini files

Last Update: K Blow 06/10/25

Contains:
read_CIF(INNAME)   # Read CIF file INNAME, return system data (note, only data relevant to PDBs is retained)
read_PBD(INNAME)   # Read PDB file INNAME, return system data and optionally title and crystallographic information
read_GRO(INNAME)   # Read GRO file INNAME, return system data and optionally title and crystallographic information

determine_element(NAME) # Determine element type from atom name

get_residues(system_data, data_type, SMILES=[], ligand_chain=False, database=df) # Get information about the residues in system_data which need to be converted (data_type gives the current data type. Add SMILES keys and determine if ligands should be separated into chains.
get_reference(NAME, database=df, Ref_loc='Ref_data/')                            # Get reference information for residue NAME. Prevent displaying information multiple times for the same moleucule, although this may need refinement.

convert_atomistic(name, inputdict, inputtype, ligand_chains=False, database=df, Ref_loc='Ref_data/') # Change ordering between atomistic orders - CCD (including SMILES strings where accepted) and CHARMM

get_CG_params(command_line, martinize, elastic, go)                                       # Determine additional flags to be passed to martinize2
to_CG(inputfile, outputfile, input_data, martinizeparams, ligands, system_data, types, convresi, prot,
          database=df, Ref_loc=Ref_data, mass=masses, newlipidome=False                   # Take in an atomic list of locations, and return a CG representation.
prot_to_CG(inputfile, basename, input_data, martinizeparams, ligands, types, database=df) # Convert protein to CG using martinize2 on the command line. 

write_PDB_atom_line(OUTFILE, counter, data)                                # For a given OUTFILE, write data corresponding to a given atom
write_PDB(OUTFILE, ordered_dict, title=[], cryst=[], ligchain = False) # For a given OUTFILE, write all of the data for the system (ordered_dict)

get_command_line_parameters(command_line, flags) # Determine additional command line parameters to pass to various commands.

build_membrane_CG(ligands, CG_output, outputfile, command_line, mempro_additional, membrane_comp, ion_conc, database=df, newlipidome=False)                                          # Builds a CG membrane around the system using MemPrO
convert_membrane_at(system_data, basename, command_line, CG2AT_additional) # Converts membrane using CG2AT

get_topology_CG(outputfile, membrane, ligands, prot, inputfile, newlipidome=False, database=df) # Write topology for CG system
get_topology_atomistic(outputfile, membrane, at_command=None, output_data=None)                 # Write topology for atomistic system

check_residue_number(ordered_dict) # Renumber residues where this exceeds max PDB can handle
convert_vectors(box_vecs)          # Convert box information from GRO format to PDB format
'''

import numpy as np
import pandas as pd
import warnings
import subprocess, os, shutil, sys
import re, glob
import ast

# Locations of packages
# ---------------------

Ref_data   = CCD2MD_dir + 'Ref_data/'
MProPath   = CCD2MD_dir + 'MemPrO/MemPrO_Script.py'
CHARMMPath = CCD2MD_dir + 'charmm36-ccd2md.ff/'
CG2ATPath  = CCD2MD_dir + 'cg2at-ccd2md/cg2at'
oldmartini = CCD2MD_dir + 'martini_v3.itp'
newmartini = CCD2MD_dir + 'martini_new_lipidome.itp'
oldinsane  = CCD2MD_dir + 'MemPrO/Insane4MemPrO.py'
newinsane  = CCD2MD_dir + 'MemPrO/Insane4MemPrO_new_lipidome.py'

os.environ['NUM_CPU'] = str(num_CPUs) # This may be best set by the user

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
               '_chem_comp_atom.pdbx_model_Cartn_y_ideal'  : 'y',     # y coordinate
                '_chem_comp_atom.pdbx_model_Cartn_z_ideal' : 'z',     # z coordinate
                '_chem_comp_atom.type_symbol'              : 'elem'  # Element symbol
                }

PDB_keywords = {'entry' : [0, 6],
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
all_flags = ['-S',   '--SMILES',   '-M',  '--martinize', '-G', '--go',  '-E',  '--elastic', '-lc', '--ligchain',
             '-mem', '--membrane', '-mp', '--mempro',    '-C', '--conc', '-gh', '--pdb2gmx', '-V', '--Version',
             '-at', '--cg2at', '-nl', '--newlipidome']

df = pd.read_csv(Ref_data+"database.csv", index_col='Name', skipinitialspace=True)

lig_names = []

base_ptms = ['CYST', 'CYSD', 'CYSP', 'CYSG', 'CYSF', 'GLYM']
PTMs = set(base_ptms + [ptm + '_user' for ptm in base_ptms])
terminal_PTMs = ['CYST', 'GLYM']

__version__ = "0.0.2"

def read_CIF(name):
    ''' Read CIF file into molecule dictionary. '''

    fle = open(name, 'r')
    CIF = fle.read().split('loop_')

    out_info = [{}]

    print('# INFO: Reading CIF file {}. Note that this assumes all fields are consistent in CIF ordering'.format(name))
    
    title = ['TITLE     Converted from {} using CCD2MD.'.format(name)]
    
    CCD_direct = False
    
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
                    print('# WARNING: This file has likely come direct from the CCD and therefore may have possible alternative labels and missing data. This may cause issues.')
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
                
                
        for i, line in enumerate(data):
            # Read in each line of data 
            line = line.split()
            # Assuming consistent ordering between blocks
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

    # Deal with missing data in CCD files
    if CCD_direct:
        # Check for the presence/absence of information
        missing = {'occ': 1.0, 'B': 0.0, 'entry': 'ATOM', 'chain': 'A'}
        for key in missing.keys():
            if list(out_info[0].keys()).count(key) == 0:
                out_info = [dict(atom, **{key: missing[key]}) for atom in out_info]

        if list(out_info[0].keys()).count('resi') == 0:
            print('# WARNING: Inferring residue IDs - these will be incorrect if there are two sequential residues of the same type')
            resi               = 1
            out_info[0]['resi'] = resi
            for i in range(1, len(out_info)):
                resi = resi+1 if out_info[i]['resnm'] != out_info[i-1]['resnm'] else resi
                out_info[i]['resi'] = resi

    else:
        print('# INFO: Co-folding predictions usually have accuracy encoded in B factor.')
    return out_info, title


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
        if line.count('REMARK') == 1:
            title.append(line)
        elif line.count('AT') == 0:
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
                else:
                    out_info[-1][key] = sect.strip(' ')
                                        
    return out_info, title, cryst

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
        
def get_residues(system_data, data_type, SMILES=[], ligand_chain=False, database=df):
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

        residues = []

        num_residues = last_residue+1-first_residue
        
        for i in range(first_residue, last_residue+1):
            try:
                residues.append(list(chain_data.loc[chain_data['resi']==i, 'resnm'])[0])
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
                new_suffix = re.compile('\d_1') # Check for \d_1 to determine if new r old Chai naming -
                                                # Prevents error with C_1 etc.
                first_smiles = False

            # SMILES are their own chain - test residue names for different co-folding programmes
            smiles_loc = np.where(np.array([r[:3] for r in residues])=='LIG')[0]
            if len(smiles_loc) == 0:
                # Test different name
                smiles_name = re.compile('l\d\d')
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
            
            for i in range(num_residues):
                res_data = chain_data.loc[chain_data['resi']==i+first_residue]
                for index, atom_data in res_data.iterrows():
                    curr          = system_data.iloc[index].to_dict()
                    curr['chain'] = ''
                    curr['resi']  = prev_max + i + 1

                    system_data.iloc[index] = curr
                    
            convresi.extend([['', int(i)+1+prev_max] for i in range(len(residues))])
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
                    for index, atom_data in res_data.iterrows():
                        curr          = system_data.iloc[index].to_dict()
                        curr['chain'] = curr_chain
                        curr['resi']  = curr_resi

                        system_data.iloc[index] = curr
                    convresi.extend([[curr_chain, curr_resi]])
                    if ligand_chain:
                        max_chnID += 1
            prev_max   = 1                  if ligand_chain else curr_resi
            prev_chain = chars[max_chnID-1] if ligand_chain else chain
            
        elif (chain != 'A' and prev_chain == ''):
            # Ligands first?
            
            for i in range(num_residues):
                res_data  = chain_data.loc[chain_data['resi']==i+first_residue]
                new_chain = chars[max_chnID] if ligand_chain else ''
                new_resi  = 1                if ligand_chain else prev_max + i + first_residue + 1
                for index, atom_data in res_data.iterrows():
                    curr          = system_data.iloc[index].to_dict()
                    curr['chain'] = new_chain
                    curr['resi']  = new_resi
                
                    system_data.iloc[index] = curr
                convresi.extend([[new_chain, int(i)+1+prev_max] for i in range(len(residues)) if (residues[i] in cres)])
                if ligand_chain:
                    max_chnID += 1

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

converted = []

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
        print('# INFO: No network information was provided. Defaulting to elastic network with default parameters.')
        martini.append('-elastic')


    SetFlags = {'-f' : [1, 'input file'], '-x' : [1, 'output file'], '-o' : [1, 'output file'],
                '-ignore' : [1, 'additional ligands']}
    
    martini = np.array(martini)

    # Remove any pre-set flags
    for flag in SetFlags.keys():
        if len(np.where(martini==flag)[0])!=0:
            print('# WARNING: You have tried to overwrite the {} passed to Martinize2. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))
            martini = np.hstack((martini[:np.where(martino == flag)[0][0]], martini[np.where(martini == flag)[0][0]+1+SetFlags[flag][0]:]))        
    return martini


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
        if int(version[1]) > 0 or int(version[1]) >= 14:
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

    martini.extend(martinizeparams)
    result = subprocess.run(martini)

    # Test output

    assert result.returncode==0, "ERROR: Failed to run martinize2, please check for errors in your input"
    
    return None


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


def get_command_line_parameters(command_line, flags):
    ''' Determine additional command line parameters to pass to various commands. '''

    CCD2MD_Flags = [f for f in all_flags if f not in flags]
    
    # User options for mempro - need to test, otherwise write defaults
    init = np.where(command_line==flags[0])[0]
    init = init[0] if len(init) != 0 else np.where(command_line==flags[1])[0][0]
    
    possargs = command_line[init+1:]

    # Possargs now contains possible arguments for the martinize2 command.
    # Need to strip out any possibilities from different flags

    for f in CCD2MD_Flags:
        if len(np.where(possargs==f)[0])!=0:
            # This flag is present, remove this
            possargs = possargs[:np.where(possargs==f)[0][0]]

    # Should now only have the relevant parameters
    build_sys = possargs

    return build_sys

def build_membrane_CG(ligands, CG_output, outputfile, command_line, mempro_additional, membrane_comp, ion_conc, database=df, newlipidome = False):
    ''' Build membrane using MemPrO. '''
    
    CG_ligands = [database.at[name, 'CGName'] for name in ligands]
    
    MemPrO = ['python', MProPath, '-f', CG_output, '-res', ','.join(set(CG_ligands)), '-o', '.'.join(outputfile.split('.')[:-1])+'_MemPrO']
    
    try:
        membrane_comp = list(get_command_line_parameters(command_line, ['-mem', '--membrane']))
    except IndexError:
        # at2mem doesn't require the same flag
        membrane_comp = []

    os.environ['PATH_TO_MARTINI'] = newmartini if newlipidome else oldmartini
    os.environ['PATH_TO_INSANE']  = newinsane  if newlipidome else oldinsane

    print(os.environ['PATH_TO_INSANE'])

    
    extra_bd_args = ''  # bd_args all need to be passed together
        
    if mempro_additional:
        build_sys = get_command_line_parameters(command_line, ['-mp', '--mempro'])
        SetFlags = {'-f' : [1, 'input file'], '--file_name' : [1, 'input file'],
                    '-res' : [1, 'additional CG ligands'], '--additional_residues' : [1, 'additional CG ligands'],
                    '-o' : [1, 'output directory'], '--output' : [1, 'output directory'],
                    '-bd_args' : [1], '--build_arguments': [1]}
        # Remove any pre-set flags - and separate out bd_args
        for flag in SetFlags.keys():
            if len(np.where(build_sys==flag)[0])!=0:
                if flag == '-bd_args' or '--build_arguments':
                    extra_bd_args = build_sys[np.where(build_sys == flag)[0][0]+1]
                else:
                    print('# WARNING: You have tried to overwrite the {} passed to MemPrO. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1])) 

                build_sys = np.hstack((build_sys[:np.where(build_sys == flag)[0][0]], build_sys[np.where(build_sys == flag)[0][0]+1+SetFlags[flag][0]:]))

        MemPrO.extend(list(build_sys))
                      
    extra_flags = {'-ni' : '15', '-bd': '1', '-ng' : '5', '-res_itp': os.environ['PATH_TO_MARTINI'],
                   '-bd_args': '-sol W -negi_c0 CL -posi_c0 NA'}

    extra_flags['-bd_args'] = extra_flags['-bd_args'] + ' ' + ' '.join(membrane_comp) if len(membrane_comp) != 0  \
        else extra_flags['-bd_args'] +' -l POPC -u POPC'

    extra_flags['-bd_args'] += ' -ion_conc '+str(ion_conc)+','+str(ion_conc)+','+str(ion_conc)
    
    extra_flags['-bd_args'] += ' ' + extra_bd_args
    
    for flag in extra_flags.keys():
        if MemPrO.count(flag)==0:
            MemPrO.extend([flag, extra_flags[flag]])    


    result = subprocess.run(MemPrO)

    # Test output
    assert result.returncode==0, "ERROR: Failed to run MemPrO, please check for errors in your input"
    
    return None

def convert_membrane_at(system_data, basename, command_line, CG2AT_additional):
    ''' Convert membrane from CG to atomistic via CG2AT. '''

    cg2at_args = [CG2ATPath, '-a', basename+'_nomem.pdb', '-c', basename+'_MemPrO/Rank_1/CG_System_rank_1/CG-system.gro', '-loc', basename+'_CG2AT']

    if os.path.isdir(basename+'_CG2AT/INPUT'):
        print('# WARNING: It appears that there is already a CG2AT directory of name {}. This will be deleted to allow for current conversion.'.format(basename+'_CG2AT/'))
        shutil.rmtree(basename+'_CG2AT/')
    
    if CG2AT_additional:
        cg2at_extra = get_command_line_parameters(command_line, ['-at', '--cg2at'])
        SetFlags = {'-a' : [1, 'atomistic input file'], '-c' : [1, 'CG input file'],
                    '-loc' : [1, 'output directory']}
        # Remove any pre-set flags
        for flag in SetFlags.keys():
            if len(np.where(cg2at_extra==flag)[0])!=0:
                cg2at_extra = np.hstack((cg2at_extra[:np.where(cg2at_extra == flag)[0][0]], cg2at_extra[np.where(cg2at_extra == flag)[0][0]+1+SetFlags[flag][0]:]))
                print('# WARNING: You have tried to overwrite the {} passed to CG2AT. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1])) 

        if len(np.where(cg2at_extra == '-fg')[0]) !=0:
            if cg2at_extra[np.where(cg2at_extra == '-fg')[0][0]+1] != 'martini_3-0_charmm36':
                print('# WARNING: There may be issues when not using martini 3')
        else:
            print('# WARNING: There may be issues when not using martini 3')
        
    else:
        cg2at_extra = ['-w', 'tip3p', '-fg', 'martini_3-0_charmm36']

    cg2at_args.extend(cg2at_extra)
        
    result = subprocess.run(cg2at_args)
    assert result.returncode==0, "ERROR: Failed to run CG2AT, please check for errors in your input"

    return None
    

def get_topology_CG(outputfile, membrane, ligands, prot, inputfile, newlipidome=False, database=df):
    ''' Write topology files for CG systems. '''

    if prot:
        martinize       = open('.'.join(outputfile.split('.')[:-1])+'_proteinCG.top', 'r').read()
        martinize_lines = martinize.split('\n')

        includes = [line for line in martinize_lines if len(line)!=0 and line[0] == '#']
        includes = [line for line in includes if line[:17] != '#include "martini']        
        files    = [line.split('"')[1] for line in includes]

    outputdir = '/'.join(outputfile.split('/')[:-1])
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
        subprocess.run(['sed', '-i', '/"martini.*itp"/d', loc])
        subprocess.run(['sed', '-i', '/"protein-cg.itp"/i #include "{}"'.format(copy_martini[1].split('/')[-1]), loc])
        subprocess.run(['sed', '-i', 's/#include "protein-cg.itp"/'+includes+'/g', loc])

        # Get correct number of molcules
        mols = martinize.split('[ molecules ]')[-1].split('\n')
        mols = [mol for mol in mols if len(mol)!=0]

        subprocess.run(['sed', '-i', 's/Protein\s*1/'+'\\n'.join(mols)+'/g', loc])

        # MemPrO as default does not include additional ligands - add here
        last_entry = mols[-1]
        for lig in set([l for l in ligands if l not in PTMs]):
            num = len([l for l in ligands if l==lig])
            CG_name = database.at[lig.strip(), 'CGName']
            subprocess.run(['sed', '-i', '/{}/a {} \t{}\n'.format(last_entry, CG_name, num), loc])
            last_entry = lig
        
        print('# INFO: Topology file has been created based on martinize2 and MemPrO outputs.')

    else:
        # Create topology from martinize2 output, if exists
        # --------------------------------------------------
        
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


def get_topology_atomistic(outputfile, membrane, at_command=None, output_data=None):
    ''' Write topology for atomistic systems. '''
    
    outputdir = '/'.join(outputfile.split('/')[:-1])+'/' if len(outputfile.split('/')[:-1])!=0 else './'
    
    if membrane:
        # Copy topology and itp files from output of CG2AT
        # -------------------------------------------------
        
        itpfiles = glob.glob('.'.join(outputfile.split('.')[:-1])+'_CG2AT/FINAL/*itp')
        subprocess.run(['scp', '.'.join(outputfile.split('.')[:-1])+'_CG2AT/FINAL/topol_final.top', outputdir+'topol.top'])
        for itp in itpfiles:
            subprocess.run(['scp', itp, outputdir])
        print('# INFO: Topology file generated by CG2AT')

    else:
        # Generate topology and add hydrogens using pdb2gmx
        # -------------------------------------------------
        
        subprocess.run(['scp', '-r', CHARMMPath, outputdir])        
        os.chdir(outputdir)
        outputname = outputfile.split('/')[-1]
        
        if at_command is None:
            # No additional parameters for pdb2gmx passed - use defaults
            
            gromacs = ['gmx',  'pdb2gmx',  '-f', outputname, '-o', '.'.join(outputname.split('.')[:-1])+'_H.pdb',
                       '-water', 'tip3p', '-ff', 'charmm36-ccd2md']

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

                gromacs = ['gmx',  'pdb2gmx',  '-f', outputname, '-o', '.'.join(outputname.split('.')[:-1])+'_H.pdb',
                           '-water', 'tip3p', '-ff', 'charmm36-ccd2md', '-ter']
            
            
                gmx = subprocess.Popen(gromacs, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                out = gmx.communicate(input=ter)
                print(out[1])
                
            else:
                gmx = subprocess.run(gromacs)
            
        else:
            # Generate command line input and pass to pdb2gmx
            # Command line options added after -gh or --pdb2gmx
            pdb2gmx_args = list(get_command_line_parameters(at_command, ['-gh', '--pdb2gmx']))
            
            SetFlags = {'-f' : [1, 'input file']} # Don't allow overwriting file name
            # Remove any pre-set flags
            for flag in SetFlags.keys():
                if len(np.where(pdb2gmx_args==flag)[0])!=0:
                    pdb2gmx_args = np.hstack((pdb2gmx_args[:np.where(pdb2gmx_args == flag)[0][0]], pdb2gmx_args[np.where(pdb2gmx_args == flag)[0][0]+1+SetFlags[flag][0]:]))
                    print('# WARNING: You have tried to overwrite the {} passed to pdb2gmx. This will cause an error so this command has been ignored.'.format(SetFlags[flag][1]))
                        
            gromacs = ['gmx', 'pdb2gmx', '-f', outputname]
            gromacs.extend(pdb2gmx_args)
            
            gmx = subprocess.run(gromacs)
            

        if gmx.returncode==0:
            print('# INFO: Topology file generated by pdb2gmx')
        else:
            print('# ERROR: pdb2gmx could not generate a topology file.')
    
    return None


def check_residue_number(ordered_dict):
    ''' Renumber residues where this exceeds max PDB can handle. '''
    
    dataframe = pd.DataFrame.from_dict(ordered_dict, orient='columns')
    
    if max(dataframe['resi']) > 9999:
        print('# WARNING: Number of residues exceeds 9999, which is a limitation of theP DB format. Renumbering all residues with high residue IDs from 0 (will repeat as often as needed).')
        
        while max(dataframe['resi']) > 9999:
            dataframe['resi'] = dataframe['resi'].apply(lambda x: (x if x < 10000 else x - 10000))
                
    return dataframe.to_dict('records')


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
    
