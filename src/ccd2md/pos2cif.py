##########################################
#                                         # 
#    Convert between CRD/PDB and mmCIF    # 
#                                         #
###########################################

# Import relevant functions
# -------------------------

import argparse
import pandas as pd
import numpy as np
import sys
import warnings
import string
from collections import Counter

def main():
    
    # Get command line arguments
    # ---------------------------
    
    parser = argparse.ArgumentParser(description='Generate user-defined CCD code(s) for use in AF3 from position and optionally additional bonding file(s). Saved as json file and cif file for each constituent ligand.')
    
    parser.add_argument('-v', '--version', action='version', version='Version 1.0.3')
    
    req = parser.add_argument_group('Required inputs')
    
    req.add_argument('-n', '--names',    help='CHARMM code for molecule(s) to convert, note this must not overlap with an existing CCD code.', nargs='+', default = [])
    req.add_argument('-r', '--rename',   help='Pair(s) of ligand names in format OLD NEW, where OLD is present in the position files and NEW is the desired new names.', nargs='+', default=[])
    req.add_argument('-f', '--files',    help='.gro/.pdb/.crd/.mol2 file(s) containing molecule(s) (note, can contain other moleucles). Optionally may also specify .rtp/.rtf/.itp file(s) containing molecule(s) (note, rtp files can contain other moleucles). If .rtp/.rtf files not provided, bonded information will be inferred from proximty.', nargs='+', default = [])
    req.add_argument('-c', '--covalent', help='Position and bonding information for ligands which are to be constructed using covalent modifications. Either mol2 file (including bonding information), or a pair of pdb and itp files.', nargs='+', default = [])
    
    opts = parser.add_argument_group('Optional cutoffs/arguments')
    
    opts.add_argument('-e', '--charge',   help='Minimum charge required for the charge to be non-0 (output integer charges only). Default 0.75', type=float, default=0.75)
    opts.add_argument('-b', '--bond',     help='Maximum distance (in A) between atoms to be considered bonded. Default 1.8', type=float, default=1.8)
    opts.add_argument('-H', '--Hydrogen', help='Retain hydrogens in mmCIF data. Default=False.', action='store_true')
    
    
    jsn = parser.add_argument_group('json file parameters')
    
    jsn.add_argument('-j', '--json',     help='Name of JSON file to write userCCD to. Default = "output.json"', default='output.json')
    jsn.add_argument('-t', '--title',    help='AF3 system title. Default = "pos2cif_system"', default='pos2cif_system')
    jsn.add_argument('-A', '--afvers',   help='AF3 version. Default = 2', default='2')
    jsn.add_argument('-s', '--seeds',    help='Model seeds - need not be comma separated. Default 1', default = ['1'], nargs='+')
    jsn.add_argument('-d', '--dialect',  help='Dialect. Default "alphafold3"', default="alphafold3")
    jsn.add_argument('-p', '--protein',  help='FASTA protein sequence(s) to add to system. For multiple of the same sequence (e.g. AACCS) can be "AACCS AACCS" or "AACCS 2".', default = [], nargs='+')
    
    args = parser.parse_args()
    
    position_files = []
    bonding_files  = []
    
    first_rtp = True ; first_rtf = True ; first_mol = True ; first_itp = True ; first_gro = True
    
    assert len(args.rename)%2==0, '#ERROR: rename must take in pairs of ligand names.'
    
    print('# WARNING: It is assumed ligand information is only present as a single ligand in one position file.')
    
    for fl in args.files:
        if fl[-3:] == 'pdb':
            position_files.append(fl)
        elif fl[-3:] == 'crd':
            position_files.append(fl)
        elif fl[-3:] == 'gro':
            position_files.append(fl)
            if first_gro:
                print('# WARNING: GRO files cannot be used for userCCD codes longer than 4 characters.')
                first_gro = False
        elif fl[-4:] == 'mol2':
            position_files.append(fl)
            bonding_files.append(fl)
            if first_mol:
                print('# WARNING: if mol2files do not contain NAME within @<TRIPOS>ATOM this must be given as the molecule name in @<TRIPOS>MOLECULE')
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
            print('# WARNING: {} was input but is not a .crd/.gro/.itp/.mol/.pdb/.rtf/.rtp file'.format(fl))
    
    print('# INFO: Writing userCCD to {}'.format(args.json))
    json = open(args.json, 'w')
    json.write('{\n\t"userCCD" : "')
    
    json_ligands = []
    
    # ================ #
    # Define functions #
    # ================ #
    
    # Charge function
    # ---------------
    
    def get_charge(charge, charge_cutoff):
        if abs(charge) < charge_cutoff:
            return 0
        # Some charge - want charge of 2 to require > 1 + charge_cutoff
        int_charge = np.floor(abs(charge))
        rem_charge = abs(charge) - int_charge
        tot_charge = int_charge + 1 if rem_charge > charge_cutoff else int_charge
        tot_charge = tot_charge if charge > 0 else -1* tot_charge
        return tot_charge
    
    # Distance function
    # -----------------
    
    def get_bonds_dist(atom_data):
        bonds = []
        for i, atomi in enumerate(atom_data[:-1]):
            posx = np.array([atomi['x'], atomi['y'], atomi['z']])
            for j, atomj in enumerate(atom_data[i+1:]):
                posy = np.array([atomj['x'], atomj['y'], atomj['z']])
                dist = np.sqrt(np.vdot((posx-posy), (posx-posy)))
                if dist < args.bond:
                    bonds.append([atomi['name'], atomj['name'], 'SING', 'N'])
        return bonds
    
    # cif information function
    # -------------------------
    
    def cif_information(nname, descript, posdata):
        cif_content = []
        cif_content.append("data_"+nname+"\\n#")
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

    # Fasta file function
    # -------------------

    def read_fasta(fastafile):
        ''' Read in and split a fasta file.'''

        fasta = open(fastafile, 'r').read().split('>')   # Split into different seqeunces                                 
        fasta = [line for line in fasta if len(line)!=0] # Trim newlines                                                  
        fasta = [entry.split('\n') for entry in fasta]
        fasta = [[line for line in sequence if len(line)!=0] for sequence in fasta] # Trim newlines                       
        protein = [''.join(sequence[1:]) for sequence in fasta]

        return protein

    
    mol2bond_map  = {'1' : ['SING', 'N'], '2': ['DOUB', 'N'], '3': ['TRIP', 'N'], 'ar': ['AROM', 'Y']}
    allchains     = list(string.ascii_uppercase)
    
    # Determine elements from names
    # -----------------------------
    
    def get_elements(name):
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
    
    # Define name changes
    # --------------------
    
    new_names = {} ; rename = np.array([])
    for i in range(0, len(args.rename), 2):
        rename                = np.append(rename, args.rename[i].upper())
        new_names[rename[-1]] = args.rename[i+1].upper()
    
    # ------------------------------------- #
    #                                       #
    # Work through list of names to convert #
    #                                       #
    # ------------------------------------- #
    
    args.names = [name.upper() for name in args.names]
    for i, name in enumerate(args.names):
        # Check renaming
        # --------------
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            # Note that there is a future warning about the behaviour of this comparison - may fail with new numpy release
            if len(np.where(rename==name)[0])!=0:
                # Rename
                nname = new_names[name]
                print('# INFO: Creating user-defined CCD code for {}, which has been renamed from {}'.format(nname, name))
                extra_info = ' (renamed from {})'.format(name)
            else:
                nname = name
                print('# INFO: Creating user-defined CCD code for {}'.format(name))
                extra_info = ''
                
        pos_info = []
        descript = "?"
    
        molfile = None
    
        # ======================== #
        # Get position information #
        # ======================== #
        
        # Scan through position files
        # ----------------------------
        for posfl in position_files:
            pos = open(posfl).read()
            lig_present   = False
            if pos.count(name) == 0 and pos.count(name.lower())==0:
                # No reference to the ligand, move on to next position file
                continue
    
            if posfl[-4:] == 'mol2':
                # Open and read mol2 file - consider both full-length and short atom information
                # -------------------------------------------------------------------------------
                keywords = {'name'   : 1,
                            'x'      : 2,
                            'y'      : 3,
                            'z'      : 4,
                            'elem'   : 5}
                resnm = 7
    
                # Test if atomic or molecular data contains residue name - avoid similar names
                atomic = pos.split('@<TRIPOS>ATOM')[1].split('@')[0].split('\n')
                molec  = pos.split('@<TRIPOS>MOLECULE')[1].split('@')[0].split('\n')
    
                atomic = [[p for p in line.split(' ') if len(p)!=0] for line in atomic if len(line)!=0]
                molec  = [dat for dat in molec  if len(dat)!=0]
                
                if len(atomic[0]) >= 8:
                    # Optional keyword of name present - scan through molecules
                    atoms = np.array(atomic)
                    ligs  = ' '.join(atoms[:, resnm])
                    if ligs.count(' {} '.format(name))==0 and ligs.count(' {} '.format(name.lower()))==0:
                        # Reject similar names
                        continue
                    else:
                        lig_present = True
                        print('# INFO: Gathering position and bonding information from {}.'.format(posfl))
                        descript = molec[0].strip().strip(';').strip() if molec[0].strip() != '[ atoms ]' else '?'
                        all_ats  = False # Check descriptor of atom before adding
                        if len(atomic[0]) >= 9:
                            keywords['charge'] = 8 # Optional keyword
                else:
                    # Check molecular descriptor
                    if molec[0].strip()==name:
                        lig_present = True
                        print('# INFO: Gathering position information from {}.'.format(posfl))
                        descript = descript if len(molec[6].strip())==0 else molec[6].strip()
                        all_ats  = True # Consider all atoms
                    else:
                        continue
                    
                molfile = posfl
                molmap  = {}
                # Look for named molecules
                for atom in atomic:
                    if not all_ats:
                        if atom[resnm].strip() != name:
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
                        else:
                            pos_info[-1][key] = float(sect.strip())            
                    
            else:
                if posfl[-3:] == 'pdb':
                    keywords = {'name'  : [12, 16],
                                'x'     : [30, 38],
                                'y'     : [38, 46],
                                'z'     : [46, 54],
                                'elem'  : [76, 78]}
                    resnm = [17, 21]
                
                elif posfl[-3:] == 'crd':    
                    keywords = {'name'  : [32,   40],
                                'x'     : [41,   60],
                                'y'     : [61,   80],
                                'z'     : [81,  100],
                                'elem'  : [101, 102]} # Empty element identifier to force self-identification
                    resnm =  [22, 30]
    
                elif posfl[-3:] == 'gro':
                    keywords = {'name'  : [10, 15],
                                'x'     : [20, 28],
                                'y'     : [28, 36],
                                'z'     : [36, 44],
                                'elem'  : [14, 15]} # Empty element identifier to force self-identification - assume 4 letter code only
                    resnm =  [5, 10]
                
                pos.split('END')[0]
                # Consider if a similar name has been used - gormatting depends on name length and file type
    
                pos           = pos.split('\n')[:-1]
                element_names = True
                
                for line in pos:
                    # Only look for named molecules
                    if line[resnm[0]:resnm[1]].strip() != name:
                        # Reject similar names
                        continue
                    if line.count('TER') == 1:
                        # Termination of residue chain, skip
                        continue
                    
                    if not lig_present:
                        print('# INFO: Gathering position information from {}. Note it is assumed ligand information is only present as a single ligand in one position file.'.format(posfl))
                        lig_present = True
                    # Append relevant information
                    pos_info.append({})
                    for i, key in enumerate(keywords.keys()):
                        sect = line[keywords[key][0]:keywords[key][1]]
                        if key == 'name':
                            pos_info[-1][key] = sect.strip()
                        elif key == 'elem':
                            pos_info[-1][key] = sect.strip()
                            if len(pos_info[-1]['elem']) == 0:
                                if element_names:
                                    print('# WARNING: Element names are missing - attempting to infer from atom names. Note that this may cause issues.')
                                    element_names = False
                                pos_info[-1]['elem'] = get_elements(pos_info[-1]['name'])
                        else:
                            pos_info[-1][key] = float(sect.strip())
            if lig_present:            
                break            
    
        position_data = pd.DataFrame.from_dict(pos_info, orient='columns')
        bonds         = []
        first_bond    = True ; first_unknown = True
        mol2bond_map  = {'1' : ['SING', 'N'], '2': ['DOUB', 'N'], '3': ['TRIP', 'N'], 'ar': ['AROM', 'Y']} 
    
        # ======================= #
        # Get bonding information #
        # ======================= #
        
        # Scan through RTP/RTF/ITP/MOL2 files
        # ------------------------------------
        for bndfl in bonding_files:
            if molfile != None:
                bondfile = open(molfile).read()
    
                # Open and read mol2 file
                # ------------------------            
                mbonds = bondfile.split('@<TRIPOS>BOND')[1].split('@')[0].split('\n')
                mbonds = [[b for b in line.split(' ') if len(b)!=0] for line in mbonds if len(line)!=0]
                
                for line in mbonds:
                    currbond = [molmap[line[1]], molmap[line[2]]]
                    try:
                        currbond.extend(mol2bond_map[line[3]])
                    except KeyError:
                        if first_unknown:
                            print('# WARNING: Allocating unknown bonds as single bonds - this should not affect output.')
                            first_unknown = False
                        currbond.extend(mol2bond_map['1'])
                    bonds.append(currbond)
                break
    
            else:
                if first_bond:
                    print('# WARNING: Allocating all bonds as single bonds - this should not affect output.')
                    first_bond = False
                    
                bondfile = open(bndfl).read()
                if bondfile.count(name) == 0 and bondfile.count(name.lower()) == 0:
                    # Note crd file from CHARMMGUI may convert to lowercase
                    continue
                
                if bndfl[-3:] == 'rtp':
                    # Open and read rtp file
                    # ----------------------
                    if bondfile.count('[ {} ]'.format(name)) == 0:
                        # Check for similar names
                        continue
                    
                    print('# INFO: Gathering bonding information from {}.'.format(bndfl))
                    RTP  = bondfile.split('[ {} ]'.format(name))[1] # Strip previous molecules - consider presence of similar names
                    RTP      = RTP.split(']')[:3]                           # Select atoms and bonds
                    descript = RTP[0].split('\n')[1].strip(';').strip()  # Name may be below the CHARMM code
    
                    # Get charges
                    atom_data = RTP[1].split('\n')[1:-1]
                    for line in atom_data:
                        if len(line)==0:
                            continue
                        info = [split for split in line.strip().split(' ') if len(split)!=0]
                        if len(info)==0 or info[0][0] == ';':
                            # Skip comments
                            continue
                        position_data.loc[position_data['name']==info[0].strip(), 'charge'] = get_charge(float(info[2].strip()), args.charge)
            
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
                    break
    
                elif bndfl[-3:] == 'rtf':
                    # Open and read rtp file
                    # ----------------------
                    if bondfile.count(' {} '.format(name)) == 0 and bondfile.count(' {} '.format(name.lower())) == 0:
                        # Check for similar names
                        continue
                
                    print('# INFO: Gathering bonding information from {}.'.format(bndfl))
                    RTF = bondfile.split('\n')
                    for line in RTF:
                        if line[:4] == 'ATOM':
                            # Get charges
                            info = [split for split in line.strip().split(' ') if len(split)!=0]
                            position_data.loc[position_data['name']==info[0], 'charge'] = get_charge(float(info[3]), args.charge)
            
                        elif line[:4] == 'BOND':
                            info = [split for split in line.strip().split(' ') if len(split)!=0]
                            bonds.append([info[1], info[2], 'SING', 'N'])
                    break
    
    
                elif bndfl[-3:] == 'itp':
                    # Open and read itp file
                    # ----------------------
                    if bondfile.count(' {} '.format(name)) == 0 and bondfile.count(' {} '.format(name.lower())) == 0 :
                        # Check for similar names
                        continue
    
                    
                    print('# INFO: Gathering bonding information from {}.'.format(bndfl))
    
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
                        position_data.loc[position_data['name']==atom[4], 'charge'] = get_charge(float(atom[6]), args.charge)
    
                    # Get bonds
                    # ---------
                    ITP = bondfile.split('[ bonds ]')[1].split('[')[0].split('\n')
                    ITP = [[b for b in line.split(' ') if len(b)!=0] for line in ITP if len(line)!=0]
    
                    for line in ITP:
                        if line[0][0] == ';':
                            # Skip comments
                            continue
                        bonds.append([molmap[line[0]], molmap[line[1]], 'SING', 'N'])
            
                    break
    
                elif bndfl[-4:] == 'mol2':
                    # Not molfile
                    continue
    
        # Generate bonding information if needed
        # --------------------------------------
                
        if len(bonds) == 0:
            # Specify bonds from pos file - either no RTP/RTF file or wrong.
            # Also do I want to get information about the bond order (single/double/triple)?
    
            # Note that AF README states that the bond order and aromacity don't matter
    
            print('# WARNING: Inferring bonding from proximity, ignoring Hs.')
            print('# WARNING: Giving all atoms a charge of 0.')
            
            if first_bond:
                print('# WARNING: Allocating all bonds as single bonds - this should not affect output.')
                first_bond = False
    
            position_data.loc[:, 'charge'] = 0
            noHdata = position_data.loc[position_data['elem'] != 'H'].to_dict(orient='records')
    
            bonds = get_bonds_dist(noHdata)
    
        if args.Hydrogen:
            posdata = position_data.to_dict(orient='records')
            Hs      = []
        else:
            posdata = position_data.loc[position_data['elem'] != 'H'].to_dict(orient='records')
            Hs      = position_data.loc[position_data['elem'] == 'H', 'name'].to_list()
            
        # ======================================== #  
        # Write mmCIF output in the desired format #
        # ======================================== #
    
        # Check renaming
        # --------------
        
        print('# INFO: Writing CCD data to {}_output.cif for residue {}'.format(nname, nname)+extra_info)
        
        # Write cif data
        # -------------
        cif_content = cif_information(nname, descript, posdata)
    
        for bond in bonds:
            if Hs.count(bond[0])==0 and Hs.count(bond[1])==0:
                cif_content.append(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}")
        cif_content.append('#')
        
        # Write to CIF file
        with open(nname+'_output.cif', 'w') as cif_file:
            cif_file.write("\\n".join(cif_content))
        cif_file.close()
            
        # Write to JSON file
        json_content = "\\n".join(cif_content)
        json.write(json_content)
        json_ligands.append(nname)
    
        
    # ------------------------------------------------------------- #
    #                                                               #
    # Work through list of molecules made of covalent modifications #
    #                                                               #
    # ------------------------------------------------------------- #
    
    if len(args.covalent)!=0:
        print('# WARNING: No other molecules can be present in files for covalently bonded molecules. If using pdb files must provide itp after it.')
    
    first_component = True
    
    # ===============
    
    # For each molecule I need to:
    # 0.   Determine the ligand name
    # 1.   Determine the different components which need to be added to userCCD
    # 1.a) Make sure that any duplicates within the molecule are ignored (take the first one?)
    # 1.b) Make sure that any duplicates outwith the molecule are ignored (take the one outside of the molecule) - also want to add this to the file for the molecule but probably take the one in the molecule
    # 2.   Determine the covalent modifications between components
    # 3.   Add the necessary components to userCCD
    # 4.   Add the necessary modifications to userCCD
    
    full_ligand_information = []
    full_modification_info  = []
    
    for Nfile, molfile in enumerate(args.covalent):
        chain = allchains.pop(0)
        if molfile[-3:]=='itp':
            continue
        print('# INFO: Creating necessary user-defined CCD codes and covalent modifications for the ligand in {}. Note that bonding information must be provided.'.format(molfile))
        
        full_bonds = []
        
        if molfile[-4:] == 'mol2':
            # Open and read mol2 file - must have substructure information
            # -------------------------------------------------------------
            keywords = {'ID'     : 0,
                        'name'   : 1,
                        'x'      : 2,
                        'y'      : 3,
                        'z'      : 4,
                        'elem'   : 5,
                        'resi'   : 6,
                        'resnm'  : 7}
    
            data = open(molfile).read()
            
            # Get full name of ligand from @<TRIPOS>MOLECULE
            # -----------------------------------------------
            ligname = data.split('@<TRIPOS>MOLECULE')[1].split('@')[0].split('\n')[1].strip()
    
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
                    if first_unknown:
                        print('# WARNING: Allocating unknown bonds as single bonds - this should not affect output.')
                        first_unknown = False
                        currbond.extend(mol2bond_map['1'])
                        full_bonds.append(currbond)
    
        elif molfile[-3:] == 'pdb':
            keywords = {'ID'    : [6, 11],
                        'name'  : [12, 16],
                        'resnm' : [17, 21],
                        'resi'  : [22, 26],
                        'x'     : [30, 38],
                        'y'     : [38, 46],
                        'z'     : [46, 54],
                        'elem'  : [76, 78]}
    
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
    
            atomic = pd.DataFrame.from_dict(pos_info, orient='columns') ; atomic = atomic.set_index('ID')
            components = atomic['resnm'].unique()
            
            # ======================= #
            # Get bonding information #
            # ======================= #
    
            assert args.covalent[Nfile+1][-3:] == 'itp', 'ERROR: missing bonding information for covalently modified ligand'     
            bondfile = open(args.covalent[Nfile+1]).read()
    
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
    
                        
        flbnds = pd.DataFrame(full_bonds, columns = ['atom1', 'resi1', 'atom2', 'resi2', 'type', 'ar'])
        chains = []
            
        # =========================== #
        # Intra-component information # 
        # =========================== #
    
        molecule_file = open(ligname+'_output.json', 'w')
        molecule_file.write('{\n\t"userCCD" : "')    
    
        resimapping = {}
        for counter, resi in enumerate(list(atomic['resi'].unique())):
            chains.append([counter+1, chain, atomic.loc[atomic['resi']==resi, 'resnm'][0]])
            # AF3 indexes from 0, add residues to the same chain
            resimapping[resi] = str(counter+1)
    
        for j, component in enumerate(components):
                
            # Determine if multiple components are the same
            resis = list(atomic.loc[atomic['resnm']==component, 'resi'].unique())
            if len(resis) > 1 and first_component:
                print('# INFO: Where multiple of the same component are present the first component present in the file is used to construct userCCD')
                first_component = False
            
            component_data        = atomic.loc[atomic['resi']==resis[0]]
            intra_component_bonds = flbnds.loc[flbnds['resi1']==resis[0]].loc[flbnds['resi2']==resis[0]]
                    
            # Write userCCD data for component only 
            # --------------------------------------
                    
            if args.Hydrogen:
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
            molecule_file.write("\\n".join(cif_content))
            # Determine if already in JSON file
            if args.names.count(component.upper())==0:
                # Write to JSON file - not already present
                json_content = "\\n".join(cif_content)
                json.write(json_content)
            first_CCD = False
    
        # =================== #
        # Write modifications # 
        # =================== #
    
        # Close userCCD section
        
        molecule_file.write('",\n  "sequences": [ \n')
    
        chains = pd.DataFrame(chains, columns=['resi', 'chain', 'resnm']) ; chains = chains.set_index('resi')
        
        # Add correct number of ligands - want to keep a record of this for JSON
        # ----------------------------------------------------------------------
        lig = '   {"ligand": {"id": ["'+chain+'"], "ccdCodes": ["'+'", "'.join(chains['resnm'])+'"]}}'
    
        full_ligand_information.append(lig)
        molecule_file.write(lig+', \n')
            
        # Add modifications - currently won't work for proteins but could extend in future
        # --------------------------------------------------------------------------------
    
        modification_info = []
        molecule_file.write('\n    ], \n   "bondedAtomPairs": [\n')
        for row, bond in flbnds.loc[flbnds['resi1']!=flbnds['resi2']].iterrows():
            # Assume all ligands at residue 1
            mod  = '    [["'+chain+'", '+resimapping[bond['resi1']]+', "'+bond['atom1']+'"],'
            mod +=      '["'+chain+'", '+resimapping[bond['resi2']]+', "'+bond['atom2']+'"]]'
            modification_info.append(mod)
    
        full_modification_info.extend(modification_info)
    
        molecule_file.write(',\n'.join(modification_info))
        molecule_file.write('],')
            
        molecule_file.close()                
    
    json.write('",\n')
    
    # Write preamble
    # --------------
    
    json.write('\t"name": "'+args.title+'",\n')
    json.write('\t"dialect": "'+args.dialect+'",\n')
    json.write('\t"version": '+args.afvers+',\n')
    seeds = [s.strip() if s[-1] != ',' else s[:-1].strip() for s in args.seeds]
    json.write('\t"modelSeeds": ['+', '.join(seeds)+'],\n')
    json.write('\t"sequences": [\n\t\t')
    
    
    print('# INFO: Assuming 1 of every ligand.')
    
    if len(args.covalent)!=0:
        # ==================== #
        # Write full JSON file # 
        # ==================== # 
    
        # Close userCCD section
        json.write(', \n'.join(full_ligand_information))
        json.write('\n    ], \n   "bondedAtomPairs": [')
        json.write(',\n'.join(full_modification_info))
        json.write('\n      ],')
    
    
    # Write protein sequences
    # -----------------------
    
    if len(args.protein)==0:
        # Write ligand information only
        pass
    else:
        if len(args.protein)==1:
            # Either a single sequence or a single fasta file   
            if '.' in args.protein[0]:
                # Fasta file  
                protein = read_fasta(args.protein[0])
            else:
                protein = args.protein
        else:
            # Split protein input to determine unique chains                                                         
            protein  = []
            in_range = True ; i = 0
            while in_range:
                try:
                    num = int(args.protein[i+1])
                    if '.' in args.protein[i]:
                        protein_i = read_fasta(args.protein[i])
                    else:
                        protein_i = [args.protein[i]]
                    protein.extend(list(protein_i)*num)
                    i += 2
                except (ValueError, IndexError) as e:
                    if '.' in args.protein[i]:
                        protein_i = read_fasta(args.protein[i])
                    else:
                        protein_i = args.protein[i]
                    protein.extend(list(protein_i))
                    i += 1
                in_range = False if i >= len(args.protein) else True        

        seqs   = Counter(protein)
        for sequence in seqs.keys():
            json.write('{"protein": {"id" : ["')
            for i in range(seqs[sequence]):
                json.write(allchains.pop(0)+'"')
                if i != seqs[sequence]-1:
                    json.write(', "')
            json.write('],\n\t\t\t"sequence": "'+sequence+'"}},\n\t\t')
        
    for i, lig in enumerate(json_ligands):        
        json.write('{"ligand": {"id": ["'+allchains.pop(0)+'"], "ccdCodes": ["'+lig+'"]}}')
        if i == len(json_ligands)-1:
            json.write(']\n}')
        else:
            json.write(',\n\t\t')
            
    json.close()


if __name__ == "__main__":
    main()
