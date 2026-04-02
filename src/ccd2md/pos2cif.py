###########################################
#                                         # 
#    Convert between CRD/PDB and mmCIF    # 
#                                         #
###########################################

# Import relevant functions
# -------------------------

import argparse
from ccd2md import FuncPos
import pandas as pd
import numpy as np
import sys, os, subprocess
import warnings
import string
from collections import Counter
from glob import glob

def main():
    
    # Get command line arguments
    # ---------------------------
    
    parser = argparse.ArgumentParser(description='Generate user-defined CCD code(s) for use in AF3 from position and optionally additional bonding file(s). Generate AF3 input files from specified parameters.', add_help=False)

    info = parser.add_argument_group('pos2cif information')
    info.add_argument('-v', '--version', action='version', version='Version '+FuncPos.__version__)
    info.add_argument('-h', '--help',    action='help',    help='Show this message and exit.')
    
    info.add_argument('-CF',  '--configuration',   help='Specifiy configuration file containing command line arguments. For an example input please see keb721/CCD2MD on GitHub.', default=None)
    info.add_argument('-ncl', '--no_command_line', help='Prevents addition to multiply defined arguments (seeds/protein/ligand/userCCDPath/name/files/rename/covalent) from the command line. The presence of these will lead to a warning and the command line input(s) will be ignored. The presence of any other flags will still cause the programme to exit with an error.', nargs='?', const=True, default=None)


    cif_generation = parser.add_argument_group('CCD code creation information')

    cif_generation.add_argument('-n', '--names',    help='CHARMM code for molecule(s) to convert, note this must not overlap with an existing CCD code.', nargs='+', default = [])
    cif_generation.add_argument('-r', '--rename',   help='Pair(s) of ligand names in format OLD NEW, where OLD is present in the position files and NEW is the desired new names. Note that this is currently only applicable to userCCD generation, not userCCD usage.', nargs='+', default=[])
    cif_generation.add_argument('-f', '--files',    help='.gro/.pdb/.crd/.mol2 file(s) containing molecule(s) (note, can contain other moleucles). Optionally may also specify .rtp/.rtf/.itp file(s) containing molecule(s) (note, rtp files can contain other moleucles). If .rtp/.rtf files not provided, bonded information will be inferred from proximty.', nargs='+', default = [])
    cif_generation.add_argument('-c', '--covalent', help='Position and bonding information for ligands which are to be constructed using covalent modifications. Either mol2 file (including bonding information), or a pair of pdb and itp files.', nargs='+', default = [])

    
    opts = parser.add_argument_group('optional arguments/cutoffs for CCD creation')

    opts.add_argument('-e', '--charge',   help='Minimum charge required for the charge to be non-0 (output integer charges only). Default 0.75', type=float, default=None)
    opts.add_argument('-b', '--bond',     help='Maximum distance (in A) between atoms to be considered bonded. Note: only utilised if no bonding information is specified in files. Default 1.8', type=float, default=None)
    opts.add_argument('-H', '--Hydrogen', help='Retain hydrogens in mmCIF data. Default=False.', nargs='?', const=True, default=None)

    
    af3_system = parser.add_argument_group('system information for AlphaFold3')
    
    af3_system.add_argument('-p', '--protein',          help='FASTA protein sequence(s) or file(s) to add to system. For multiple of the same sequence (e.g. AACCS) can be "AACCS AACCS" or "AACCS 2". Fasta files may contain multiple sequences but each sequence must start with an information line beginning ">". "-p Test.fasta 3" will insert three copies of all sequences within "Test.fasta". A mixture of sequences and files may be used.', default = [], nargs='+')
    af3_system.add_argument('-ptm', '--post_trans_mod',  help='Add post translational modifications. Each PTM should be specified in the style "A 12 LYSM" (i.e., chain resID CCD). Protein chains are labelled in alphabetical order starting from A, resIDs start from 1, and CCD codes must be specified.', default = [], nargs='+')
    af3_system.add_argument('-l', '--ligand',            help='Ligands (and numbers) to add to the system. These can be userCCD codes or CCD codes. Note that if unset, one copy of every converted ligand will be added.', nargs='+', default = [])
    af3_system.add_argument('-u', '--userCCDPath',       help='Locations of userCCD files to add. Note that CCD2MD.cif is added automatically unless disabled via the "-nC/--no_CCD2MD" flag.', nargs='+', default = [])
    af3_system.add_argument('-nC', '--no_CCD2MD',        help='Prevent use of CCD2MD.cif.', nargs='?', const=True, default=None)    

    
    af3_setup = parser.add_argument_group('setup information for AlphaFold3')

    af3_setup.add_argument('-j', '--json',     help='Name of JSON file to write output to. Default = "output.json"', default=None)
    af3_setup.add_argument('-t', '--title',    help='AF3 system title. Default = "pos2cif_system". Will also create a files containing userCCD codes, "{title}_CCD.cif", in current directory', default=None)
    af3_setup.add_argument('-A', '--afvers',   help='AF3 version. Default = 4', default=None)
    af3_setup.add_argument('-s', '--seeds',    help='Model seeds - need not be comma separated. Default 1', default = [], nargs='+')
    af3_setup.add_argument('-d', '--dialect',  help='Dialect. Default "alphafold3"', default=None)
    
 

    args = parser.parse_args()
    
    if args.configuration != None:
        args = FuncPos.read_configuration_file(args.configuration, args)


    # Restore defaults
    # -----------------
    
    default_dict = {'charge'   : 0.75,
                    'bond'     : 1.8,
                    'Hydrogen' : False,
                    'json'     : 'output.json',
                    'title'    : 'pos2cif_system',
                    'afvers'   : '4', # Required for CCDPath
                    'dialect'  : 'alphafold3'}


    for entry in default_dict.keys():
        if args.__dict__[entry] == None:
            args.__dict__[entry] = default_dict[entry]
            

    if args.seeds == []:
        args.seeds = ['1']

    curr_ligands = []
        
    # Input parameters parsed

    # Determine information for userCCD generation
    # ---------------------------------------------

    if len(args.covalent) != 0 or len(args.names) != 0:
        # Requesting conversion - check input files

        # Sense checks on inputs
        # -----------------------

        position_files, bonding_files = FuncPos.get_files(args.files)
        if len(args.names) != 0 and len(position_files) == 0 :
            print('ERROR: Cannot generate userCCD codes for {} - no position files provided.'.format(args.names))
            exit()
            
        assert len(args.rename)%2==0, 'ERROR: rename must take in pairs of ligand names.'

        print('# WARNING: It is assumed ligand information is only present as a single ligand in one position file.')

        # Define name changes
        # --------------------
        
        new_names = {} 
        for i in range(0, len(args.rename), 2):
            new_names[args.rename[i]] = args.rename[i+1].upper()
            
        # ======================================================== #
        # Generate userCCD codes from position/bonding information #
        # ======================================================== #
                
        args.names = [name.upper() for name in args.names]
        for i, name in enumerate(args.names):
            
            # Check renaming
            # --------------

            nname = new_names.get(name, name)
            if nname == name:
                print('# INFO: Creating user-defined CCD code for {}'.format(name))
                extra_info = ''
            else:
                print('# INFO: Creating user-defined CCD code for {}, which has been renamed from {}'.format(nname, name))
                extra_info = ' (renamed from {})'.format(name)

            descript, pos, bonds, Hs = FuncPos.CCD_from_name(name, position_files, bonding_files, args.bond, args.charge)

            # ======================================== #  
            # Write mmCIF output in the desired format #
            # ======================================== #
    
            # Check renaming
            # --------------
        
            print('# INFO: Writing CCD data to {}_output.cif for residue {}'.format(nname, nname)+extra_info)
        
            # Write cif data
            # -------------
            cif_content = FuncPos.cif_information(nname, descript, pos)
    
            for bond in bonds:
                if Hs.count(bond[0])==0 and Hs.count(bond[1])==0:
                    cif_content.append(f"{bond[0]} {bond[1]} {bond[2]} {bond[3]}")
            cif_content.append('#')
        
            # Write to CIF file
            with open(nname+'_output.cif', 'w') as cif_file:
                cif_file.write("\n".join(cif_content))
            cif_file.close()

            curr_ligands.append(nname)
            
        # Single molecule conversion complete

    # ============================ #
    #  Covalently modified ligands #
    # ============================ #
    
    if len(args.covalent)!=0:
        print('# WARNING: No other molecules can be present in files for covalently bonded molecules. If using pdb files must provide itp after it.')
        
        for Nfile, molfile in enumerate(args.covalent):
            chain = 'A'
            if molfile[-3:]=='itp':
                continue
            print('# INFO: Creating necessary user-defined CCD codes and covalent modifications for the ligand in {}. Note that bonding information must be provided.'.format(molfile))
        
            if molfile[-4:] == 'mol2':
                atomic, full_bonds, components = FuncPos.read_covalent_mol2(molfile)
                
            elif molfile[-3:] == 'pdb':
                # Note currently only applicable to pdb/itp format - may change in future

                assert args.covalent[Nfile+1][-3:] == 'itp', 'ERROR: missing bonding information for covalently modified ligand'     

                atomic, components = FuncPos.get_covalent_position(molfile)
                atomic, full_bonds = FuncPos.get_covalent_bonding_itp(args.covalent[Nfile+1], atomic)
    
                        
            flbnds = pd.DataFrame(full_bonds, columns = ['atom1', 'resi1', 'atom2', 'resi2', 'type', 'ar'])
            
            FuncPos.covalent_cif(ligname, atomic, full_bonds, components, args.Hydrogen)
            curr_ligands.append(ligname)
            
    if len(args.files) != 0 and len(args.names)==0:
            print('# WARNING: files ({}) for creation of userCCD codes have been added but will not be parsed. To generate userCCD codes please add the flag -n NAME for NAME in files.'.format(args.files))
    if len(args.names) ==0 and len(args.rename) != 0:
        print('# WARNING: renaming is only applicable to userCCD generation. The specified renaming(s) ({}) will not be performed when generating AF3 input.'.format(args.rename))
            
    # ================ #
    #    WRITE JSON    #
    # ================ #

    print('# INFO: Writing output json, {}'.format(args.json))
    json = open(args.json, 'w')
    json.write('{\n')
            
    # Write preamble
    # --------------
    
    json.write('\t"name": "'+args.title+'",\n')
    json.write('\t"dialect": "'+args.dialect+'",\n')
    json.write('\t"version": '+args.afvers+',\n')
    seeds = [s.strip() if s[-1] != ',' else s[:-1].strip() for s in args.seeds]
    json.write('\t"modelSeeds": ['+', '.join(seeds)+'],\n')
    json.write('\t"userCCDPath": "{}_CCD.cif",\n'.format(args.title))
    json.write('\t"sequences": [\n\t\t')

    # Write protein sequences and modifications 
    # -----------------------------------------

    add_PTMs = []

    base_chains = np.array(list(string.ascii_uppercase))

    allchains = list(string.ascii_uppercase)
    for entry in base_chains:
        allchains.extend(list(entry+base_chains))

    allchains = [str(item) for item in allchains]
    
    
    if len(args.protein)==0:
        # Write ligand information only
        if len(args.post_trans_mod) != 0:
            print('# WARNING: PTMs were sepcified but there are no proteins to apply this to.')
        pass
    else:
        if len(args.protein)==1:
            # Either a single sequence or a single fasta file   
            if '.' in args.protein[0]:
                # Fasta file  
                protein = FuncPos.read_fasta(args.protein[0])
            else:
                protein = args.protein
        else:
            # Split protein input to determine unique chains
            protein = FuncPos.get_components(args.protein)
            
        seqs      = Counter(protein)

        if len(args.post_trans_mod)!=0:
        
            assert len(args.post_trans_mod)%3==0, 'ERROR: PTMs must be triplets of chain ID, residue number and CCD code.'

            PTMs = np.array(args.post_trans_mod).reshape(int(len(args.post_trans_mod)/3), 3)
            PTMs = pd.DataFrame(PTMs, columns = ['chain', 'resID', 'CCD'])  # May not be in order

            if len(PTMs)!=0:
                print('# WARNING: No checks are made between specified PTM and amino acid at specified sequence position.')
                print('# WARNING: Covalently bonded PTMs are currently not supported.')
                
            for sequence in seqs.keys():
                        
                # Determine all relevant chains

                protein_chains = allchains[:seqs[sequence]]
                allchains      = allchains[seqs[sequence]:]

                while len(protein_chains) > 0:
                    json.write('{"protein": {"id" : ["')
        
                    firstchain = protein_chains.pop(0)
                    json.write(firstchain+'"')
                    
                    firstPTM = PTMs.loc[PTMs['chain']==firstchain]
                    
                    reduced_chain = []
                
                    for i in range(len(protein_chains)):
                    
                        currchain = protein_chains[i]

                        diffPTM = pd.merge(firstPTM,PTMs.loc[PTMs['chain']==currchain],how='outer',on=['resID','CCD'],indicator=True)

                        # All both indicates the same modifications
                        
                        if sum(diffPTM['_merge'] == 'both') == len(diffPTM):
                            json.write(', "'+currchain+'"')
                        else:
                            reduced_chain.append(currchain)


                    # Gone through same modifications as first chain
                    json.write('],\n\t\t\t"sequence": "'+sequence+'"')
                    if len(firstPTM) != 0:
                        json.write(',\n\t\t"modifications": [')
                        numPTMs = len(firstPTM)
                        for i in range(numPTMs):
                            json.write('\n\t\t\t{"ptmType": "'+firstPTM.iloc[i]['CCD']+'", "ptmPosition": '+firstPTM.iloc[i]['resID'])
                            add_PTMs.append(firstPTM.iloc[i]['CCD'])
                            if i == numPTMs-1:
                                json.write(']')
                            else:
                                json.write(', ')
                
                    json.write('}},\n\t\t')                
                
                    protein_chains = reduced_chain

        else:
            for sequence in seqs.keys():
                        
                # Determine all relevant chains

                protein_chains = allchains[:seqs[sequence]]
                allchains      = allchains[seqs[sequence]:]

                json.write('{"protein": {"id" : ["'+'", "'.join(protein_chains)+'"],')
                json.write('\n\t\t\t"sequence": "'+sequence+'"')
                json.write('}},\n\t\t')                
                
    # Completed writing proteins

    # Check no PTMs applied to non-protein
    if len(args.post_trans_mod)!=0:
        assert len(PTMs.loc[PTMs['chain'].isin(allchains)])==0, 'ERROR: PTM applied to a chain which does not correspond to a protein'

    # ============= #
    # WRITE LIGANDS #
    # ============= #

    if args.ligand == []:
        args.ligand = curr_ligands
    
    curr_ligands = pd.DataFrame(curr_ligands[::-1], columns=['CCD'])

    # Open file for CCD library
    # -------------------------

    exists = os.path.exists('{}_CCD.cif'.format(args.title)) # Check if file already exists
    if exists:
        backup = glob('#{}_CCD.cif.*#'.format(args.title))
        subprocess.run(['scp', '{}_CCD.cif'.format(args.title), '#{}_CCD.cif.{}#'.format(args.title, len(backup)+1)])

    CIF = open('{}_CCD.cif'.format(args.title), 'w')

    bonded_pairs = []
                        
    # First need to consider if any of the ligands to add are covalently bonded as these are the most difficult!

    ligs     = FuncPos.get_components(args.ligand)
    ligands  = Counter(ligs)
    for PTM in set(add_PTMs):
        ligands[PTM] = 0
    add_json = {lig: True for lig in list(ligands.keys())}
        
    inCIF   = {}

    base = '# INFO: CCD codes will be pulled from '
    if args.no_CCD2MD:
        base = base + 'the ligands just converted.' if len(args.userCCDPath)==0 else base + '(in order): the ligands just converted; and the ligands in {}.'.format(', '.join(args.userCCDPath))
    else:
        base += '(in order): the ligands just converted;'
        if len(args.userCCDPath)!=0:
            base += ' the ligands in {}; '.format(', '.join(args.userCCDPath))
        base += ' and the ligands in CCD2MD.cif.' 
        
    base += ' Ligands present in multiple locations will be added from the first location. Any ligands not present will be assumed to be conventional CCD codes.'

    print(base)                                                                                                                                                                                                                 
    
    # Read inputs
    # Consider in the order:
    # 1. Ligands defined in this session (start with covalent)
    # 2. Ligands in path specified locations
    # 3. Ligands in CCD2MD.cif

    # Check ligands defined in this session
    # --------------------------------------

    this_session = list(curr_ligands.loc[curr_ligands['CCD'].isin(add_json)].to_numpy().flatten())
        
    for lig in this_session:
        currlig = open('{}_output.cif'.format(lig), 'r').read()

        numligs = ligands.pop(lig)
        
        lig_chains = allchains[:numligs]
        allchains  = allchains[numligs:]

        if currlig[0] == '{':
            
            json, CIF, inCIF, bonding_info = FuncPos.extract_covalent(currlig, lig, json, inCIF, CIF,
                                                                          args.title, lig_chains, True)
            bonded_pairs.append(bonding_info)

        else:
            # Will only be a single ligand
            if inCIF.get(lig) == None:
                inCIF[lig] = '{}_output.cif (this session)'.format(lig)
                CIF.write(currlig)
            else:
                print('# WARNING: {} is defined in {}_output.cif but {}_CCD.cif contains a version from ,'.format(lig, lig))

            
            if len(lig_chains)!=0:    
                json.write('{"ligand": {"id": ["'+', '.join(lig_chains)+'"], "ccdCodes": ["{}"]}}'.format(lig))
                
        if len(ligands) != 0:
            json.write(',\n\t\t')                
            
        add_json[lig] = False

            
    # Ok, now want to go through the rest of the files and determine what needs to be added
    # --------------------------------------------------------------------------------------

    for ligfile in args.userCCDPath:
        currlig = open(ligfile, 'r').read()

        if currlig[0] == '{':
            # Determine if contents are covalent or not - note assumptions here on name etc.

            test_ligname = ligfile.split('_output.cif')[0]
            print('# WARNING: Assuming that {} contains a covalently bonded molecule of the name {}'.format(ligfile, test_ligname))
                
            # Test if this ligand is to be considered
            if add_json.get(test_ligname)!=None:
                if add_json[test_ligname]:

                    numligs = ligands.pop(test_ligname)
                    
                    lig_chains = allchains[:numligs]
                    allchains  = allchains[numligs:]
                    
                    json, CIF, inCIF, bonding_info = FuncPos.extract_covalent(currlig, test_ligname, json, inCIF,
                                                                          CIF, args.title, lig_chains)
                    bonded_pairs.append(bonding_info)

                    add_json[test_ligname] = False
                
                    if len(ligands) != 0:
                        json.write(',\n\t\t')
                else:
                    print('# WARNING: {} is defined in {} but {}_CCD.cif contains a version from {}'.format(test_ligname, ligfile, title, inCIF[ligname]))
          
        else:
                
            json, add_json, ligands, allchains, CIF, inCIF = FuncPos.extract_ligand(currlig, add_json, allchains,
                                                                                    inCIF, CIF, json, args.title,
                                                                                    ligands, ligfile)

                
    if not args.no_CCD2MD:
        # Check CCD2MD.cif - no covalent ligands

        # currlig = open(os.path.dirname(ccd2md.__file__)+'/CCD2MD.cif', 'r').read()  # Location of CCD2MD
        currlig = open('CCD2MD.cif', 'r').read()  # Location of CCD2MD

        json, add_json, ligands, allchains, CIF, inCIF = FuncPos.extract_ligand(currlig, add_json, allchains, inCIF,
                                                                                CIF, json, args.title, ligands,
                                                                                'CCD2MD.cif')
        

    # Add the other ligands, which may represent conventional CCD codes - check if they are not 2/3 letter codes.
    if len(ligands) != 0:
        missed_CCD = [code for code in ligands.keys() if len(code) > 3]
        if len(missed_CCD):
            print('# WARNING: The CCD codes {} have not been defined - these are presumed to be userCCD codes.')
        else:
            print('# INFO: 2/3 letter CCD codes undefined are presumed to be ions/conventional CCD codes.')
        conventional = list(ligands.keys())
        for lig in conventional:

            numligs = ligands[lig]
                        
            lig_chains = allchains[:numligs]
            allchains  = allchains[numligs:]    

            if numligs!=0:
                json.write('{"ligand": {"id": ["'+', '.join(lig_chains)+'"], "ccdCodes": ["{}"]}}'.format(lig))
            if lig == conventional[-1]:
                json.write(',\n\t\t')                

    if len(bonded_pairs)!=0:
        json.write('],\n\t"bondedAtomPairs": [')
        json.write(',\n'.join(bonded_pairs))
    
    json.write(']\n}')            
    json.close()


if __name__ == "__main__":
    main()
