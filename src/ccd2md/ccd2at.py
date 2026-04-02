#######################################
#                                     # 
#    Convert between CCD and CHARMM   # 
#                                     #
#######################################

# Import relevant functions
# -------------------------

from ccd2md import FuncConv
# import FuncConv
import argparse
import numpy as np
import pandas as pd
from copy import deepcopy
import sys, subprocess, os

def main():

    # Get command line arguments
    # ---------------------------
    parser = argparse.ArgumentParser(description='Convert from output from co-folding (specified CCD codes and SMILES strings, and userCCD codes) to CHARMM. Please see keb721/CCD2MD on github for a list of the allowed SMILES strings and CCD codes.', add_help=False)
    
    parser.add_argument('inputfile',  help='Input file name - .cif, .pdb or .gro.', nargs='?')
    parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.', nargs='?')
    
    parser.add_argument('-L', '--ligchain',   help='Output ligands in their own chains - default is off. Only applicable if NOT embedding the system in a membrane.', action='store_true')
    parser.add_argument('-S',  '--SMILES',    help='Used SMILES strings, list the order of the name of the ligand used. Note that when multiple of the same ligand are used this can be written either e.g. "POPE POPE" or "POPE 2".', nargs='+', default=None)
    parser.add_argument('-T',  '--title',     help='Optionally specify the title of the generated pdb file.', default = None)
    parser.add_argument('-ndx', '--make_ndx', help='Make a gromacs index file either at the end of the CCD2MD run, or before attempting to create an equilibration tpr file. Please note that no additional arguments can be added after this command and index file creation is currently required to be interactive. Please also note that index file creation may be necessary for correct equilibration of the system.', action='store_true')

    parser.add_argument('-C',  '--conc',  help='Concentration of ions in the system - default = 0.15. Default ions are Na and Cl. For membrane embedded proteins, charge balance is maintained; specify different ions via -mem to pass to Insane4MemPrO. For non-membrane proteins charge balance is not maintained; specify different ions via -I/--ions.', default=None)
    
    parser.add_argument('-g', '--gmx', help='Path to GROMACS exectuable - this will also be passed to CG2AT-lite. Default = "gmx".')
    parser.add_argument('-gh', '--pdb2gmx', help='Override ALL defaults of pdb2gmx and optionally pass extra arguments - note that this may require interactivity, and may be necessary for a starting MET. Default is topology in topol.top, OUTPUTNAME_H.pdb, TIP3P water, charmm36-ccd2md forcefield, and charged termini (excepting starting CYST or GLYM which are set to None). Only applicable if NOT embedding the system in a membrane; --cg2at may provide some functionality in this case.', action='store_true')

    
    sys_opts = parser.add_argument_group('membrane-embedded system options')
    sys_opts.add_argument('-mem', '--membrane',   help='Embed the converted system into a membrane. Note that this will lead to minor rearragements of any ligands. After this flag, it is possible to add the arguments to be passed to Insane4MemPrO - most notably specifying the composition of the upper and lower leaflet in the form "-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1" where the codes represent lipids and the numbers represent the ratio between them. The default is two leaflets of pure POPC.', action='store_true')
    sys_opts.add_argument('-mp',  '--mempro',     help='Additional arguments for embedding the protein in the membrane using MemPrO. Add any additional arguments after this flag - default 5 grid points and 15 minimisation operations.', action='store_true')
    sys_opts.add_argument('-mdef', '--memprod',   help='Additional arguments for calculating the deformation of the membrane around the protein using MemPrOD. If not included, no deformations will be calculated. Otherwise any additional parameters for MemPrOD may be passed after this flag - otherwise MemPrOD defaults are used.', action='store_true')
    sys_opts.add_argument('-ncpu', '--num_cpus',  help='Number of CPUs to use for membrane embedding. Default = 1.', default=None)
    sys_opts.add_argument('-at', '--cg2at',       help='Additional arguments to pass to CG2AT-lite. Add additional arguments after this flag - note that this may require interactivity as it removes assumption abpout fragment location.', action='store_true')
    # sys_opts.add_argument('-sol', '--water',    help='The water model to use - currently only TIP3P') 

    
    glob_opts = parser.add_argument_group('non-membrane protein system')
    glob_opts.add_argument('-B',  '--box',      help='Create a box for a non-membrane protein. The presence of this flag alone creates a box, solvates, and (if a non-zero concentration is specified) adds ions. After this flag, it is possible to add the arguments to be passed to gmx editconf. The default is "-c -d 2".', action='store_true')
    glob_opts.add_argument('-SV', '--solvate',  help='Solvate a box for a non-membrane protein. The presence of this flag alone creates a box, solvates, and (if a non-zero concentration is specified) adds ions. After this flag, it is possible to add the arguments to be passed to gmx solvate. The default is "-cs".', action='store_true')
    glob_opts.add_argument('-I',  '--ions',     help='Add ions for a non-membrane protein. The presence of this flag alone creates a box, solvates, and (if a non-zero concentration is specified) adds ions. After this flag, it is possible to add the arguments to be passed to gmx genion. The default is "-neutral -c {conc}" for conc specified via -C/--conc.', action='store_true')
    

    cgmd_opts = parser.add_argument_group('coarse-grained membrane-embedded system MD options')
    cgmd_opts.add_argument('-nl', '--newlipidome',               help='Use updated Martini 3 mappings for CG lipids as in DOI: 10.1021/acscentsci.5c00755. Default off.', action='store_true')
    cgmd_opts.add_argument('-CGEM', '--CG_energy_minimise',      help='Generate an energy minimisation for the CG system before converting back to atomistic. Additional arguments for grompp may be added after this flag - if none are added the default energy minimisation script will be used.', action='store_true')
    cgmd_opts.add_argument('-rCGEM', '--run_CG_energy_minimise', help='Add additional arguments for running CG energy minimistaion (via gmx mdrun) before conversion back to atomstic.', action='store_true')
    cgmd_opts.add_argument('-CGeq', '--CG_equil',                help='Generate an equilibration for the CG system before converting back to atomistic. Additional arguments for grompp may be added after this flag - if none are added the default CG equilibration script will be used.', action='store_true')
    cgmd_opts.add_argument('-rCGeq', '--run_CG_equil',           help='Add additional arguments for running CG equilibration (via gmx mdrun) before conversion back to atomstic.', action='store_true')

    
    CG_opts = parser.add_argument_group('martinize2 options - note of limited relevance outside CG MD')    
    # Add user-extensible martinize2 options
    protein_network = CG_opts.add_mutually_exclusive_group()
    protein_network.add_argument('-E', '--elastic', help='Use the elastic network for Martini when converting an atomisitic protein to coarse-grained. When not present, or only the flag is present only -elastic is passed to martinize2. Include any desired elastic network commands to be passed to martinize2 after this argument (this may but does not need to include -elastic).', action='store_true')
    CG_opts.add_argument('-M', '--martinize',       help='Add additional arguments to martinize2.', action='store_true')

    steer_opts = parser.add_argument_group('steering options - any flag will initiate steering')
    steer_opts.add_argument('-Stnm', '--steer_name', help='CHARMM names of ligands to steer. If not set, all identifiable ligands will be steered. If not set, the reference ligand will be attempted to be taken from CCD2MD.cif. The presence of this flag alone will run steering.', nargs='+')
    steer_opts.add_argument('-Stref', '--steer_ref', help='Reference ligand file(s) for steering. Must be a .pdb/.gro/.cif file. These files must contain only a single ligand. If not set, the reference ligand will be attempted to be taken from CCD2MD.cif. The presence of this flag alone will run steering.', nargs='+')
    steer_opts.add_argument('-St',    '--steer',     help='Run steering.', action='store_true')

    aamd_opts = parser.add_argument_group('MD options')
    aamd_opts.add_argument('-AAEM', '--AA_energy_minimise',      help='Generate an energy minimisation tpr file for the final AA system. Additional arguments for grompp may be added after this flag - if none are added the default energy minimistation script will be used.', action='store_true')
    aamd_opts.add_argument('-rAAEM', '--run_AA_energy_minimise', help='Run energy minimisation for the final AA system (will also happen if equilibration or production steps are added). Additional arguments for mdrun may be added after this flag.', action='store_true')
    aamd_opts.add_argument('-AAeq', '--AA_equil',                help='Generate an equilibration tpr file for the final AA system. Additional arguments for grompp may be added after this flag - if none are addded the default AA equilibration script will be used.', action='store_true')
    aamd_opts.add_argument('-rAAeq', '--run_AA_equil',           help='Run equilibration for the final AA system (will be run if production tpr to be generated). Additional arguments for mdrun may be added after this flag.', action='store_true')
    aamd_opts.add_argument('-AAprd', '--AA_prod',                help='Generate a production tpr file for the final AA system. Additional arguments for grompp may be added after this flag - if none are added the default AA production script will be used.', action='store_true')

    info = parser.add_argument_group()
    info.add_argument('-V', '--Version', action='version', version='Version '+FuncConv.__version__)
    info.add_argument('-h', '--help',    action='help',    help='Show this message and exit.')
    info.add_argument('-CF', '--configuration', help='Specifiy configuration file. For an example input please see keb721/CCD2MD on GitHub.', default=None)
    
    args, unknownargs = parser.parse_known_args()
    command_line      = np.array(sys.argv)

    if args.configuration != None:    
        args, command_line = FuncConv.read_configuration_file(args.configuration, args, command_line)
    if args.conc == None:
        args.conc = 0.15
    if args.num_cpus == None:
        args.num_cpus = 1
    if args.SMILES == None:
        args.SMILES = []
    if args.gmx == None:
        args.gmx = 'gmx'

    # Test for --run_MD flag where --MD is not set
    # ---------------------------------------------
    
    for flag in ['CG_energy_minimise', 'CG_equil', 'AA_energy_minimise', 'AA_equil']:
        if args.__dict__[flag] == False and args.__dict__['run_'+flag] == True:
            print('# WARNING: Flag --run_{} has been set but flag --{} has not. Adding default mdp file for tpr generation.'.format(flag, flag))

            args.__dict__[flag] = True

            command_line = np.append(command_line, ['--{}'.format(flag)])
            
    # Read input data
    # ---------------
    
    input_data, title, cryst = FuncConv.read_in(args.inputfile, args.title)
    
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
    
    if args.membrane and args.ligchain:
        print('# WARNING: Ligand chains are currently incompatible with membrane insertion. Ligands will be appended to protein.')
        args.ligchain = False
    
    
    # Find residues to reorder
    # ------------------------
    
    ligands, atoms, IDs, types, locs = FuncConv.get_residues(input_data, 'CCD', SMILES, args.ligchain)

    # Reorder residues
    # ----------------
            
    output_data = input_data.to_dict('records')

    user_CCD_converted = []
    lig_CHARMM         = []
    
    for i, lig in enumerate(ligands):
        if types[i] != 'CHARMM':
            output_residues = FuncConv.convert_atomistic(lig, atoms[i], 'CCD', args.ligchain)
            lig_CHARMM.append(output_residues[0]['resnm'])
            output_data[IDs[i][0]:IDs[i][-1]+1] = output_residues
        else:
            if len(np.where(np.array(user_CCD_converted)==lig)[0])==0:
                nm = lig[:-5] + ' (userCCD)'
                print('# INFO: Gathering database information for molecule name {}'.format(nm))
                user_CCD_converted.append(lig)
                lig_CHARMM.append(lig)
    
    # Write output
    # -------------
    
    if args.outputfile.count('/') != 0:
        output_dir = '/'.join(args.outputfile.split('/')[:-1])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = ''    
            
    basename = '.'.join(args.outputfile.split('.')[:-1])
    
    PDBfile = basename+'_nomem.pdb' if args.membrane else args.outputfile
    FuncConv.write_PDB(PDBfile, output_data, title=title, cryst=cryst, ligand_chains=args.ligchain)
    
    # =================== #
    # Embed into membrane #
    # =================== # 
    
    if args.membrane:
        # Convert system to CG
        # --------------------
    
        if len(IDs) == 0:
            # Protein only
            prot = True
        else:
            prot = True if len(np.concatenate([np.array(i) for i in IDs])) != len(input_data) else False


        # Allow for warning flags
        
        mart_extras = FuncConv.get_CG_params(np.array(command_line), args.martinize, args.elastic, False) if prot else []
        check_mart  = np.array(mart_extras)
        
        if len(np.where(check_mart=='-maxwarn')[0])==0:
            # Check for secondary structure warning
            if len(np.where(check_mart=='-dssp')[0])==0 and len(np.where(check_mart=='-ss')[0])==0:
                # No additional secondary structure specification
            
                mart_v = subprocess.check_output(['martinize2', '-V'], universal_newlines = True)
                mart_v = mart_v.split()[-1].split('.')
        
                if int(mart_v[0]) > 0 or int(mart_v[1]) >= 15:
                    # Warning for secondary structure introduced in martinize 0.15.0
                    # Ignore this warning
                    print('# INFO: Ignoring martinize2 secondary structure prediction warning.')
                    mart_extras = np.hstack((mart_extras, np.array(['-maxwarn', '1'])))
            
                
        output_data = FuncConv.to_CG(args.inputfile, args.outputfile, input_data, mart_extras, ligands, input_data,
                                 types, locs, prot, newlipidome=args.newlipidome)
    
        CG_output = basename+'_CG_system.pdb'
    
        FuncConv.write_PDB(CG_output, output_data, title=title, cryst=cryst, ligand_chains=False)
    
        # Embed in membrane
        # -----------------
        
        FuncConv.build_membrane_CG(ligands, CG_output, args.outputfile, command_line, args.mempro, args.memprod, args.membrane, args.conc, num_CPUs=args.num_cpus, newlipidome=args.newlipidome) 

        CG_name = 'CG-system.gro' # For subsequent conversion
        
        # Optionally, energy minimise and/or equilibrate CG system
        # ---------------------------------------------------------
        
        if args.CG_energy_minimise or args.CG_equil:

            if not args.CG_energy_minimise:
                # Check for just equilibration
                command_line = np.append(command_line, ['-CGEM', '-rCGEM'])
                args.CG_energy_minimise = True ; args.run_CG_energy_minimise = True
                print('# WARNING: CG energy minimisation will be performed in addition to CG equilibration. The default energy minimisation script and parameters will be used.')
            

            Insane_dir = '.'.join(args.outputfile.split('.')[:-1])+'_MemPrO/Rank_1/CG_System_rank_1/'
            
            # NOTE: This will change output_data downstream!
            # Update CG_name as needed.

            ndx = None
            
            if args.make_ndx:
                # Making an ndx file - using CG information
                subprocess.run([args.gmx, 'make_ndx', '-f', Insane_dir+CG_name, '-o', Insane_dir+'CG-system.ndx'])
                ndx = 'CG-system.ndx'
                
            # Generate topology
            subprocess.run(['scp', Insane_dir+'topol.top', Insane_dir+'MemPrO_topology.top'])
            FuncConv.get_topology_CG(args.outputfile, args.membrane, ligands, prot, args.inputfile, args.newlipidome, topoldir=Insane_dir)
                           
            # Generate and run CG energy minimisation file
            em_name = FuncConv.make_gmx(Insane_dir+CG_name, command_line, 'CGEM', args.gmx, ndx=ndx, run=True)
            CG_name = em_name
            if args.CG_equil:
                # Equilibrate
                if not args.run_CG_equil:
                    command_line = np.append(command_line, ['-rCGeq'])
                    args.CG_equil = True
                eq_name = FuncConv.make_gmx(Insane_dir+em_name, command_line, 'CGeq_lip', args.gmx, ndx=ndx, run=True)
                CG_name = eq_name
                            
        # Convert back to atomistic
        # -------------------------
    
        FuncConv.convert_membrane_at(output_data, basename, command_line, args.cg2at, CG_name, args.gmx,
                                     args.newlipidome)
           
        subprocess.run(['scp', basename+'_CG2AT/FINAL/final_cg2at_aligned.pdb', args.outputfile])
    
    if args.membrane or not args.pdb2gmx:
        final, topol = FuncConv.get_topology_atomistic(args.outputfile, args.membrane, args.gmx,
                                                       output_data=pd.DataFrame.from_dict(output_data))
    else:
        final, topol = FuncConv.get_topology_atomistic(args.outputfile, args.membrane, at_command=command_line)

    # NOTE: At this point in output directory

    # Perform steering
    # ----------------

    if args.steer_ref != None or args.steer or args.steer_name != None:
        output_name = '.'.join(args.outputfile.split('/')[-1].split('.')[:-1])
        full_outdir = os.getcwd()

        # Create architecture
        if os.path.exists(output_name+'_steer'):
            print('# WARNING: Directory {} exists. This will be overwritten.'.format(full_outdir+'/'+output_name+'_steer'))
        else:
            os.makedirs(output_name+'_steer')

        # Identify the ligands to steer
        st_ligs = set(lig_CHARMM) if args.steer_name is None else args.steer_name
        

        FuncConv.steer_reference(args.steer_ref, st_ligs, output_name+'_steer', args.gmx)


        # Read in final position (with Hs)
        final_H, final_title, final_cryst = FuncConv.read_PDB(final)
        final_Steer = final_H
        
        for ligand in st_ligs:
            steer_loc         = output_name+'_steer/steer_{}/'.format(ligand)
            to_steer, lig_IDs = FuncConv.get_ligand(final_H, ligand, steer_loc)
            final_Steer       = FuncConv.steer_ligand(ligand, steer_loc, to_steer, lig_IDs, final_Steer, args.gmx)
            
        FuncConv.write_PDB(output_name+'_steered.pdb', final_Steer, title=final_title, cryst=final_cryst,
                           ligand_chains=args.ligchain)            
            
        final = output_name+'_steered.pdb'
        
    ndx = None
            
    if args.make_ndx:
        # Making an ndx file - using AA information
        ndx = '.'.join(final.split('.')[:-1])+'_ndx.ndx'
        subprocess.run([args.gmx, 'make_ndx', '-f', final, '-o', ndx])

        
    # For globular proteins, add to box if specified OR performing MD
    # ----------------------------------------------------------------

    make_box  = False
    run_AA_MD = False
    if args.AA_energy_minimise or args.AA_equil or args.AA_prod:
        run_AA_MD = True
        if not args.membrane:
            make_box = True
        
        
    if args.box or args.solvate or args.ions or make_box:
        if args.membrane:
            print('ERROR: Adding a box is not applicable for membrane-embedded proteins - this is done via Insane4MemPrO. Ignoring')
        else:
            # Test for box creation flags
            # ----------------------------

            print('# INFO: Creating solvation box.')
            
            for flag in ['box', 'solvate', 'ions']:
                if args.__dict__[flag] == False:
                    args.__dict__[flag] = True
                    command_line = np.append(command_line, ['--{}'.format(flag)])
            
            final = FuncConv.make_glob(final, command_line, args.gmx, 'AA', args.conc, topol=topol)

        
    # Optionally, energy minimise and/or equilibrate AA system and/or make production tpr
    # ------------------------------------------------------------------------------------
        
    if run_AA_MD:
        if not args.membrane and (args.AA_equil or args.AA_prod):
            print('# WARNING: Currently ccd2at has not created a box for the converted system. This is likely to cause issues for MD simulations other than energy minimisation.')
        if not args.AA_energy_minimise:
            # Check for just equilibration/production
            command_line = np.append(command_line, ['-AAEM', '-rAAEM'])
            args.run_AA_energy_minimise = True ; args.AA_energy_minimise = True
            if not args.AA_equil:
                # Just production
                print('# WARNING: AA energy minimisation will be performed in addition to generating AA production. Please note that equilibration is recommended. The default energy minimisation script and parameters will be used.')
            else:
                print('# WARNING: AA energy minimisation will be performed in addition to equilibration and generation of production tpr. The default energy minimisation script and parameters will be used.')
                

        output_dir = '/'.join(args.outputfile.split('/')[:-1])
        suffix     = '_lip' if args.membrane else '_prot'
        
        # Generate and run AA energy minimisation file
        
        em_name = FuncConv.make_gmx(final, command_line, 'AAEM', args.gmx, ndx=ndx, run=args.run_AA_energy_minimise, topol=topol)
        AA_name = em_name

        if args.AA_equil:
            # Equilibrate
            if args.AA_prod and not args.run_AA_equil:
                args.run_AA_equil = True
                command_line = np.append(command_line, ['-rAAeq'])
                
            eq_name = FuncConv.make_gmx(em_name, command_line, 'AAeq'+suffix, args.gmx, ndx=ndx, run=args.run_AA_equil, topol=topol)
            AA_name = eq_name


        if args.AA_prod:
            # Generate production data
            _ = FuncConv.make_gmx(AA_name, command_line, 'AAeq'+suffix, args.gmx, ndx=ndx, run=False, topol=topol)

            
	
	
if __name__=="__main__":
    main()
