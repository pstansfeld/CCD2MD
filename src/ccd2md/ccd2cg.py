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
import numpy as np
from copy import deepcopy
import sys, subprocess

def main():

    # Get command line arguments
    # ---------------------------
    
    parser = argparse.ArgumentParser(description='Convert from output from co-folding (specified CCD codes and SMILES strings and userCCD codes) to Martini 3 CG representation. Please see keb721/CCD2MD on GitHub for a list of the allowed SMILES strings and CCD codes.', add_help=False)
    
    parser.add_argument('inputfile',  help='Input file name - .cif, .pdb or .gro.', nargs='?')
    parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.', nargs='?')

    parser.add_argument('-T',  '--title',         help='Optionally specify the title of the generated pdb file. Default is the input PDB title (if present), or the name of the input file.', default=None)
    parser.add_argument('-S',  '--SMILES',        help='Used SMILES strings, list the order of the name of the ligand used. Note that when multiple of the same ligand are used this can be written either e.g. "POES POES" or "POES 2".', default=None)
    parser.add_argument('-nl', '--newlipidome',   help='Use updated Martini 3 mappings for lipids as in DOI: 10.1021/acscentsci.5c00755. Default off.', nargs='?', const=True, default=None)
    parser.add_argument('-g', '--gmx', help='Path to GROMACS exectuable - this will also be passed to CG2AT-lite. Default = "gmx".')

    parser.add_argument('-C',  '--conc',  help='Concentration of ions in the system - default = 0.15. Default ions are Na and Cl. For membrane embedded proteins, charge balance is maintained; specify different ions via -mem to pass to Insane4MemPrO. For non-membrane proteins charge balance is not maintained; specify different ions via -I/--ions.', default=None)
    
    CG_opts = parser.add_argument_group('martinize2 options')
    
    # Add user-extensible martinize2 options
    protein_network = CG_opts.add_mutually_exclusive_group()
    protein_network.add_argument('-E', '--elastic', help='Use the elastic network for Martini when converting an atomisitic protein to coarse-grained. When not present, or only the flag is present only -elastic is passed to martinize2. Include any desired elastic network commands to be passed to martinize2 after this argument (this may but does not need to include -elastic).', action='store_true')
    protein_network.add_argument('-G', '--go',      help='Use the Go network for Martini when converting an atomistic protein to coarse-grained. Include any required/desired commands to be passed to martinize2 after this argument.', action='store_true')
    
    CG_opts.add_argument('-M', '--martinize',       help='Add additional arguments to martinize2.', action='store_true')
    
    
    sys_opts = parser.add_argument_group('membrane-embedded system options')
    
    sys_opts.add_argument('-mem', '--membrane',   help='Embed the converted system into a membrane. Note that this will lead to minor rearragements of any ligands. After this flag, it is possible to add the arguments to be passed to Insane4MemPrO - most notably specifying the composition of the upper and lower leaflet in the form "-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1" where the codes represent lipids and the numbers represent the ratio between them. The default is two leaflets of pure POPC.', action='store_true')
    sys_opts.add_argument('-mp',  '--mempro',     help='Additional arguments for embedding the protein in the membrane using MemPrO. Add any additional arguments after this flag - default 5 grid points and 15 minimisation operations.', action='store_true')
    sys_opts.add_argument('-mdef', '--memprod',   help='Additional arguments for calculating the deformation of the membrane around the protein using MemPrOD. If not included, no deformations will be calculated. Otherwise any additional parameters for MemPrOD may be passed after this flag - otherwise MemPrOD defaults are used.', action='store_true')
    sys_opts.add_argument('-ncpu', '--num_cpus',  help='Number of CPUs to use for membrane embedding. Default = 1.', default=None, type=int)


    glob_opts = parser.add_argument_group('non-membrane protein system')
    glob_opts.add_argument('-B',  '--box',      help='Create a box for a non-membrane protein. The presence of this flag alone creates a box, solvates, and (if a non-zero concentration is specified) adds ions. After this flag, it is possible to add the arguments to be passed to gmx editconf. The default is "-c -d 2".', action='store_true')
    glob_opts.add_argument('-SV', '--solvate',  help='Solvate a box for a non-membrane protein. The presence of this flag alone creates a box, solvates, and (if a non-zero concentration is specified) adds ions. After this flag, it is possible to add the arguments to be passed to gmx solvate. The default is "-cs mdp_files/martini_water.gro -radius 0.21"; note that these are removed if any other flag is passed.', action='store_true')
    glob_opts.add_argument('-I',  '--ions',     help='Add ions for a non-membrane protein. The presence of this flag alone creates a box, solvates, and (if a non-zero concentration is specified) adds ions. After this flag, it is possible to add the arguments to be passed to gmx genion. The default is "-neutral -c {conc}" for conc specified via -C/--conc.', action='store_true')

    
    cgmd_opts = parser.add_argument_group('MD options')
    cgmd_opts.add_argument('-ndx', '--make_ndx', help='Make a gromacs index file either at the end of the CCD2MD run, or before attempting to create an equilibration tpr file. Please note that no additional arguments can be added after this command and index file creation is currently required to be interactive. Please also note that index file creation may be necessary for correct equilibration of the system.', action='store_true')
    cgmd_opts.add_argument('-CGEM', '--CG_energy_minimise',      help='Generate an energy minimisation for the final CG system. Additional arguments for grompp may be added after this flag - if none are added the default energy minimisation script will be used.', action='store_true')
    cgmd_opts.add_argument('-rCGEM', '--run_CG_energy_minimise', help='Run energy minimisation for the final CG system (will also happen if equilibration or production steps are added). Additional arguments for mdrun may be added after this flag.', action='store_true')
    cgmd_opts.add_argument('-CGeq', '--CG_equil',                help='Generate an equilibration for the final CG system. Additional arguments for grompp may be added after this flag - if none are added the default CG equilibration script will be used.', action='store_true')
    cgmd_opts.add_argument('-rCGeq', '--run_CG_equil',           help='Add additional arguments for running CG equilibration via gmx mdrun (will be run if production tpr generation specified).', action='store_true')
    cgmd_opts.add_argument('-CGprd', '--CG_prod',                help='Generate a production tpr file for the final CG system. Additional arguments for grompp may be added after this flag - if none are added the default CG production script will be used.', action='store_true')

    info = parser.add_argument_group()
    info.add_argument('-V', '--Version', action='version', version='Version '+FuncConv.__version__)
    info.add_argument('-h', '--help',    action='help',    help='Show this message and exit.')
    info.add_argument('-CF', '--configuration', help='Specifiy configuration file. For an example input please see keb721/CCD2MD on GitHub.', default=None)
    
    
    args, unknownargs = parser.parse_known_args()
    command_line      = np.array(sys.argv)
    

    if args.configuration != None:
        args, command_line = FuncConv.read_configuration_file(args.configuration, args, command_line)
        
    # Restore defaults - not set to determine command line presence
        
    if args.conc == None:
        args.conc = 0.15
    if args.num_cpus == None:
        args.num_cpus = 1
    if args.SMILES == None:
        args.SMILES = []

        
    # Test for --run_MD flag where --MD is not set
    # ---------------------------------------------
    
    for flag in ['CG_energy_minimise', 'CG_equil']:
        if args.__dict__[flag] == False and args.__dict__['run_'+flag] == True:
            print('# WARNING: flag --run_{} has been set but flag --{} has not. Adding default mdp file for tpr generation.'.format(flag, flag))

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
    
    # Find residues to reorder
    # ------------------------
        
    ligands, atoms, IDs, types, locs = FuncConv.get_residues(input_data, 'CCD', SMILES, False)
    
    # Convert system to CG
    if len(IDs) == 0:
        # Protein only
        prot = True
    else:
        prot = True if len(np.concatenate([np.array(i) for i in IDs])) != len(input_data) else False    
    
    # Generate the martinize2 parameters
    mart_params = FuncConv.get_CG_params(np.array(command_line), args.martinize, args.elastic, args.go) if prot else []
    
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
        FuncConv.build_membrane_CG(ligands, PDBfile, args.outputfile, command_line, args.mempro, args.memprod,
                                   args.membrane, args.conc, newlipidome=args.newlipidome, num_CPUs=args.num_cpus)
        mempro_output, title, cryst = FuncConv.read_GRO(basename+'_MemPrO/Rank_1/CG_System_rank_1/CG-system.gro')
        
        FuncConv.write_PDB(args.outputfile, mempro_output, title=title, cryst=cryst, ligand_chains=False)
    
    FuncConv.get_topology_CG(args.outputfile, args.membrane, ligands, prot, args.inputfile, newlipidome=args.newlipidome)

    # For globular proteins, add to box if specified OR performing MD
    # ----------------------------------------------------------------
    
    run_CG_MD = False
    make_box  = False
    if args.CG_energy_minimise or args.CG_equil or args.CG_prod:
        run_CG_MD = True
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
            
            final = FuncConv.make_glob(final, command_line, args.gmx, 'CG', args.conc, topol=topol)
    
    ndx = None
        
    if args.make_ndx:
        # Making an ndx file - using CG information
        subprocess.run([args.gmx, 'make_ndx', '-f', Insane_dir+CG_name, '-o', Insane_dir+'CG-system.ndx'])
        ndx = 'CG-system.ndx'

    # Optionally, energy minimise and/or equilibrate CG system and/or make production tpr
    # ------------------------------------------------------------------------------------
        
    if run_CG_MD:
        if not args.CG_energy_minimise:
            # Check for just equilibration/production
            command_line = np.append(command_line, ['-CGEM', '-rCGEM'])
            args.run_CG_energy_minimise = True ; args.CG_energy_minimise = True
            if not args.CG_equil:
                # Just production
                print('# WARNING: CG energy minimisation will be performed in addition to generating CG production. Please note that equilibration is recommended. The default energy minimisation script and parameters will be used.')
            else:
                print('# WARNING: CG energy minimisation will be performed in addition to CG equilibration and generation of production tpr. The default energy minimisation script and parameters will be used.')
                

        output_dir = '/'.join(args.outputfile.split('/')[:-1])
        suffix     = '_lip' if args.membrane else '_prot'
        
        # Generate and run AA energy minimisation file

        if args.CG_equil or args.CG_prod:
            command_line = np.append(command_line, ['-rCGEM'])
            args.run_CG_energy_minimise = True
        
        em_name = FuncConv.make_gmx(final, command_line, 'CGEM', args.gmx, ndx=ndx,
                                    run=args.run_AA_energy_minimise, topol=topol)
        CG_name = em_name

        if args.CG_equil:
            # Equilibrate
            if args.CG_prod and not args.run_CG_equil:
                args.run_CG_equil = True
                command_line = np.append(command_line, ['-rCGeq'])
            
            eq_name = FuncConv.make_gmx(em_name, command_line, 'CGeq'+suffix, args.gmx,
                                        ndx=ndx, run=args.run_AA_equil, topol=topol)
            CG_name = eq_name


        if args.CG_prod:
            # Generate production data
            _ = FuncConv.make_gmx(CG_name, command_line, 'CGeq'+suffix, args.gmx, ndx=ndx, run=False, topol=topol)


    
if __name__ == "__main__":
    main()
