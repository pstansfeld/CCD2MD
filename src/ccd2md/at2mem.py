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
    
    parser = argparse.ArgumentParser(description='Embed a CHARMM system in a membrane - please note that this requires coarse-graining. Please see keb721/CCD2AT on github for a list of the allowed (CG-mapped) molecules.', add_help=False)
    
    parser.add_argument('inputfile',  help='Input file name - .cif, .pdb or .gro.', nargs='?')
    parser.add_argument('outputfile', help='Output file name - will be written in .pdb format.', nargs='?')

    parser.add_argument('-g', '--gmx', help='Path to GROMACS exectuable - this will also be passed to CG2AT-lite. Default = "gmx".')
    parser.add_argument('-T',  '--title',     help='Optionally specify the title of the generated pdb file.', default = None)
    parser.add_argument('-ndx', '--make_ndx', help='Make a gromacs index file either at the end of the CCD2MD run, or before attempting to create an equilibration tpr file. Please note that no additional arguments can be added after this command and index file creation is currently required to be interactive. Please also note that index file creation may be necessary for correct equilibration of the system.', action='store_true')
    
    sys_opts = parser.add_argument_group('membrane-embedded system options')
    
    sys_opts.add_argument('-mem', '--membrane',   help='Embed the converted system into a membrane. Note that this will lead to minor rearragements of any ligands. After this flag, it is possible to add the arguments to be passed to Insane4MemPrO - most notably specifying the composition of the upper and lower leaflet in the form "-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1" where the codes represent lipids and the numbers represent the ratio between them. The default is two leaflets of pure POPC.', action='store_true')
    sys_opts.add_argument('-C',  '--conc',        help='Concentration of NaCl in system - charge balance is maintained, Default = 0.15.', default=None)
    sys_opts.add_argument('-mp',  '--mempro',     help='Additional arguments for embedding the protein in the membrane using MemPrO. Add any additional arguments after this flag - default 5 grid points and 15 minimisation operations.', action='store_true')
    sys_opts.add_argument('-mdef', '--memprod',   help='Additional arguments for calculating the deformation of the membrane around the protein using MemPrOD. If not included, no deformations will be calculated. Otherwise any additional parameters for MemPrOD may be passed after this flag - otherwise MemPrOD defaults are used.', action='store_true')
    sys_opts.add_argument('-ncpu', '--num_cpus',  help='Number of CPUs to use for membrane embedding. Default = 1.', default=None)
    sys_opts.add_argument('-at', '--cg2at',       help='Additional arguments to pass to CG2AT. Add additional arguments after this flag.', action='store_true')
    # sys_opts.add_argument('-sol', '--water',    help='The water model to use - currently only TIP3P') 



    cgmd_opts = parser.add_argument_group('coarse-grained membrane-embedded system MD options')
    cgmd_opts.add_argument('-nl', '--newlipidome',               help='Use updated Martini 3 mappings for CG lipids as in DOI: 10.1021/acscentsci.5c00755. Default off.', action='store_true')
    cgmd_opts.add_argument('-CGEM', '--CG_energy_minimise',      help='Generate an energy minimisation for the CG system before converting back to atomistic. Additional arguments for grompp may be added after this flag - if none are added the default energy minimisation script will be used.', action='store_true')
    cgmd_opts.add_argument('-rCGEM', '--run_CG_energy_minimise', help='Add additional arguments for running CG energy minimistaion (via gmx mdrun) before conversion back to atomstic.', action='store_true')
    cgmd_opts.add_argument('-CGeq', '--CG_equil',                help='Generate an equilibration for the CG system before converting back to atomistic. Additional arguments for grompp may be added after this flag - if none are added the default CG equilibration script will be used.', action='store_true')
    cgmd_opts.add_argument('-rCGeq', '--run_CG_equil',           help='Add additional arguments for running CG equilibration (via gmx mdrun) before conversion back to atomstic.', action='store_true')

    
    CG_opts = parser.add_argument_group('martinize2 options - note only relevant for CG MD')    
    # Add user-extensible martinize2 options
    protein_network = CG_opts.add_mutually_exclusive_group()
    protein_network.add_argument('-E', '--elastic', help='Use the elastic network for Martini when converting an atomisitic protein to coarse-grained. When not present, or only the flag is present only -elastic is passed to martinize2. Include any desired elastic network commands to be passed to martinize2 after this argument (this may but does not need to include -elastic).', action='store_true')
    CG_opts.add_argument('-M', '--martinize',       help='Add additional arguments to martinize2.', action='store_true')

    
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

    # Test for --run_MD flag where --MD is not set
    # ---------------------------------------------

    for flag in ['CG_energy_minimise', 'CG_equil', 'AA_energy_minimise', 'AA_equil']:
        if args.__dict__[flag] == False and args.__dict__['run_'+flag] == True:
            print('# WARNING: flag --run_{} has been set but flag --{} has not. Adding default mdp file for tpr generation.'.format(flag, flag))

            args.__dict__[flag] = True

            command_line = np.append(command_line, ['--{}'.format(flag)])

    
    # Read input data
    # ---------------

    input_data, title, cryst = FuncConv.read_in(args.inputfile, args.title)
    
    # Find ligands to omit
    # ---------------------
        
    ligands, atoms, IDs, types, locs = FuncConv.get_residues(input_data, 'CHARMM', [], False)
    
    output_data = input_data.to_dict('records')
    
    
    # Write output
    # -------------
    
    if args.outputfile.count('/') != 0:
        output_dir = '/'.join(args.outputfile.split('/')[:-1])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
                
    
    basename = '.'.join(args.outputfile.split('.')[:-1])
    
    PDBfile = basename+'_nomem.pdb'
    FuncConv.write_PDB(PDBfile, output_data, title=title, cryst=cryst, ligand_chains=False)
    
    # =================== #
    # Embed into membrane #
    # =================== # 
    
    # Convert system to CG
    # --------------------

    if len(IDs) == 0:
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
                mart_extras.extend(['-maxwarn', '1'])
            
            
    output_data = FuncConv.to_CG(args.inputfile, args.outputfile, input_data, mart_extras, ligands, input_data,
                                 types, locs, prot)
    
    CG_output = basename+'_CG_system.pdb'
    
    FuncConv.write_PDB(CG_output, output_data, title=title, cryst=cryst, ligand_chains=False)
    
    # Embed in membrane
    # -----------------
        
    FuncConv.build_membrane_CG(ligands, CG_output, args.outputfile, command_line, args.mempro, args.memprod, args.membrane, args.conc, num_CPUs=args.num_cpus) 

    CG_name = 'CG-system.gro' # For subsequent conversion


    # Optionally, energy minimise and/or equilibrate CG system
    # ---------------------------------------------------------

    if args.CG_energy_minimise or args.CG_equil:

        if not args.CG_energy_minimise:
            # Check for just equilibration
            command_line = np.append(command_line, ['-CGEM', '-rCGEM'])
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
            eq_name = FuncConv.make_gmx(Insane_dir+em_name, command_line, 'CGeq_lip', args.gmx, ndx=ndx, run=True)
            CG_name = eq_name
    
    # Convert back to atomistic
    # -------------------------
        
    FuncConv.convert_membrane_at(output_data, basename, command_line, args.cg2at, args.newlipidome)
    
    subprocess.run(['scp', args.outputfile.split('.')[0]+'_CG2AT/FINAL/final_cg2at_aligned.pdb', args.outputfile])
    
    
    FuncConv.get_topology_atomistic(args.outputfile, True)

    ndx = None
            
    if args.make_ndx:
        # Making an ndx file - using AA information
        ndx = '.'.join(final.split('.')[:-1])+'_ndx.ndx'
        subprocess.run([args.gmx, 'make_ndx', '-f', final, '-o', ndx])

        
    # Optionally, energy minimise and/or equilibrate AA system and/or make production tpr
    # ------------------------------------------------------------------------------------
        
    if args.AA_energy_minimise or args.AA_equil or args.AA_prod:
        if not args.membrane and (args.AA_equil or args.AA_prod):
            print('# WARNING: Currently ccd2at has not created a box for the converted system. This is likely to cause issues for MD simulations other than energy minimisation.')
        if not args.AA_energy_minimise:
            # Check for just equilibration/production
            command_line = np.append(command_line, ['-AAEM', '-rAAEM'])
            args.run_AA_energy_minimise = True
            if not args.AA_equil:
                # Just production
                print('# WARNING: AA energy minimisation will be performed in addition to generating AA production. Please note that equilibration is recommended. The default energy minimisation script and parameters will be used.')
            else:
                print('# WARNING: AA energy minimisation will be performed in addition to equilibration and generation of production tpr. The default energy minimisation script and parameters will be used.')
                
        output_dir = '/'.join(args.outputfile.split('/')[:-1])
        suffix     = '_lip'
        
        # Generate and run AA energy minimisation file

        if args.AA_equil or args.AA_prod:
            command_line = np.append(command_line, ['-rAAEM'])
            args.run_AA_energy_minimise = True
        
        em_name = FuncConv.make_gmx(final, command_line, 'AAEM', args.gmx, ndx=ndx,
                                    run=args.run_AA_energy_minimise, topol=topol)
        AA_name = em_name

        if args.AA_equil:
            # Equilibrate
            if args.AA_prod and not args.run_AA_equil:
                args.run_AA_equil = True
                command_line = np.append(command_line, ['-rAAeq'])
            
            eq_name = FuncConv.make_gmx(em_name, command_line, 'AAeq'+suffix, args.gmx, ndx=ndx,
                                        run=args.run_AA_equil, topol=topol)
            AA_name = eq_name


        if args.AA_prod:
            # Generate production data
            _ = FuncConv.make_gmx(AA_name, command_line, 'AAeq'+suffix, args.gmx, ndx=ndx, run=False, topol=topol)


    

if __name__ == "__main__":
    main()
