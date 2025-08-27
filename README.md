# CCD2MD codes for conversion between cofolding outputs and simulations

For more information please contact Kat Blow (katarina.blow@warwick.ac.uk)

Ref_data/ contains reference data
FuncConv.py contains functions to enable conversion

Functions are converted using:
- ccd2at.py  (cofolding to CHARMM)
- at2ccd.py  (CHARMM to CCD/SMILES)
- ccd2cg.py  (cofolding to Martini 3, where possible)
- at2cg.py   (CHARMM to Martini 3, where possible)
- at2mem.py  (embed CHARMM structure in membrane)
- pos2cif.py (structure file(s) to userCCD)

Note that there are limitations for the ligands available for conversion to Martini 3

### Usage

pdb2cif is self contained, otherwise (currently) need to ensure that `Ref_data/` and `FuncConv.py` are in the directory of use. Dependencies are numpy and pandas with vermouth-martinize (martinize2) required for protein conversion, MemPrO required for protein membrane embedding and gromacs (specifically pdb2gmx) required for atomistic topology generation. Topology files are generated except for at2ccd.py and pos2cif.py

#### ccd2at

> python ccd2at.py INPUT_FILE OUTPUT_FILE [-lc] [-S &lt;smiles&gt; ...] [-mem [&lt;membrane&gt;]] [-C &lt;conc&gt;] [-mp &lt;mempro&gt;]

`INPUT_FILE`  name of the file to convert - the extension must be either .cif, .pdb or .gro

`OUTPUT_FILE` name of the converted file. This will be written as a pdb file - with an OUTPUT_FILE_H.pdb version with added hydrogens if not membrane embedded.

*Optional arguments*

`-lc/--ligchain` include ligands as their own chains - default off (binary switch)

`-S/--SMILES <smiles> ...`  indicates that some ligands have been input as SMILES strings - note that other than for POPE (see below) this requires user-defined ligand mapping. To use, list the **CCDName** of the ligands represented by SMILES strings in order of use (e.g. POES for a POPE smiles string, ideally should not overlap with CCD or CHARMM names). Note that when multiple of the same ligand are used this can be written either e.g. `-S POES POES` or `-S POES 2`. 

*Optional arguments -- membrane embedding*

`-mem/--membrane [<membrane>]` acts as a flag to embed the system in a membrane. Optionally, it may be followed by the membrane composition in the form `-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1` where the numbers give the ratio of lipids, `-u` represents the upper leaflet and `-l` represents the lower leaflet. If not specified, both leaflets are pure POPC. Note that membrane embedding will lead to minor rearrangement of bound ligands.

`-C/--conc <conc>` passes a single number giving the concentration of NaCl in the system, where the charge is balanced. 

`-mp/--mempro <mempro>` passes optional arguments to [MemPrO](https://github.com/ShufflerBardOnTheEdge/MemPrO). Default is 5 grid points and 15 minimisation operations 


#### at2ccd

> python at2ccd.py INPUT_FILE OUTPUT_FILE [-lc] 


`INPUT_FILE`  name of the file to convert - the extension must be either .cif, .pdb or .gro

`OUTPUT_FILE` name of the converted file. This will be written as a pdb file

*Optional arguments*

`-lc/--ligchain` include ligands as their own chains - default off (binary switch)


#### ccd2cg

> python ccd2cg.py INPUT_FILE OUTPUT_FILE [-S &lt;smiles&gt; ...] [-nl] [(-E [&lt;elastic&gt;]) | (-G &lt;go&gt;)] [-M &lt;martinize&gt;] [-mem [&lt;membrane&gt;]] [-C &lt;conc&gt;] [-mp &lt;mempro&gt;]

`INPUT_FILE`  name of the file to convert - the extension must be either .cif, .pdb or .gro

`OUTPUT_FILE` name of the converted file. This will be written as a pdb file

*Optional arguments*

`-S/--SMILES <smiles> ...`  indicates that some ligands have been input as SMILES strings - note that other than for POPE (see below) this requires user-defined ligand mapping. To use, list the **CCDName** of the ligands represented by SMILES strings in order of use (e.g. POES for a POPE smiles string, ideally should not overlap with CCD or CHARMM names). Note that when multiple of the same ligand are used this can be written either e.g. `-S POES POES` or `-S POES 2`. 

`-nl/--newlipidome` is a flag indicating that the new mappings for the Martini 3 lipidome (from [DOI: 10.1021/acscentsci.5c00755](https://doi.org/10.1021/acscentsci.5c00755)) should be used.

*Optional arguments -- protein coarse-graining*

`-E/--elastic [<elastic>]` specifies that protein secondary structure should be biased using an elastic network (default). Additional options for the elastic network in martinize may be added added this flag; if omitted the defaults will be used.

`-G/--go <go>` specifies that protein secondary structure should be biased using a go network. Required and optional parameters should be added after this flag.

`-M/--martinize <martinize>` is used for any additonal optional arguments which should be passed to martinize2 for protein CG conversion.


*Optional arguments -- membrane embedding*

`-mem/--membrane [<membrane>]` acts as a flag to embed the system in a membrane. Optionally, it may be followed by the membrane composition in the form `-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1` where the numbers give the ratio of lipids, `-u` represents the upper leaflet and `-l` represents the lower leaflet. If not specified, both leaflets are pure POPC.

`-C/--conc <conc>` passes a single number giving the concentration of NaCl in the system, where the charge is balanced.

`-mp/--mempro <mempro>` passes optional arguments to [MemPrO](https://github.com/ShufflerBardOnTheEdge/MemPrO). Default is 5 grid points and 15 minimisation operations 

#### at2cg

> python at2cg.py INPUT_FILE OUTPUT_FILE [(-E [&lt;elastic&gt;]) | (-G &lt;go&gt;)] [-M &lt;martinize&gt;] [-mem [&lt;membrane&gt;]] [-C &lt;conc&gt;] [-mp &lt;mempro&gt;]

`INPUT_FILE`  name of the file to convert - the extension must be either .cif, .pdb or .gro

`OUTPUT_FILE` name of the converted file. This will be written as a pdb file

*Optional arugments*

`-nl/--newlipidome` is a flag indicating that the new mappings for the Martini 3 lipidome (from [DOI: 10.1021/acscentsci.5c00755](https://doi.org/10.1021/acscentsci.5c00755) should be used.

*Optional arguments -- protein coarse-graining*

`-E/--elastic [<elastic>]` specifies that protein secondary structure should be biased using an elastic network (default). Additional options for the elastic network in martinize may be added added this flag; if omitted the defaults will be used.

`-G/--go <go>` specifies that protein secondary structure should be biased using a go network. Required and optional parameters should be added after this flag.

`-M/--martinize <martinize>` is used for any additonal optional arguments which should be passed to martinize2 for protein CG conversion.


*Optional arguments -- membrane embedding*

`-mem/--membrane [<membrane>]` acts as a flag to embed the system in a membrane. Optionally, it may be followed by the membrane composition in the form `-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1` where the numbers give the ratio of lipids, `-u` represents the upper leaflet and `-l` represents the lower leaflet. If not specified, both leaflets are pure POPC.

`-C/--conc <conc>` passes a single number giving the concentration of NaCl in the system, where the charge is balanced.

`-mp/--mempro <mempro>` passes optional arguments to [MemPrO](https://github.com/ShufflerBardOnTheEdge/MemPrO). Default is 5 grid points and 15 minimisation operations 


#### at2mem

> python at2mem.py INPUT_FILE OUTPUT_FILE [-lc]  [-M &lt;martinize&gt;] [-mem [&lt;membrane&gt;]] [-C &lt;conc&gt;] [-mp &lt;mempro&gt;]

`INPUT_FILE`  name of the file to convert - the extension must be either .cif, .pdb or .gro

`OUTPUT_FILE` name of the converted file. This will be written as a pdb file


*Optional arguments*

`-lc/--ligchain` include ligands as their own chains - default off (binary switch)


*Optional arguments -- membrane embedding*

`-mem/--membrane [<membrane>]` may be followed by the membrane composition in the form `-u POPE:7 -u POPG:2 -u CARD:1 -l POPE:7 -l POPG:2 -l CARD:1` where the numbers give the ratio of lipids, `-u` represents the upper leaflet and `-l` represents the lower leaflet. If not specified, both leaflets are pure POPC. Note that membrane embedding will lead to minor rearrangement of bound ligands.

`-C/--conc <conc>` passes a single number giving the concentration of NaCl in the system, where the charge is balanced.

`-mp/--mempro <mempro>` passes optional arguments to [MemPrO](https://github.com/ShufflerBardOnTheEdge/MemPrO). Default is 5 grid points and 15 minimisation operations 


#### pos2cif

pdb2cif allows for the creation of userCCD codes which can be input into AF3. When the atomic naming and orderings match that of the desired forcefield then no additional post-processing is required. Modified amino acids can also be input and the resulting CCD code manually added as a PTM for the correct location in the protein. Ligands composed of chains of multiple residues can be added, although the quality of the joins may be variable.

> python pos2cif.py [-n &lt;ligand&gt; ... -f &lt;files&gt; ...] [-r (&lt;old&gt; &lt;new&gt;) ...] [-c ((&lt;pdb&gt; &lt;itp&gt;) | &lt;mol2&gt;) ...] [-e &lt;charge&gt;] [-b &lt;bond&gt;] [-H] [-j &lt;json&gt;] [-t &lt;title&gt;] [-A &lt;afvers&gt;] [-s &lt;seed&gt; ...] [-d &lt;dialect&gt;] [-p (&lt;protein&gt; | (&lt;protein N &gt;)) ...]  

*Specifying ligands* 

`-n/--names <ligand> ...` gives the names of the ligands to extract from files and produce cif files for
 
`-f/--files <files> ...` gives the name of position (required) and bonding (optional) files for the ligands to be converted. Position files may be .pdb, .crd, .gro, or .mol2 with additional bonding information provided by .rtp, .rtf, or .itp files. Position files may contain information for several ligands, but should only include one copy of the ligand of interest.

`-c/--covalent ((<pdb> <itp>) | <mol2>) ...` gives either a .pdb and .itp file (note, that this must follow the .pdb then .itp order) or a .mol2 file for every covalently bonded ligand to be added to the system. Note that in contrast to the files for single component ligands, these files should include information for the covalently bonded ligand only.

*Optional arguments*

`-r/--rename (<old> <new>) ...` allows ligands to have a different userCCD name than in the original file. For each ligand to be renamed, `old` is the name in the position file, and `new` is the desired output name

`-e/--charge <charge>` is an optional cutoff for increments of the integer charge in the cif files. For each partial integer charge above the cutoff the integer charge output will be changed by 1. The default is 0.75 $e$.

`-b/--bond <bond>` is an optional cutoff for generating bonding information from positional information -- note that this is not utilised where bonding information is given in a bonding file. If the distance between two atoms is less than the cutoff, a bond between these atoms is output. The default is 1.4 \AA.

`-H/--Hydrogens` is a flag to output hydrogen position and bonding information.

*json file parameters* 

`-j/--json <json>` gives the option to change the name of the json file containing the required userCCD field for AlphaFold3 inputs (default output.json). For each ligand, an additional file containing the cif data is output to NAME_output.cif)

`-t/--title <title>` gives the name of the system to be passed to the output json file. The default is pos2cif_system

`-A/--afvers <version>` gives the version number to be output in the json file. The default is 2.

`-s/--seeds <model seed> ...` gives the model seeds to pass to AF3 - these need not be comma separated. The default seed is 1.

`-d/--dialect <dialect>` gives the dialect to be output in the json file. The default is alphafold3.

`-p/--protein (<protein> | (<protein> <N>)) ...` allows input of FASTA protein sequence(s) for the json file. Sequences can be put in in FASTA form, numbers can be inserted to detail how many of the sequence should be added.


### Current ligands
 
| Name                       | CCDName | CHARMMName | CGName | Note                                 |   
| -------------------------- | ------- | ---------- | ------ | ------------------------------------ |
| Cardiolipin                |         |    CARD    |  CARD  | Palmitoyl-oleoyl cardiolipin (POCL2) |
| Cholesterol                |   CLR   |    CHL1    |  CHL1  |                                      |
| DGPC                       |         |    DGPC    |  DGPC  |                                      |
| DGPE                       |         |    DGPE    |  DGPE  |                                      |
| DLPC                       |         |    DLPC    |  DLPC  |                                      |
| DLPE                       |   WNZ   |    DLPE    |  DLPE  |                                      |
| DMPC                       |   MC3   |    DMPC    |  DMPC  |                                      |
| DMPE                       |   46E   |    DMPE    |  DMPE  |                                      |
| DMPI                       |         |    DMPI    |  DMPI  |                                      |
| DMPS                       |         |    DMPS    |  DMPS  |                                      |
| DNPC                       |         |    DNPC    |  DNPC  |                                      |
| DNPE                       |         |    DNPE    |  DNPE  |                                      |
| dodecylBDmaltoside         |   LMT   |    BDDM    |        |                                      |
| DOPC                       |         |    DOPC    |  DOPC  |                                      |
| DOPE                       |         |    DOPE    |  DOPE  |                                      |
| DOPS                       |   17F   |    DOPS    |  DOPS  |                                      |
| DPG3                       |         |    DPG3    |  DPG3  | Non-native in CHARMM36               |
| DPPC				         |   PCF   |    DPPC    |  DPPC  |                                      |
| DPPE				         |         |    DPPE    |  DPPE  |                                      |
| DYPC				         |         |    DYPC    |  DYPC  |                                      |
| DYPE				         |         |    DYPE    |  DYPE  |                                      |
| laurylBMNglycol            |   LMN   |    BLMN    |        | BLMNG in CHARMM                      |
| LIP1                       |         |    LIP1    |  LIP1  | Non-native in CHARMM36               |
| LIP2                       |         |    LIP2    |  LIP2  | Non-native in CHARMM36               |
| LIP3                       |         |    LIP3    |  LIP3  | Non-native in CHARMM36               |
| LIPA                       |         |    LIPA    |  LIPA  | Non-native in CHARMM36               |
| OBDglucopyranoside         |   BOG   |    BOG1    |        | BOG in CHARMM36                      |
| POPI 3 phosphate           |         |    POP1	|  POP1  | POPI13 in CHARMM36				    |
| POPI 3,4 bisphosphate      |         |    POP2	|  POP2  | POPI2D in CHARMM36				    |
| POPI 3,4,5 trisphosphate   |         |    POP3	|  POP3  | POPI34 in CHARMM36				    |
| POPI 4 phosphate           |         |    POP4	|  POP4  | POPI14 in CHARMM36				    |
| POPI 5 phosphate           |         |    POP5	|  POP5  | POPI15 in CHARMM36				    |
| POPI 4,5 bisphosphate      |         |    POP6	|  POP6  | POPI24 in CHARMM36				    |
| POPI 3,5 bisphosphate      |         |    POP7	|  POP7  | POPI2A in CHARMM36				    |	
| POPA                       |   D21   |    POPA    |  POPA  |                                      |
| POPC                       |   POV   |    POPC    |  POPC  |                                      |
| POPE                       |   PEV   |    POPE    |  POPE  |                                      |
| POPE_SMILES                |   POES  |    POPE    |  POPE  | From SMILES string                   |
| POPG                       |   PGW   |    POPG    |  POPG  |                                      |
| POPI   				     |         |    POPI    |  POPI  |                                      |
| POPS                       |   D39   |    POPS    |  POPS  |                                      |
| RAMP                       |         |    RAMP    |  RAMP  | Non-native in CHARMM36               |
| Kdo2-lipidA                |   KDL   |    REMP    |  REMP  |                                      |
| sphingomyelin              |         |    SSM1    |  SSM1  | SSM in CHARMM36                      |
| TMM                        |         |    TMM1    |  TMM1  | Non-native in CHARMM36               |
| TMMA                       |         |    TMMA    |  TMMA  | Non-native in CHARMM36               |
| undecaprenyl phosphate     |   5TR   |    UDP1    |  UDP1  |                                      |
| undecaprenyl pyrophosphate |         |    UDP2    |  UDP2  | UNDPP in CHARMM36                    |


### SMILES strings

Note any modification which does not affect the order of atoms can be made (e.g. changes to chiral centres/charges/double bonds). The name of the SMILES string must differ from the CHARMM name if wishing to utilise userCCD codes, and should not be more than 4 charaters if wishing to convert to CG. Ideally, it should also differ from utilised CCD codes.

| POPE_SMILES                |   POES  |  POPE  |```CCCCCCCC\C=C/CCCCCCCC(=O)O[C@H](COC(=O)CCCCCCCCCCCCCCC)COP(O)(=O)OCCN``` |


For issues and suggestions, please feel free to add to the GitHub or contact Kat Blow at katarina.blow[at]warwick.ac.uk