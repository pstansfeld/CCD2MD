# Current ligands

A list of the ligands currently natively available for conversion within CCD2MD. For addition of other molecules, the use of [pos2cif](https://github.com/keb721/CCD2MD/blob/main/docs/pos2cif_docs.md) is recommended if AF3 can be used. Otherwise, additional mappings can be created, as documented in the [CCD2MD paper](https://doi.org/10.1021/acs.jcim.5c02066).

| Name                       | CCDName | CHARMMName | CGName | Note                                                                 |  
| -------------------------- | ------- | ---------- | ------ | -------------------------------------------------------------------- |
| ADP                        |   ADP   |     ADP    |        |                                                                      |
| ATP                        |   ATP   |     ATP    |        |                                                                      |
| Cardiolipin                |         |    CARD    |  CARD  | Palmitoyl-oleoyl cardiolipin (POCL2), new lipidome mapping available |
| Cholesterol                |   CLR   |    CHL1    |  CHL1  | New lipidome mapping available                                       |
| DGPC                       |         |    DGPC    |  DGPC  | New lipidome mapping available                                       |
| DGPE                       |         |    DGPE    |  DGPE  | New lipidome mapping available                                       |
| DLPC                       |         |    DLPC    |  DLPC  | New lipidome mapping available                                       |
| DLPE                       |   WNZ   |    DLPE    |  DLPE  | New lipidome mapping available                                       |
| DMPC                       |   MC3   |    DMPC    |  DMPC  | New lipidome mapping available                                       |
| DMPE                       |   46E   |    DMPE    |  DMPE  | New lipidome mapping available                                       |
| DMPI                       |         |    DMPI    |  DMPI  | New lipidome mapping available                                       |
| DMPS                       |         |    DMPS    |  DMPS  | New lipidome mapping available                                       |
| DNPC                       |         |    DNPC    |  DNPC  | New lipidome mapping available                                       |
| DNPE                       |         |    DNPE    |  DNPE  | New lipidome mapping available                                       |
| dodecylBDmaltoside         |   LMT   |    BDDM    |        |                                                                      |
| DOPC                       |         |    DOPC    |  DOPC  | New lipidome mapping available                                       |
| DOPE                       |         |    DOPE    |  DOPE  | New lipidome mapping available                                       |
| DOPS                       |   17F   |    DOPS    |  DOPS  | New lipidome mapping available                                       |
| DPG3                       |         |    DPG3    |  DPG3  | Non-native in CHARMM36                                               |
| DPPC				         |   PCF   |    DPPC    |  DPPC  | New lipidome mapping available                                       |
| DPPE				         |         |    DPPE    |  DPPE  | New lipidome mapping available                                       |
| DYPC				         |         |    DYPC    |  DYPC  | New lipidome mapping available                                       |
| DYPE				         |         |    DYPE    |  DYPE  | New lipidome mapping available                                       |
| laurylBMNglycol            |   LMN   |    BLMN    |        | BLMNG in CHARMM                                                      |
| LIP1                       |         |    LIP1    |  LIP1  | Non-native in CHARMM36                                               |
| LIP2                       |         |    LIP2    |  LIP2  | Non-native in CHARMM36                                               |
| LIP3                       |         |    LIP3    |  LIP3  | Non-native in CHARMM36                                               |
| LIPA                       |         |    LIPA    |  LIPA  | Non-native in CHARMM36                                               |
| OBDglucopyranoside         |   BOG   |    BOG1    |        | BOG in CHARMM36                                                      |
| POPI 3 phosphate           |         |    POP1	|  POP1  | POPI13 in CHARMM36, new lipidome mapping available        	        |
| POPI 3,4 bisphosphate      |         |    POP2	|  POP2  | POPI2D in CHARMM36, new lipidome mapping available				    |
| POPI 3,4,5 trisphosphate   |         |    POP3	|  POP3  | POPI34 in CHARMM36, new lipidome mapping available				    |
| POPI 4 phosphate           |         |    POP4	|  POP4  | POPI14 in CHARMM36, new lipidome mapping available				    |
| POPI 5 phosphate           |         |    POP5	|  POP5  | POPI15 in CHARMM36, new lipidome mapping available				    |
| POPI 4,5 bisphosphate      |         |    POP6	|  POP6  | POPI24 in CHARMM36, new lipidome mapping available				    |
| POPI 3,5 bisphosphate      |         |    POP7	|  POP7  | POPI2A in CHARMM36, new lipidome mapping available				    |
| POPA                       |   D21   |    POPA    |  POPA  | New lipidome mappi, new lipidome mapping available                   |
| POPC                       |   POV   |    POPC    |  POPC  | New lipidome mappi, new lipidome mapping available                   |
| POPE                       |   6OU   |    POPE    |  POPE  | New lipidome mappi, new lipidome mapping available                   |
| POPE_SMILES                |   POES  |    POPE    |  POPE  | From SMILES string, new lipidome mapping available                   |
| POPG                       |   PGW   |    POPG    |  POPG  | New lipidome mapping available                                       |
| POPI   				     |         |    POPI    |  POPI  | New lipidome mapping available                                       |
| POPS                       |   D39   |    POPS    |  POPS  | New lipidome mapping available                                       |
| RAMP                       |         |    RAMP    |  RAMP  | Non-native in CHARMM36                                               |
| REMP                       |   KDL   |    REMP    |  REMP  |                                                                      |
| SSM1		                 |         |    SSM1    |  SSM1  | SSM in CHARMM36, new lipidome mapping available                      |
| TMM                        |         |    TMM1    |  TMM1  | Non-native in CHARMM36                                               |
| TMMA                       |         |    TMMA    |  TMMA  | Non-native in CHARMM36                                               |
| undecaprenyl phosphate     |   5TR   |    UNP1    |  UDP1  |                                                                      |
| undecaprenyl pyrophosphate |         |    UDP2    |  UDP2  | UNDPP in CHARMM36                                                    |


### Post-translational Modifications

| Name  | CCDName | CHARMMName | CGName  |  
| ----- | ------- | ---------- | ------- |
| CYSD	|         | CYSD       | CYSD    |
| CYSF	|         | CYSF       | CYSF    |
| CYSG	|         | CYSG       | CYSG    |
| CYSP	|         | CYSP       | CYSP    |
| CYST	|         | CYST       | CYST    |
| GLYM	|         | GLYM       | GLYM    |



### SMILES strings

Note any modification which does not affect the order of atoms can be made (e.g. changes to chiral centres/charges/double bonds). The name of the SMILES string must differ from the CHARMM name if wishing to utilise userCCD codes, and should not be more than 4 charaters if wishing to convert to CG. Ideally, it should also differ from utilised CCD codes.

| POPE_SMILES                |   POES  |  POPE  |```CCCCCCCC\C=C/CCCCCCCC(=O)O[C@H](COC(=O)CCCCCCCCCCCCCCC)COP(O)(=O)OCCN``` |


For issues and suggestions, please feel free to add these via GitHub or contact Kat Blow at katarina.blow[at]warwick.ac.uk