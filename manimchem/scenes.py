import hashlib
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Dict, Union

import numpy as np
import pandas as pd
from manimlib.animation.fading import FadeOut, FadeIn
from manimlib.imports import Scene, Circle, ShowCreation, COLOR_MAP, Line, Write, MovingCameraScene, VMobject, RIGHT, \
    VGroup, LEFT, LEFT_SIDE, RIGHT_SIDE, \
    GrowFromCenter, DOWN, IntegerMatrix, Arrow, ReplacementTransform, UP, ThreeDScene, Transform, ApplyMethod, Matrix
from manimlib.mobject.svg.tex_mobject import TextMobject
from rdkit.Chem import AllChem

from manimchem.rdkit_utils import to_rdkit_mol

# --- Data
#
# * ARTEMISININ
#     https://pubchem.ncbi.nlm.nih.gov/compound/68827
#     https://en.wikipedia.org/wiki/Artemisinin
#   We should also get its conformers from Pubchem3D.
#   Tu You You, the discovery of artemisinin and Nobel Prize in Physiology or Medicine:
#     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4966551/
#     https://www.nobelprize.org/prizes/medicine/2015/tu/facts/
#   We could also work with dihydroartemisinin.
#

DATA_DIR = Path(__file__).parent.parent / 'data'

ARTEMISININ_PUBCHEM_SMILES = 'CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C'
ARTEMISININ_PUBCHEM_ISOMERIC_SMILES = 'C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C'
ARTEMISININ_INCHI = 'InChI=1S/C15H22O5/c1-8-4-5-11-9(2)12(16)17-13-15(11)10(8)6-7-14(3,18-13)19-20-15/' \
                    'h8-11,13H,4-7H2,1-3H3/t8-,9-,10+,11+,13-,14-,15-/m1/s1'
ARTEMISININ_PUBCHEM_2D = DATA_DIR / 'artemisinin' / 'Structure2D_CID_68827.sdf'
ARTEMISININ_PUBCHEM_3D = DATA_DIR / 'artemisinin' / 'Conformer3D_CID_68827.sdf'
ARTEMISININ_PUBCHEM_JSON = DATA_DIR / 'artemisinin' / 'compound_CID_68827.json'

SMILES_ACTIVITIES = (
    ('CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C', True),
    ('FC(F)(F)c1ccc(nc1)N2CCOCC2', False),
    ('CCc1nnc(NC(=O)C(C)Sc2ccccc2)s1', False),
    ('CCCOC(=O)C1=C(C)NC2=C(C1c3ccc(cc3)[N+](=O)[O-])C(=O)CC(C2)c4ccccc4OC', True),
    ('OC(=O)CCc1nnc(SCC(=O)Nc2cc(ccc2Cl)C(F)(F)F)nc1O', False),
    ('Cc1nc2cc(NC(=O)Nc3cccc(Br)c3)ccc2n1C', False),
    ('COC(=O)c1ccc(C)c(NC(=O)CN(c2ccc(C)c(Cl)c2)S(=O)(=O)c3ccccc3)c1', False),
    ('CCOc1ccc(cc1)N2C(=O)C3=C(CCS3)N=C2SCC(=O)Nc4ccccc4', False),
    ('Cc1ccc(cc1)\\N=C\\c2ccc(O)cc2', False),
    ('Cc1ccc(cc1)C2SCCN2S(=O)(=O)C', False),
    ('CCc1ccc(cc1)C(=O)CC2=Cc3cccc(OC)c3OC2=O', False),
    ('COc1cccc2C=C(C(=O)Oc12)c3csc(Nc4ccccc4Br)n3', False),
    ('COc1ccc(CN2CCN(CC2)c3ncnc4c3cnn4c5ccccc5)cc1', False),
    ('CCCC(=O)OC1=COC(=CC1=O)CSc2nnc(NC(=O)C3CC3)s2', False),
    ('NC1=NC(c2ccc(O)cc2)n3c(N1)nc4ccccc34', False),
    ('Cc1c(cnn1c2ccccc2)C(=O)N3CCCN3C(=O)CCC=C', False),
    ('CCc1ccccc1NC(=O)c2cc(cn2C)S(=O)(=O)N3CCc4ccccc34', False),
    ('COc1ccc(cc1)C2N3C(=NC(=C2C(=O)OC(C)C)C)S\\C(=C/c4cc(OC)ccc4OC)\\C3=O', False),
    ('NC(=N)NC(=N)N(c1ccccc1)c2ccccc2', False),
    ('Cc1cccc(c1)C(=O)N2CCN=C2SCc3ccccc3C', False),
    ('CCCC1=NN2C(=N)\\C(=C\\c3ccc(O)c(OCC)c3)\\C(=O)N=C2S1', False),
    ('Cc1sc(NC(=O)C2CC2)c(C#N)c1C', False),
    ('CN(C)CCCNC(=O)c1ccc2[nH]c3CCCCc3c2c1', False),
    ('CCSc1nc(O)c2c[nH]nc2n1', False),
    ('OC(=O)CCc1nc2cccnc2n1Cc3ccccc3', False),
    ('CCc1oc(CCC(=O)Nc2cccc3ccccc23)cc1', False),
    ('COc1cc(CN2CCNC(=O)C2CC(=O)NCC=C)cc(OC)c1', False),
    ('CC1CC(C)(C)NC(=S)N1CC(=O)Nc2cccc(Cl)c2', False),
    ('CC(CO)N(Cc1ccccc1)C(=O)c2csc(Cc3ccc(C)cc3)n2', False),
    ('Cc1onc(NS(=O)(=O)c2ccccc2[N+](=O)[O-])c1', False),
    ('CCc1ccc(NC(=O)CSCC(=O)Nc2ccc(OC)cc2)cc1', False),
    ('CCCCCNC(=O)c1ccc(cc1)N(Cc2ccc(C)cc2)S(=O)(=O)C', False),
    ('O=C1S\\C(=C\\c2ccccc2)\\C(=O)N1CC#C', False),
    ('FC(F)(F)c1cccc(NC(=O)CC2Sc3ccccc3NC2=O)c1', False),
    ('COc1ccc(NC(=O)c2sc3nc4CCCCCc4cc3c2N)cc1', False),
    ('COc1ccc(cc1OC)C(=O)\\C=C\\Nc2ccc(C)c(C)c2', False),
    ('CC(=C)CSc1nnc2N(Cc3ccccc3)C(=O)c4c5CCCc5sc4n12', False),
    ('Fc1ccc(cc1)S(=O)(=O)N2CCCCC2CCNC(=O)C(=O)Nc3ccc(Cl)cc3', False),
    ('CCC1Sc2ccc(cc2NC1=O)S(=O)(=O)CCC(=O)Nc3ccc(cc3)C(=O)OC', False),
    ('CCN(CC)CCNC(=O)COCc1cc(on1)c2occc2', False),
    ('CN(Cc1cc2ccccc2[nH]1)Cc3nnc(SCc4ccccc4C)n3CC5CCCO5', False),
    ('CCCCOC(=O)c1cccc(c1)c2oc(\\C=C/3\\C(=NN(C3=O)c4ccc(cc4)S(=O)(=O)N)C)cc2', False),
    ('Cc1ncc(Cc2ccccc2)c(n1)C3CCCN(C3)C(=O)C4CC4', False),
    ('Fc1cccc(Cl)c1CSCCNC(=O)c2ccc(Cl)cc2', False),
    ('CCC1CCCCN1C(=O)C2=CN(Cc3cccc(OC)c3)C=C(C(=O)NCC(C)C)C2=O', False),
    ('CC(NC(=O)Nc1ccc(Cl)c(Cl)c1)c2ccccc2', False),
    ('O=C(NC1CC1)C2=CN(Cc3ccccc3)C=C(C(=O)N4CC[C@@H]5CCCC[C@H]5C4)C2=O', False),
    ('Cc1nc(Cl)nc(Nc2ccccc2)c1[N+](=O)[O-]', False),
    ('COc1ccc(cc1)c2ccc(CCC(=O)O)n2c3ccccc3OC', False),
    ('CCOC(=O)N1N=C(C)CC1(O)c2ccc(Br)cc2', False),
    ('FC(F)(F)c1nnc2sc(nn12)c3ccncc3', False),
    ('CCOc1ccccc1N2CCN(CC2)S(=O)(=O)c3cc(OC)c(Br)cc3OC', False),
    ('COc1cc(CNC2CCCC2)cc(Cl)c1OCC=C', False),
    ('Cn1ccnc1SCc2cc(O)nc3ccccc23', False),
    ('COC(=O)c1cc(NC(=O)C2COc3ccccc3O2)ccc1Cl', False),
    ('COCCCNC(=O)CN(Cc1cccc(Cl)c1)S(=O)(=O)c2ccc(C)cc2', False),
    ('CCOc1cc(cc(OCC)c1OCC)C(=O)Nc2sc3CCCCCc3c2C(=O)N', False),
    ('Cc1c2C=NN(CC(=O)Nc3cccc(F)c3)C(=O)c2c(C)n1Cc4ccc(F)cc4', False),
    ('Oc1ccccc1CNN2C(=O)c3ccccc3N=C2c4ccc(F)cc4', False),
    ('OC(=O)\\C(=C/c1oc(cc1)c2ccc(Cl)cc2)\\C#N', False),
    ('COC(=O)c1cc(NC(=O)c2oc(cc2)c3ccc(cc3)[N+](=O)[O-])cc(c1)C(=O)OC', False),
    ('CCc1ccc(CNC(=O)C(=O)Nc2cccc(Cl)c2Cl)cc1', False),
    ('CC(C)n1c(nc2N(C)C(=O)NC(=O)c12)N3CCC(C)CC3', False),
    ('[O-][N+](=O)c1cccc(CN2CCN(C\\C=C\\c3ccccc3)CC2)c1', False),
    ('COc1ccccc1NC(=O)C(=O)NCc2ccccc2C', False),
    ('COC(=O)c1c2CCCCc2sc1NC(=O)\\C=C\\c3cccc(c3)[N+](=O)[O-]', False),
    ('COCCOC(=O)c1ccc(NC(=O)c2ccc(cc2)C(C)(C)C)cc1', False),
    ('CC(C)C(NC(=O)c1cccc(OCC=C)c1)c2ccccc2', False),
    ('CCc1ccc(NC(=O)C(NS(=O)(=O)c2cccc3nsnc23)C(C)C)cc1', False),
    ('COc1ccc(OC)c(NC(=O)C2CCN(CC2)S(=O)(=O)c3ccc(C)cc3)c1', False),
    ('CCOc1ccccc1Nc2cc(C)nc3ncnn23', False),
    ('CC(N1N2C(=CC(=O)N=C2c3ccccc13)C)C(=O)NCc4ccc(C)cc4', False),
    ('CC(C)c1ccc(cc1)N(CC(=O)NCc2ccc(C)cc2)S(=O)(=O)c3c(C)n[nH]c3C', False),
    ('Fc1ccc(cc1)S(=O)(=O)N2CCC(CC2)C(=O)Nc3ccc(F)c(Cl)c3', False),
    ('COc1ccc(cc1OC)N2C(=O)C3=C(CCS3)N=C2SCc4cccc(c4)[N+](=O)[O-]', False),
    ('COc1ccccc1NC(=O)C2Sc3nnc(c4ccccc4)n3NC2c5ccccc5', False),
    ('OC(COc1ccc(Cl)cc1)CN2CCc3ccccc3C2', False),
    ('CN(CC(=O)Nc1cc(ccc1C)[N+](=O)[O-])S(=O)(=O)c2ccc(Cl)cc2', False),
    ('CN1CCN(CC2CCC=CC2)CC1', False),
    ('COc1ccc(cc1)S(=O)(=O)C2=CN(Cc3ccccc3Cl)c4cc(OC)c(OC)cc4C2=O', False),
    ('CCNC(=S)N\\N=C\\c1cccc(OCc2ccc(F)cc2)c1', False),
    ('CCOC(=O)C1=C(CN2CCN(CC2)c3cc(Cl)ccc3C)NC(=O)NC1c4ccc(F)cc4', False),
    ('CCCC(=O)NC(c1occc1)c2c(O)ccc3ccccc23', False),
    ('CC(=O)Nc1ccc(CN2CCC3(CCC(CNC(=O)C4CCCO4)O3)CC2)cc1', False),
    ('Clc1ccc(NC(=O)CN2C(=O)N(Cc3occc3)C(=O)c4sc5ccccc5c24)cc1', False),
    ('CCCOc1ccc(cc1)C2N(Cc3ccncc3)C(=O)C(=C2C(=O)c4ccc(C)c(F)c4)O', False),
    ('CC(=O)Nc1ccc(N\\C=C\\C(=O)c2occc2)cc1', False),
    ('CCOc1ccc(OCC(=O)Nc2c(oc3ccccc23)C(=O)Nc4cccc(F)c4)cc1', False),
    ('COc1ccccc1OC2=C(\\C=C(/C#N)\\S(=O)(=O)c3ccc(Cl)cc3)C(=O)N4C=CC=CC4=N2', False),
    ('Cc1ccc(N2C(=O)c3nccnc3C2=O)c(C)c1', False),
    ('CCn1c(COc2ccc(C)cc2)nnc1SCC(=O)Nc3ccc(cc3)[N+](=O)[O-]', False),
    ('CCOc1ccc(cc1OCC)C(=O)Nc2nnc(s2)c3ccccc3', False),
    ('CCc1ccc(cc1)C2NC(=O)NC(=C2C(=O)OCC3CCCO3)C', False),
    ('CCOc1ccc(cc1)N(C2CS(=O)(=O)C=C2)C(=O)c3cc(OC)c(OC)c(OC)c3', False),
    ('CCOC(=O)C(C)Oc1ccc(cc1)C(=O)Nc2scc(c3ccccc3)c2C(=O)OCC', False),
    ('Fc1cccc(NC(=O)c2cc3C(=O)N4C=CC=CC4=Nc3s2)c1', False),
    ('Cc1ccc(cc1)N2C(Nc3ccc(C)c(C)c3)c4ncccc4C2=O', False),
    ('CCCOc1ccc(\\C=C\\C(=O)Nc2cc(C)ccc2OC)cc1OC', False),
    ('Fc1ccccc1OCCOc2cccc3cccnc23', False),
    ('CCS(=O)(=O)c1nc(c(NCCCOC)s1)S(=O)(=O)c2ccc(Cl)cc2', False),
    ('Cc1ccc(Nc2nc(Nc3cc(Cl)cc(Cl)c3)nc(n2)N4CCOCC4)cc1', False),
    ('CCN(CC)c1ccc(NC(=O)c2cccc3ccccc23)c(C)c1', False),
    ('Cc1ccc(Cl)cc1NC(=O)C2CCN(CC2)C(=O)c3ccc(F)cc3', False),
    ('Ic1ccc(cc1)c2cc3ccccn3c2', False),
    ('Clc1ccc(NC(=O)c2ccc(cc2)S(=O)(=O)N3CCc4ccccc4C3)cc1', False),
    ('NC(=O)c1ccccc1NC(=O)Nc2ccccc2', False),
    ('CC(C)NC(=O)CSc1oc(CNC(=O)c2ccc(F)cc2)nn1', False),
    ('Cc1ccc(cc1)C(=O)\\C=C\\Nc2cc(Cl)ccc2O', False),
    ('CSCCNCc1c(nc2ccc(Cl)cn12)C(=O)N3CCCCCCC3', False),
    ('CCOC(=O)c1[nH]c2ccccc2c1NC(=O)c3ccccc3[N+](=O)[O-]', False),
    ('CCOC(=O)COc1ccc(Br)cc1\\C=C\\c2nc(O)nc(O)c2[N+](=O)[O-]', False),
    ('CCN1CCN(CC1)C(=O)C23CC4CC(C)(CC(C)(C4)C2)C3', False),
    ('CSc1cccc(NC(=O)CN(Cc2ccc(Cl)cc2)S(=O)(=O)c3ccc(Cl)cc3)c1', False),
    ('CCCc1cc(O)nc2nnc(SCC(=O)Nc3ccc(F)cc3F)n12', False),
    ('COc1cc(cc(OC)c1OC)C(=O)N2CCN=C2SCc3cccc(Cl)c3', False),
    ('Cc1ccc2nc(c(Nc3c(C)cccc3C)n2c1)c4ccccc4F', False),
    ('CC(C)Oc1ccc(cc1)C(=O)Nc2ccc(O)c(c2)c3nc4ccccc4s3', False),
    ('CC(=O)Oc1cccc(c1)C(=O)Nc2cnc(O)nc2O', False),
    ('CCN1C(=O)C(O)(CC(=O)c2ccccc2)c3cc(Br)ccc13', False),
    ('Cc1ccc(cc1I)\\N=C\\c2ccc(cc2)C#N', False),
    ('CCCc1cc2C(=O)C(=C(C)Oc2cc1OCC(=O)O)c3ccc4OCOc4c3', False),
    ('NS(=O)(=O)c1ccc(CCNC(=O)c2cccc(Br)c2)cc1', False),
    ('Cc1cccc(NC(=O)CSc2ccc3nnc(c4cccnc4)n3n2)c1', False),
    ('Cc1oc(nc1CSCC(=O)NCc2cccnc2)c3cccs3', False),
    ('OC(=O)CCC(=O)N(Cc1ccccc1)Cc2ccccc2', False),
    ('CC(C)C(=O)Nc1cccc(c1)c2cn3c(C)cccc3n2', False),
    ('COCCCNc1nc2N(C)C(=O)NC(=O)c2n1CC(O)COc3ccccc3C', False),
    ('CCOc1ccccc1C2CC(Nc3nc(N)nn23)c4ccc(Cl)cc4', False),
    ('COc1ccc(cc1OC)C(=O)OC2=COC(=CC2=O)CSc3nccc(C)n3', False),
    ('COc1ccc(Cn2c(nc3ccccc23)c4nonc4NC(=O)C)cc1', False),
    ('CCC(=O)Nc1ccc(cc1)S(=O)(=O)Nc2nccs2', False),
    ('COc1ccc(OC)c(CN2C(SCC2=O)c3ccccc3OC)c1', False),
    ('C1CCC(C1)Sc2nnc(c3ccccc3)c(n2)c4ccccc4', False),
    ('Cc1[nH]c2ccc(cc2c1C)C(=O)N(Cc3ccc(OCC4CCCO4)cc3)Cc5cccnc5', False),
    ('CCCCN1C(=O)NC(=O)C(=C2CC(Sc3ccccc3N2)c4ccc(Cl)cc4)C1=O', False),
    ('COc1ccc(cc1OC)C(=O)\\N=C\\2/Sc3cc(ccc3N2CC#C)S(=O)(=O)N', False),
    ('CSc1ccccc1c2oc(NC(=O)c3ccc(Cl)s3)nn2', False),
    ('Cc1ccc(\\C=N\\NC(=O)C(=O)Nc2cccc(C)c2)cc1', False),
    ('OCCN1C(=Nc2ccccc2C1=O)\\C=C\\c3cccc(Br)c3', False),
    ('COC(=O)C1=C(CCS1)NC(=O)c2ccccc2F', False),
    ('Clc1cccc(NC(=O)C2=CC(=O)c3ccccc3O2)c1', False),
    ('COc1cc(OC)c(OC)cc1CNc2cccc(C)c2', False),
    ('Cc1ncc2CN(CCc2c1CNC(=O)c3ccsc3)C(=O)C(=O)c4ccccc4', False),
    ('COc1ccc(NC(=O)c2nnn(c2C)c3ccc(C)c(F)c3)cc1Cl', False),
    ('COC(=O)CN1C=Nc2c(nc3CCCCCn23)C1=O', False),
    ('Cc1ccc(cc1)N2CCN(CC2)\\N=C/c3cccc(Br)c3', False),
    ('COc1ccc(cc1C(=O)N)S(=O)(=O)N2CCN(Cc3ccccc3)CC2', False),
    ('CSc1nnc(NC(=O)c2nc(ncc2Cl)S(=O)(=O)Cc3ccc(C)cc3)s1', False),
    ('COc1ccc(NC(=O)CC2(CC(=O)O)CCCC2)cc1', False),
    ('COc1cc(cc2CN(Cc3ccc[nH]3)CCOc12)c4ccc(SC)cc4', False),
    ('CC1CCc2c(C1)sc3NC(NC(=O)c23)c4cccc(Br)c4', False),
    ('CC1CCCN(C1)S(=O)(=O)c2ccc3nc(O)cc(C(=O)NCc4ccccc4Cl)c3c2', False),
    ('Fc1ccc2SC=C(N3CCN(CC3)c4cccc(c4)C(F)(F)F)C(=O)c2c1', False),
    ('CCOc1ccc(cc1Cl)C(=O)Nc2ccccc2N3CCN(CC3)S(=O)(=O)C', False),
    ('COc1ccc(OC)c(c1)S(=O)(=O)N2CCC(CC2)C(=O)Nc3ccccc3F', False),
    ('COc1ccc(NC(=O)Cn2cc(C(=O)C3CCCCC3)c4ccccc24)cc1', False),
    ('COc1ccc(CN2CCN(Cc3ccccc3F)CC2)c(OC)c1', False),
    ('COc1cc(ccc1Cl)S(=O)(=O)NCCCn2ccnc2', False),
    ('CCc1ccc(cc1)C(Nc2ccccn2)c3ccc4ccc(C)nc4c3O', False),
    ('CN1C(=C(C(=O)N(C)C1=O)S(=O)(=O)Nc2ccc(OC(F)(F)F)cc2)C', False),
    ('COC[C@@H]1CCCN1S(=O)(=O)c2cc(OC)ccc2OC', False),
    ('COc1ccc(cc1S(=O)(=O)N2CCCC2)C(=O)N3CCN(CC3)c4ccccc4', False),
    ('CC1=C(C)CC2C(C1)C(=O)C3=C(C4c5ccccc5C3c6ccccc46)C2=O', True),
    ('CN1\\C(=N\\C(=O)c2cccc(c2)[N+](=O)[O-])\\Sc3c(C)cc(C)cc13', False),
    ('Fc1ccc(cc1)N2CCN(CC2)C(=O)Cn3ncc4COc5ccccc5c34', False),
    ('CC(Oc1cccc(C)c1)c2nc3ccccc3n2CC=C', False),
    ('COc1cc(NC(=O)c2ccccc2)c(OC)cc1Cl', False),
    ('COc1ccc(\\C=N/NC(=O)c2ccc(C)nc2)cc1', False),
    ('COc1ccc(cc1OC)N(CC(=O)NC2CCCC2)S(=O)(=O)C', False),
    ('COC(=O)c1c(C)[nH]c2c3ccccc3c(O)c2c1c4cc(Br)c(OCC#C)c(OC)c4', False),
    ('Cc1cc(C)c(c(C)c1)S(=O)(=O)Nc2ccc(Cl)cc2Cl', False),
    ('COc1cc(CN2CCN(CC2)C(=O)c3cccc(F)c3)cc(OC)c1OC', False),
    ('CCCCNC(=O)CN1C(=O)\\C(=C/c2ccc(F)cc2)\\Oc3ccccc13', False),
    ('ON1N(C(=O)Nc2ccccc12)c3cccc(c3)[N+](=O)[O-]', False),
    ('Clc1ccc(cc1)S(=O)(=O)NCCSCc2occc2', False),
    ('Cc1cc(Cl)ccc1OCC(=O)Nc2sc(Cc3ccccc3)c(C)c2C(=O)N', False),
    ('[O-][N+](=O)c1ccc(cc1)S(=O)(=O)Nc2cc3CCCN4C(=O)CCc(c2)c34', False),
    ('COc1ccc(cc1[N+](=O)[O-])C(=O)Nc2c(C)cc(C)cc2C', False),
    ('COC(=O)c1c2CCC(C)Cc2sc1NC(=O)COc3cc(Cl)ccc3Cl', False),
    ('CSCCN1C(SCC1=O)c2cccnc2', False),
    ('COc1ccccc1CNCc2c(C)n(Cc3ccccc3C)c(C)c2C(=O)O', False),
    ('COc1ccccc1N2C(=O)NC(=O)C(CCc3ccncc3)(CCc4ccncc4)C2=O', False),
    ('COc1ccc(NC(=O)c2cnc3c(c(C)nn3c4ccccc4)c2Cl)c(OC)c1', False),
    ('Cc1ccc2nc(sc2c1)c3ccc(NC(=O)CN4CCN(CC4)C(=O)C5CCCO5)cc3', False),
    ('CCCN(c1ccc(cc1)C(C)C)S(=O)(=O)c2c(C)nc(O)nc2O', False),
    ('CN(C)S(=O)(=O)c1ccc(cc1)C(=O)C2=C(O)C(=O)N(Cc3cccnc3)C2c4cccnc4', False),
    ('CCCCC(=O)Nc1nc2ccc3nc(C)sc3c2s1', False),
    ('Cc1cc(C)nc(NS(=O)(=O)c2ccc(N\\C=C\\C(=O)c3ccc(Br)cc3)cc2)n1', False),
    ('CCN(CC)C(=O)CSc1nc2c3ccccc3nc2c(O)n1CC=C', False),
    ('CCCOc1ccc(\\C=C\\C(=O)NCCOC)cc1', False),
    ('COc1ccc(NC(=O)CSc2cn(CCNC(=O)c3ccc(OC)cc3)c4ccccc24)cc1', False),
    ('Cc1ccccc1NC(=O)c2c(C)c(Cl)c(C)nc2Cl', False),
    ('COc1ccc(Cl)cc1NC(=O)CCCN2C(=O)c3ccccc3C2=O', False),
    ('CCN1C(=NN/C/1=C\\2/C=CC=CC2=O)SCC(=O)Nc3nc(cs3)c4ccccc4', False),
    ('Brc1ccc(NC(=O)c2ccc(NS(=O)(=O)c3ccccc3)cc2)cc1', False),
    ('COP(=O)(OC)c1nc(oc1NCCc2ccccc2)c3ccc(cc3)[N+](=O)[O-]', False),
    ('CCOC(=O)C1=NN(C(=O)c2c(NC(=O)Cc3ccc(F)cc3)scc12)c4ccc(OC)cc4', False),
    ('OC1(CC(=O)c2cc(Cl)ccc2O1)C(F)(F)C(F)F', False),
    ('COc1ccccc1NC(=O)CN(c2ccccc2)S(=O)(=O)C', False),
    ('Oc1ccc(\\C=N\\NC(=O)C2COc3ccccc3O2)cc1', False),
    ('CC(=O)Nc1cccc2C(=O)N(C(=O)c12)c3ccccc3OC(=O)C', False),
    ('Brc1ccc(OCC(=O)N2CCC(Cc3ccccc3)CC2)cc1', False),
    ('Fc1ccc(cc1)c2[nH]c(nc2SCC(=O)Nc3ccc4OCCOc4c3)c5ccccc5', False),
    ('COc1cc(NC(=O)COc2ccccc2[N+](=O)[O-])ccc1NC(=O)c3cccs3', False),
    ('Clc1ccc(cc1)C(=O)Nc2nc(ns2)c3ccccc3', False),
    ('CN(CC(=O)NC1CC1)S(=O)(=O)c2ccc(Cl)cc2', False),
    ('CCOC(=O)C1(CC2CCCCO2)CCN(CC1)C(=O)Nc3ccccc3SC', False),
    ('COc1cccc(NC(=O)C2=Cc3sc4CCCCc4c3CS2)c1', False),
    ('CCCn1c(nc2ccccc12)C(C)NC(=O)c3ccc(F)cc3', False),
    ('Cc1ccc(cc1)S(=O)(=O)N(CC(=O)Nc2ccc(C)c(C)c2)Cc3ccccc3', False),
    ('Cc1cc(C)c(NC(=O)COc2cccc(Cl)c2)c(C)c1', False),
    ('COc1ccc(OC)c(c1)N(CC(=O)Nc2cccc(c2)C(F)(F)F)S(=O)(=O)C', False),
    ('CCN1C(=O)C=C(OCC(=O)Nc2cccc(Cl)c2)c3ccccc13', False),
    ('COc1ccc(Br)cc1C=C2C(=O)NC(=O)NC2=O', False),
    ('COc1cc(\\C=C\\2/C(=O)N=C3C=C(C)ON3C2=N)cc(CC=C)c1O', False),
    ('CCOc1ccc(cc1)n2c(C)cc(\\C=C/3\\NC(=O)N(Cc4ccccc4F)C3=O)c2C', False),
    ('CCc1cccc(NC(=O)NCc2ccc(F)c(F)c2)c1', False),
    ('Fc1ccc(NC2=C(Cl)C(=O)N(Cc3ccc(Cl)cc3)C2=O)cc1Cl', False),
    ('CCCN(CC1CC1)S(=O)(=O)c2sc3CN(CCc3c2C(=O)OC)C(=O)CCn4cncn4', False),
    ('Cc1ccccc1CNC(=O)CN2C(=O)c3cccnc3Oc4ccc(Cl)cc24', False),
    ('Fc1ccc(cc1)N2CCN(CC2)C(=O)c3ccc4NC(CS(=O)(=O)Cc5ccccc5Cl)C(=O)Nc4c3', False),
    ('COc1ccc(Cl)cc1S(=O)(=O)N2CCCC(C2)C(=O)O', False),
    ('COc1ccc(Cl)cc1NC(=O)c2c(C)onc2c3c(Cl)cccc3Cl', False),
    ('COc1ccc(NC(=O)CN(C)S(=O)(=O)c2ccc3N(C)C(=O)N(C)C(=O)c3c2)cc1Cl', False),
    ('COc1ccc(cc1OC)c2noc(n2)c3ccc(Cl)c(Cl)c3', False),
    ('CCC(C(=O)Nc1nc2ccc(C)cc2s1)c3ccccc3', False),
    ('CCCCCCCCCn1c(nc2N(C)C(=O)NC(=O)c12)N3CCN(CC3)c4ccccc4', False),
    ('COc1ccccc1CNC(=O)C(=O)NCCC2CCCCN2S(=O)(=O)c3ccc(F)cc3', False),
    ('Clc1cc(NC(=O)COc2ccc(cc2)C3CCCCC3)ccc1N4CCOCC4', False),
    ('OC(=O)c1ccc(NC(=O)CSc2nnnn2c3cccc4ccccc34)cc1', False),
    ('CC(C)(C)c1ccc(cc1)C2=NNC(=S)N2C3CCCCC3', False),
    ('COc1ccc(cc1OC)N2C=CN(Cc3cccc(F)c3)C(=O)C2=O', False),
    ('CC(C)CCNC(=O)C(=O)Nc1c2CS(=O)(=O)Cc2nn1c3cccc(Cl)c3', False),
    ('OC(=O)\\C=C\\C(=O)Nc1ccc(cc1)c2ccc(NC(=O)\\C=C\\C(=O)O)cc2', False),
    ('Cc1ccc(N\\C=C(/C#N)\\c2nc(cs2)c3ccc(cc3)[N+](=O)[O-])cc1C', False),
    ('CC(=O)NCCc1ccc(cc1)S(=O)(=O)N2CCCC2', False),
    ('CC1=C(C(C2=C(CC(C)(C)CC2=O)N1)c3cccc(Cl)c3)C(=O)OCC=C', False),
    ('Cc1ccc(CNC(=O)CCCN2C(=O)N(Cc3cccc(F)c3)c4ccccc4C2=O)cc1', False),
    ('CCN1C(=O)C(C2CC1(C)Oc3ccccc23)C(=O)Nc4ccc(Br)cc4', False),
    ('CCN1C(=O)CC(Sc2ccccc2N)C1=O', False),
    ('COC(=O)c1cccc(c1)n2c(C)cc(C=C3C(=O)NC(=O)NC3=O)c2C', False),
    ('CCOc1cccc(OCCCOc2ccc(Cl)c(C)c2)c1', False),
    ('Clc1cccc(c1)C(=O)N\\N=C\\c2ccc(Cl)c(Cl)c2', False),
    ('COc1ccc(NC(=O)c2ccc(cc2)N(C)S(=O)(=O)c3ccc(OC)c(OC)c3)c(OC)c1', False),
    ('CC(C)c1ccc(C)cc1OCC(=O)N\\N=C\\c2cc(I)cc(I)c2O', False),
    ('COC(=O)c1cc(c(Cl)cc1Cl)S(=O)(=O)n2nc(C)cc2C', False),
    ('Cc1ccccc1OCC(=O)N2CCC3(CC2)CCN(C3)C(=O)CC(=O)OC(C)(C)C', False),
    ('Cc1cc(ccc1NC(=O)CSc2nnc(c3ccccc3)n2C)[N+](=O)[O-]', False),
    ('COc1cc(OC)c(cc1Br)C2C(=C(N)OC3=C2C(=O)CC(C)(C)C3)C#N', False),
    ('COc1ccc(CCNC(=O)CC2=C(C)c3c(O)cc(C)cc3OC2=O)cc1', False),
    ('[O-][N+](=O)c1ccccc1C(=O)N2CCN(CC2)c3ccc(F)cc3', False),
    ('CCCCN(C)C(=O)NC1CN(C(=O)C1)c2ccc3OCCOc3c2', False),
    ('COc1ccc(C)cc1NC(=O)c2ccc(OCC(=O)Nc3cc(Cl)ccc3Cl)cc2', False),
    ('O=C(Nc1ccc(cc1)S(=O)(=O)N2CCCCC2c3cccnc3)c4ccc(cc4)S(=O)(=O)N5CCCCC5', False),
    ('COc1ccc(NC(=S)N2CCN(CC2)c3ccccc3SC)c(C)c1', False),
    ('CC(C)OCCOC(=O)c1[nH]c2CC(CC(=O)c2c1C)c3ccccc3Cl', False),
    ('CCc1ccc(NC(=O)C(C)SCc2ccccc2)cc1', False),
    ('Cc1ccccc1N2CC(CC2=O)c3nc4ccccc4n3Cc5ccc(F)cc5', False),
    ('Cc1cccc(NC(=O)CS(=O)(=O)c2cn(CC(=O)N3CCCCC3)c4ccccc24)c1', False),
    ('CCN(CC)c1ccc(\\C=C\\2/N=C(OC2=O)c3cc(OC)c(OC)c(OC)c3)cc1', False),
    ('CC(N1C=Nc2sc3CCCCc3c2C1=O)C(=O)N\\N=C/c4ccc(O)cc4', False),
    ('Nc1ccc(NC(=O)c2cccc(c2)C(=O)Nc3ccc(N)cc3)cc1', False),
    ('Cc1nn(C(=O)CCCC(=O)NCc2ccccc2Cl)c3ccccc13', False),
    ('CCC(=O)N1CCN(CC1)c2ccc(NC(=O)c3cccc(C)c3)cc2', False),
    ('CN1C(=O)N(C(=O)c2c1c3cc(C)ccc3n2C)c4cccc(Cl)c4', False),
    ('COc1ccc(CNCc2ccccc2C)cc1', False),
    ('CCOC(=O)c1cnc2c(OC)cccc2c1Nc3ccc(F)cc3Cl', False),
    ('COCCCN(Cc1occc1)C(=O)c2cc3c(O)nc4ccc(C)cc4c3s2', False),
    ('Brc1ccc2N(CC#C)\\C(=N\\C(=O)C3CC3)\\Sc2c1', False),
    ('OC(=O)c1cccc(c1)S(=O)(=O)Nc2ccc(I)cc2', False),
    ('CN(c1ccc(Cl)cc1)S(=O)(=O)c2ccc3ncc(C(=O)N4CCCCCC4)c(O)c3c2', False),
    ('O=C1CC(NC2CCCCC2)C(=O)N1c3ccccc3', False),
    ('COC(=O)CN1C(=O)C2(C(=C(N)Oc3[nH]nc(C)c23)C#N)c4ccccc14', False),
    ('Fc1ccc(cc1)c2oc(\\C=C(/C#N)\\c3ccccc3)cc2', False),
    ('CCOCCCN1CN(c2nc3ccccc3nc12)S(=O)(=O)c4ccc(Cl)cc4', False),
    ('CCOC(=O)c1ccsc1NC(=O)c2ccccc2', False),
    ('O=C(\\C=C\\c1ccc2OCOc2c1)\\C=C\\c3ccc4OCOc4c3', False),
    ('COc1cccc(CN2CC(CCC2=O)C(=O)NCC(O)(CC=C)CC=C)c1', False),
    ('COc1ccc(CN2CCC(CC2)n3nccc3NC(=O)CC(C)C)c(F)c1', False),
    ('Cc1ccc(cc1)C(N2CCN(CC2)c3ccc(F)cc3)c4sc5ncnn5c4O', False),
    ('Cc1c(CCNS(=O)(=O)c2ccc(F)cc2)sc3nc(nn13)c4ccc(F)cc4', False),
    ('CCOCCN1CCN(CC1)S(=O)(=O)c2ccccc2', False),
    ('CC(=O)N1C(=C(Sc2nnc(c3ccccc3)n12)C(=O)C)c4ccc(OC(F)F)cc4', False),
    ('Cc1ccc(cc1)n2c(C)cc(C(=O)CSc3nnc4ccccn34)c2C', False),
    ('CCCCOc1ccc(NC(=O)C2CCN(CC2)c3cnccn3)cc1', False),
    ('Cc1cc(nn1CC(=O)Nc2cc(F)ccc2F)[N+](=O)[O-]', False),
    ('COC(=O)CN1\\C(=N\\C(=O)c2ccc3OCCOc3c2)\\Sc4cc(ccc14)C(=O)OC', False),
    ('CCS(=O)(=O)c1ccc(nn1)c2cccc(NC(=O)c3ccc(F)c(F)c3)c2', False),
    ('CC1CC(C)CN(C1)C(=O)C2CCN(CC2)S(=O)(=O)c3cccnc3', False),
    ('CC(=O)c1cccc(NC(=O)C(=O)c2cn(CC(=O)N3CCCCC3)c4ccccc24)c1', False),
    ('CCOC(=O)CN1CCN(CC1)S(=O)(=O)c2cccc(c2)C#N', False),
    ('COc1ccc(Cl)cc1NC(=O)c2ccc(CN3CCOCC3)cc2', False),
    ('CN(Cc1ccc(F)c(Br)c1)C2CCN(C)CC2', False),
    ('Oc1nc(nc2ccccc12)\\C(=C\\c3oc(cc3)c4cccc(Cl)c4)\\C#N', False),
    ('CCCCCCOc1ccc(cc1)C(=O)N\\N=C\\c2ccc(OC)c(OC)c2', False),
    ('OC1=C(C(N(Cc2cccnc2)C1=O)c3ccc(O)cc3)C(=O)c4ccc(Cl)cc4', False),
    ('CS(=O)(=O)N(CC(=O)NCc1ccc(Cl)cc1)c2cccc(F)c2', False),
    ('Brc1ccc(Nc2nnc(s2)C3=Cc4c(OC3=O)ccc5ccccc45)cc1', False),
    ('Fc1cccc(c1)C(=O)Nc2ccc3N(CCCc3c2)C(=O)C4CC4', False),
    ('CCNc1oc(\\C=C\\c2cccc(OC)c2)nc1C#N', False),
    ('[O-][N+](=O)c1cccc(NC(=O)CSc2ccc(nn2)c3cccnc3)c1', False),
    ('CCOC(=O)Cn1c(C)c(C)n2c3C(=O)N(Cc4cc(C)ccc4C)C(=O)N(C)c3nc12', False),
    ('CCOC(=O)C1=C(CN(C)Cc2ccccc2)NC(=O)NC1c3ccc(C)cc3', False),
)

# --- Aesthetics

NO_COLLISION_COLOR = COLOR_MAP['BLUE_E']
COLLISION_COLOR = COLOR_MAP['RED_E']
SELECTED_SUBSTRUCTURE_COLOR = COLOR_MAP['YELLOW_E']
ACTIVE_COLOR = COLOR_MAP['GREEN_E']
INACTIVE_COLOR = COLOR_MAP['RED_E']


def rgb2hex(r, g, b, percentage=True):
    if percentage:
        r, g, b = int(round(r * 255)), int(round(g * 255)), int(round(b * 255))
    r, g, b = max(0, min(r, 255)), max(0, min(g, 255)), max(0, min(b, 255))
    return "#{0:02x}{1:02x}{2:02x}".format(r, g, b)


#
# Ripped from MIT Licensed
#   https://github.com/patrickfuller/blender-chemicals/blob/master/blender_chemicals/atoms.json
#
# Copyright (C) 2012 Patrick Fuller, patrick-fuller.com
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
ATOM_INFO = {
    'Ac': {'color': [0.439216, 0.670588, 0.980392], 'radius': 1.114285},
    'Ag': {'color': [0.752941, 0.752941, 0.752941], 'radius': 0.914285},
    'Al': {'color': [0.74902, 0.65098, 0.65098], 'radius': 0.714285},
    'Am': {'color': [0.329412, 0.360784, 0.94902], 'radius': 1.0},
    'Ar': {'color': [0.501961, 0.819608, 0.890196], 'radius': 0.4057145},
    'As': {'color': [0.741176, 0.501961, 0.890196], 'radius': 0.657145},
    'Au': {'color': [1, 0.819608, 0.137255], 'radius': 0.77143},
    'B': {'color': [1, 0.709804, 0.709804], 'radius': 0.4857145},
    'Ba': {'color': [0, 0.788235, 0], 'radius': 1.22857},
    'Be': {'color': [0.760784, 1, 0], 'radius': 0.6},
    'Bi': {'color': [0.619608, 0.309804, 0.709804], 'radius': 0.914285},
    'Br': {'color': [0.65098, 0.160784, 0.160784], 'radius': 0.657145},
    'C': {'color': [0.564706, 0.564706, 0.564706], 'radius': 0.4},
    'Ca': {'color': [0.239216, 1, 0], 'radius': 1.02857},
    'Cd': {'color': [1, 0.85098, 0.560784], 'radius': 0.885715},
    'Ce': {'color': [1, 1, 0.780392], 'radius': 1.057145},
    'Cl': {'color': [0.121569, 0.941176, 0.121569], 'radius': 0.57143},
    'Co': {'color': [0.941176, 0.564706, 0.627451], 'radius': 0.77143},
    'Cr': {'color': [0.541176, 0.6, 0.780392], 'radius': 0.8},
    'Cs': {'color': [0.341176, 0.0901961, 0.560784], 'radius': 1.485715},
    'Cu': {'color': [0.784314, 0.501961, 0.2], 'radius': 0.77143},
    'Dy': {'color': [0.121569, 1, 0.780392], 'radius': 1.0},
    'Er': {'color': [0, 0.901961, 0.458824], 'radius': 1.0},
    'Eu': {'color': [0.380392, 1, 0.780392], 'radius': 1.057145},
    'F': {'color': [0.564706, 0.878431, 0.313725], 'radius': 0.2857145},
    'Fe': {'color': [0.878431, 0.4, 0.2], 'radius': 0.8},
    'Ga': {'color': [0.760784, 0.560784, 0.560784], 'radius': 0.742855},
    'Gd': {'color': [0.270588, 1, 0.780392], 'radius': 1.02857},
    'Ge': {'color': [0.4, 0.560784, 0.560784], 'radius': 0.714285},
    'H': {'color': [1, 1, 1], 'radius': 0.142857},
    'Hf': {'color': [0.301961, 0.760784, 1], 'radius': 0.885715},
    'Hg': {'color': [0.721569, 0.721569, 0.815686], 'radius': 0.857145},
    'Ho': {'color': [0, 1, 0.611765], 'radius': 1.0},
    'I': {'color': [0.580392, 0, 0.580392], 'radius': 0.8},
    'In': {'color': [0.65098, 0.458824, 0.45098], 'radius': 0.885715},
    'Ir': {'color': [0.0901961, 0.329412, 0.529412], 'radius': 0.77143},
    'K': {'color': [0.560784, 0.25098, 0.831373], 'radius': 1.257145},
    'La': {'color': [0.439216, 0.831373, 1], 'radius': 1.114285},
    'Li': {'color': [0.8, 0.501961, 1], 'radius': 0.82857},
    'Lu': {'color': [0, 0.670588, 0.141176], 'radius': 1.0},
    'Mg': {'color': [0.541176, 1, 0], 'radius': 0.857145},
    'Mn': {'color': [0.611765, 0.478431, 0.780392], 'radius': 0.8},
    'Mo': {'color': [0.329412, 0.709804, 0.709804], 'radius': 0.82857},
    'N': {'color': [0.188235, 0.313725, 0.972549], 'radius': 0.3714285},
    'Na': {'color': [0.670588, 0.360784, 0.94902], 'radius': 1.02857},
    'Nb': {'color': [0.45098, 0.760784, 0.788235], 'radius': 0.82857},
    'Nd': {'color': [0.780392, 1, 0.780392], 'radius': 1.057145},
    'Ni': {'color': [0.313725, 0.815686, 0.313725], 'radius': 0.77143},
    'Np': {'color': [0, 0.501961, 1], 'radius': 1.0},
    'O': {'color': [1, 0.0509804, 0.0509804], 'radius': 0.342857},
    'Os': {'color': [0.14902, 0.4, 0.588235], 'radius': 0.742855},
    'P': {'color': [1, 0.501961, 0], 'radius': 0.57143},
    'Pa': {'color': [0, 0.631373, 1], 'radius': 1.02857},
    'Pb': {'color': [0.341176, 0.34902, 0.380392], 'radius': 1.02857},
    'Pd': {'color': [0, 0.411765, 0.521569], 'radius': 0.8},
    'Pm': {'color': [0.639216, 1, 0.780392], 'radius': 1.057145},
    'Po': {'color': [0.670588, 0.360784, 0], 'radius': 1.085715},
    'Pr': {'color': [0.85098, 1, 0.780392], 'radius': 1.057145},
    'Pt': {'color': [0.815686, 0.815686, 0.878431], 'radius': 0.77143},
    'Pu': {'color': [0, 0.419608, 1], 'radius': 1.0},
    'Ra': {'color': [0, 0.490196, 0], 'radius': 1.22857},
    'Rb': {'color': [0.439216, 0.180392, 0.690196], 'radius': 1.342855},
    'Re': {'color': [0.14902, 0.490196, 0.670588], 'radius': 0.77143},
    'Rh': {'color': [0.0392157, 0.490196, 0.54902], 'radius': 0.77143},
    'Ru': {'color': [0.141176, 0.560784, 0.560784], 'radius': 0.742855},
    'S': {'color': [1, 1, 0.188235], 'radius': 0.57143},
    'Sb': {'color': [0.619608, 0.388235, 0.709804], 'radius': 0.82857},
    'Sc': {'color': [0.901961, 0.901961, 0.901961], 'radius': 0.914285},
    'Se': {'color': [1, 0.631373, 0], 'radius': 0.657145},
    'Si': {'color': [0.941176, 0.784314, 0.627451], 'radius': 0.62857},
    'Sm': {'color': [0.560784, 1, 0.780392], 'radius': 1.057145},
    'Sn': {'color': [0.4, 0.501961, 0.501961], 'radius': 0.82857},
    'Sr': {'color': [0, 1, 0], 'radius': 1.142855},
    'Ta': {'color': [0.301961, 0.65098, 1], 'radius': 0.82857},
    'Tb': {'color': [0.188235, 1, 0.780392], 'radius': 1.0},
    'Tc': {'color': [0.231373, 0.619608, 0.619608], 'radius': 0.77143},
    'Te': {'color': [0.831373, 0.478431, 0], 'radius': 0.8},
    'Th': {'color': [0, 0.729412, 1], 'radius': 1.02857},
    'Ti': {'color': [0.74902, 0.760784, 0.780392], 'radius': 0.8},
    'Tl': {'color': [0.65098, 0.329412, 0.301961], 'radius': 1.085715},
    'Tm': {'color': [0, 0.831373, 0.321569], 'radius': 1.0},
    'U': {'color': [0, 0.560784, 1], 'radius': 1.0},
    'V': {'color': [0.65098, 0.65098, 0.670588], 'radius': 0.77143},
    'W': {'color': [0.129412, 0.580392, 0.839216], 'radius': 0.77143},
    'Y': {'color': [0.580392, 1, 1], 'radius': 1.02857},
    'Yb': {'color': [0, 0.74902, 0.219608], 'radius': 1.0},
    'Zn': {'color': [0.490196, 0.501961, 0.690196], 'radius': 0.77143},
    'Zr': {'color': [0.580392, 0.878431, 0.878431], 'radius': 0.885715},
    'undefined': {'color': [0, 0, 0], 'radius': 0.405},
    'bond': {'color': [0.05, 0.05, 0.05], 'radius': 0.103}
}
for d in ATOM_INFO.values():
    # noinspection PyTypeChecker
    d['color'] = rgb2hex(*d['color'])


# --- Helpers

def stable_hash(string, hasher=hashlib.md5, fold_to=None):
    shash = int(hasher(string.encode('utf-8')).hexdigest(), 16)
    return shash if not fold_to else shash % fold_to


def random_tdt_train_data_selection(p=0.001, seed=0):
    df = (pd.
          read_csv(DATA_DIR / 'tdt' / 'training' / 'malariahts_trainingset.txt', sep='\t').
          query('Pf3D7_ps_hit != "ambiguous"')
          [['Canonical_Smiles', 'Pf3D7_ps_hit']].
          sample(frac=p, replace=False, random_state=seed))
    print('(')
    for smiles, active in zip(df['Canonical_Smiles'], df['Pf3D7_ps_hit']):
        active = 'True' if active == 'true' else 'False'
        smiles = smiles.replace('\\', '\\\\')
        print(f"    ('{smiles}', {active}),")
    print(')')
# random_tdt_train_data_selection()
# exit(22)


# --- Mobjects

class Molecule(VMobject):

    #
    # Bring inspiration:
    #   https://patrickfuller.github.io/molecules-from-smiles-molfiles-in-blender/
    #   https://github.com/Helpsypoo/primer
    #   https://github.com/patrickfuller/blender-chemicals
    # This in 3D in manim is pretty slow (well, sphere rendering is obviously slowest)
    #
    # The API for node_repr and edge_repr should instead take mol and index,
    # then we can do much fancier things.
    #

    def __init__(self,
                 molecule: AllChem.Mol,
                 conformer: int = 0,
                 node_repr=Circle,
                 edge_repr=Line,
                 sub_center: int = None,
                 sub_radius: int = 2,
                 **kwargs):

        super().__init__(**kwargs)

        self.molecule = molecule
        self.conformer = conformer
        self._conformer = self.molecule.GetConformer(self.conformer)
        self.node_repr = node_repr
        self.edge_repr = edge_repr
        self._atom_index_to_node: Dict[int, VMobject] = {}
        self._bond_index_to_edge: Dict[int, VMobject] = {}

        # submolecule info
        self.atom_center = sub_center
        self.radius = sub_radius
        if self.atom_center is not None:
            self.sub_atom_map, self.submol, self.sub_smiles = self.atom_environment(atom_index=self.atom_center,
                                                                                    radius=self.radius)
        else:
            self.sub_atom_map, self.submol, self.sub_smiles = {}, None, ''

        # noinspection PyArgumentList
        for atom_index in range(self.molecule.GetNumAtoms()):
            self.add(self.node(atom_index))
        # noinspection PyArgumentList
        for bond_index in range(self.molecule.GetNumBonds()):
            self.add(self.edge(bond_index))

    def node(self, atom_index: int):
        try:
            return self._atom_index_to_node[atom_index]
        except KeyError:
            atom = self.molecule.GetAtomWithIdx(atom_index)
            pos = self._conformer.GetAtomPosition(atom_index)
            node = self.node_repr()
            node.set_x(pos.x)
            node.set_y(pos.y)
            node.set_z(pos.z)
            node.set_color(ATOM_INFO[atom.GetSymbol()]['color'])
            # node.set_width(0.6 if atom.GetSymbol() == 'C' else 0.8)
            node.set_width(0.5)
            if self.atom_in_submol(atom):
                # node.set_width(1.0)
                node.set_opacity(1)
            if atom_index == self.atom_center:
                node.scale(1.7)
            self._atom_index_to_node[atom_index] = node
        return self._atom_index_to_node[atom_index]

    def edge(self, bond_index: int):
        try:
            return self._bond_index_to_edge[bond_index]
        except KeyError:
            bond = self.molecule.GetBondWithIdx(bond_index)
            node1 = self._conformer.GetAtomPosition(bond.GetBeginAtomIdx())
            node2 = self._conformer.GetAtomPosition(bond.GetEndAtomIdx())
            edge = self.edge_repr()
            edge.set_points_as_corners([[node1.x, node1.y, node1.z],
                                        [node2.x, node2.y, node2.z]])
            if self.bond_in_submol(bond):
                edge.set_stroke(color=SELECTED_SUBSTRUCTURE_COLOR)
            self._bond_index_to_edge[bond_index] = edge
            # FIXME: of course, missing bond type (double and the like)
        return self._bond_index_to_edge[bond_index]

    def atom_in_submol(self, atom: Union[int, AllChem.Atom]):
        if not isinstance(atom, int):
            # noinspection PyArgumentList
            atom = atom.GetIdx()
        return atom in self.sub_atom_map

    def bond_in_submol(self, bond: Union[int, AllChem.Bond]):
        if isinstance(bond, int):
            bond = self.molecule.GetBondWithIdx(bond)
        # noinspection PyArgumentList
        return self.atom_in_submol(bond.GetBeginAtom()) and self.atom_in_submol(bond.GetEndAtom())

    def atom_environment(self, atom_index: int, radius: int):
        # https://www.rdkit.org/docs/GettingStartedInPython.html#explaining-bits-from-morgan-fingerprints
        atom_environment = AllChem.FindAtomEnvironmentOfRadiusN(self.molecule, radius, atom_index, False)
        atom_map = {}  # atom_index_molecule -> atom_index_submolecule
        submol = AllChem.PathToSubmol(self.molecule, atom_environment, atomMap=atom_map)
        smiles = AllChem.MolToSmiles(submol, rootedAtAtom=atom_map[atom_index], canonical=False)
        return atom_map, submol, smiles

    def sub_molecule(self):
        if not self.submol:
            return self
        return Molecule(self.submol,
                        conformer=self.conformer,
                        node_repr=self.node_repr,
                        edge_repr=self.edge_repr,
                        sub_center=self.sub_atom_map[self.atom_center],
                        sub_radius=self.radius)

    def submol_graph(self, copy=True):
        nodes = [deepcopy(node) if copy else node
                 for atom_index, node in self._atom_index_to_node.items()
                 if self.atom_in_submol(atom_index)]
        edges = [deepcopy(edge) if copy else edge
                 for bond_index, edge in self._bond_index_to_edge.items()
                 if self.bond_in_submol(bond_index)]
        return VGroup(*nodes + edges)

    @property
    def smiles(self):
        return AllChem.MolToSmiles(self.molecule)


# --- Scenes

# Good scene from https://talkingphysics.wordpress.com/2018/06/14/creating-text-manim-series-part-4/
# noinspection PyShadowingNames
def einstein_quotes(scene):
    quote = TextMobject("Imagination is more important than knowledge")
    quote.set_color(COLOR_MAP['RED_B'])
    quote.to_edge(UP)
    quote2 = TextMobject("A person who never made a mistake never tried anything new")
    quote2.set_color(COLOR_MAP['YELLOW_E'])
    author = TextMobject("-Albert Einstein")
    author.scale(0.75)
    author.next_to(quote.get_corner(DOWN + RIGHT), DOWN)

    scene.add(quote)
    scene.add(author)
    scene.wait(2)
    scene.play(Transform(quote, quote2),
               ApplyMethod(author.move_to, quote2.get_corner(DOWN + RIGHT) + DOWN + 2 * LEFT))
    scene.play(ApplyMethod(author.match_color, quote2), Transform(author, author.scale(1)))
    scene.wait(2)
    scene.play(FadeOut(quote), FadeOut(author))


class MorganFingerprint(MovingCameraScene):

    def __init__(self, molecule=None, centers=(1, 3, 5, 7), radii=(1, 2, 3), conformer=0, **kwargs):
        if molecule is None:
            molecule = next(AllChem.SDMolSupplier(str(ARTEMISININ_PUBCHEM_2D)))
        self.molecule = molecule
        self.conformer = conformer
        self.centers = centers
        self.radii = radii
        super().__init__(**kwargs)

    def construct(self):

        # Keep track of substructures
        seen_substructures = defaultdict(set)

        # Display the original molecule
        original_molecule = Molecule(self.molecule,
                                     conformer=self.conformer,
                                     node_repr=Circle,
                                     edge_repr=Line,
                                     sub_center=None,
                                     sub_radius=0)
        original_molecule.scale(0.8)
        original_molecule.next_to(LEFT_SIDE, RIGHT)
        molecule_name = TextMobject('Artemisinin(active)', tex_to_color_map={'(active)': ACTIVE_COLOR})
        molecule_name.next_to(original_molecule, DOWN)
        self.play(
            ShowCreation(original_molecule),
            Write(molecule_name)
        )
        self.wait(7)

        # Display fingerprint
        matrix = IntegerMatrix(np.zeros((8, 1), dtype=int))
        matrix.next_to(RIGHT_SIDE, 7 * LEFT)
        matrix_name = TextMobject('Fingerprint (size=8)')
        matrix_name.next_to(matrix, DOWN)
        self.play(Write(matrix), Write(matrix_name))
        self.wait(7)

        # Animate fingerprinting algorithm
        current_molecule = None
        for center in self.centers:

            # Initialize molecule without highl
            if current_molecule is not None:
                # Initialize molecule without highlighting
                original_molecule = Molecule(self.molecule,
                                             conformer=self.conformer,
                                             node_repr=Circle,
                                             edge_repr=Line,
                                             sub_center=None,
                                             sub_radius=0)
                original_molecule.scale(0.8)
                original_molecule.next_to(LEFT_SIDE, RIGHT)
                self.play(ReplacementTransform(current_molecule, original_molecule))

            current_molecule = original_molecule

            for radius in self.radii:

                # Replace current molecule with a version with the current submol highlighted
                highlighted_molecule = Molecule(self.molecule,
                                                conformer=self.conformer,
                                                node_repr=Circle,
                                                edge_repr=Line,
                                                sub_center=center,
                                                sub_radius=radius)
                highlighted_molecule.scale(0.8)
                highlighted_molecule.next_to(LEFT_SIDE, RIGHT)
                radius_text = TextMobject(f'atom {center}, radius {radius}')
                radius_text.next_to(highlighted_molecule, UP)
                self.play(ReplacementTransform(current_molecule, highlighted_molecule),
                          GrowFromCenter(radius_text))
                current_molecule = highlighted_molecule

                # Animate extracting the submol
                submol_origin = highlighted_molecule.submol_graph(copy=True)
                submol_target = highlighted_molecule.submol_graph(copy=True)
                submol_target.scale(0.8)
                submol_target.next_to(highlighted_molecule, 4 * RIGHT)
                self.play(ReplacementTransform(submol_origin, submol_target))

                # Animate assigning the submol to a fingerprint bucket
                # FIXME: actually we should simply use rdkit hash...
                submol_hash = stable_hash(highlighted_molecule.sub_smiles,
                                          fold_to=len(matrix.get_entries()))
                bucket_set = seen_substructures[submol_hash]
                entry = matrix.get_entries()[submol_hash]
                is_new_collision = len(bucket_set) and (highlighted_molecule.sub_smiles not in bucket_set)
                seen_substructures[submol_hash].add(highlighted_molecule.sub_smiles)
                color = NO_COLLISION_COLOR if not is_new_collision else COLLISION_COLOR
                if is_new_collision:
                    collision_text = TextMobject('Collision!', color=COLLISION_COLOR)
                else:
                    collision_text = TextMobject('No collision', color=NO_COLLISION_COLOR)
                collision_text.next_to(matrix, UP)
                submol_to_bucket_arrow = Arrow(submol_target, entry)
                entry.set_value(1)
                entry.set_color(color)
                self.play(
                    entry.set_value, 1,
                    entry.set_color, color,
                    Write(submol_to_bucket_arrow),
                    Write(collision_text),
                    ReplacementTransform(submol_target.copy(), entry)
                )

                # Get rid of all this "bit"
                self.play(FadeOut(submol_to_bucket_arrow),
                          FadeOut(submol_origin),
                          FadeOut(submol_target),
                          FadeOut(radius_text),
                          FadeOut(collision_text))


def transpose_vector(vector: Matrix):
    # FIXME: Nasty in a hurry
    clazz = Matrix.__class__
    config = Matrix.CONFIG
    array = np.vectorize(lambda x: x.value)

    if vector.get_mob_matrix().ndim == 2:
        vector.mob_matrix = vector.get_mob_matrix().flatten()
    else:
        vector.mob_matrix = np.atleast_2d(vector.mob_matrix)

    return vector


# Create fingerprint vectors
def growing_fingerprints():
    matrices = [IntegerMatrix(np.zeros((num_columns, 1), dtype=int))
                for num_columns in (8,)]  # 2, 4, 8, 16, 32, 64, 128, 256
    for matrix in matrices:
        matrix.next_to(RIGHT_SIDE, 5*LEFT)

    # Transform vectors to make a point of number of Collisions
    current_matrix = matrices[0]
    self.play(Write(current_matrix))
    for matrix in matrices[1:]:
        height = (matrix.get_height() if matrix.get_height() > self.camera_frame.get_height()
                  else self.camera_frame.get_height())
        self.play(
            ReplacementTransform(current_matrix, matrix),
            self.camera_frame.set_height, height
        )
        current_matrix = matrix

        # transpose_vector(matrix)
        # matrix.flip()
        # matrix.get_image()
        # matrix.get_mob_matrix()


class FeatureMatrix(MovingCameraScene):

    def construct(self):

        # Some inspiring quotes
        # einstein_quotes(self)

        # Show the original molecule / fingerprint pairs
        original_molecule = Molecule(next(AllChem.SDMolSupplier(str(ARTEMISININ_PUBCHEM_2D))),
                                     conformer=0,
                                     node_repr=Circle,
                                     edge_repr=Line,
                                     sub_center=None,
                                     sub_radius=0)
        original_molecule.scale(0.8)
        original_molecule.next_to(LEFT_SIDE, RIGHT)
        original_molecule_name = TextMobject('Artemisinin(active)', tex_to_color_map={'(active)': ACTIVE_COLOR})
        original_molecule_name.next_to(original_molecule, DOWN)
        original_matrix = IntegerMatrix(np.array([0, 1, 1, 0, 1, 0, 1, 1], dtype=int))
        original_matrix.next_to(RIGHT_SIDE, 7 * LEFT)
        original_matrix_name = TextMobject('Fingerprint (size=8)')
        original_matrix_name.next_to(original_matrix, DOWN)
        self.play(
            ShowCreation(original_molecule),
            Write(original_molecule_name),
            Write(original_matrix),
            Write(original_matrix_name)
        )

        corner_molecule = original_molecule.copy()
        corner_molecule.scale(0.2)
        corner_molecule.to_corner(UP + LEFT)
        corner_molecule_activity = TextMobject('A', color=ACTIVE_COLOR)
        corner_molecule_activity.scale(0.6)
        corner_molecule_activity.next_to(corner_molecule, RIGHT)
        # FIXME: implement transpose for manim matrices
        corner_molecule_fingerprint = IntegerMatrix(np.array([0, 1, 1, 0, 1, 0, 1, 1], dtype=int).reshape(1, -1))
        corner_molecule_fingerprint.set_height(0.7 * corner_molecule.get_height())
        corner_molecule_fingerprint.next_to(corner_molecule_activity, 2 * RIGHT)

        self.play(
            ReplacementTransform(original_molecule, corner_molecule),
            ReplacementTransform(original_molecule_name, corner_molecule_activity),
            ReplacementTransform(original_matrix, corner_molecule_fingerprint),
            FadeOut(original_matrix_name),
        )

        # Exercising matrix transpose

        # self.play(original_molecule.scale, 0.2,
        #           original_molecule.to_corner, LEFT + UP)
        self.wait(1)

        # We will generate some random fingerprints
        rng = np.random.RandomState(0)

        # Add many other molecules
        mols = [corner_molecule]
        mol_activities = [corner_molecule_activity]
        mol_fingerprints = [corner_molecule_fingerprint]
        for smiles, active in SMILES_ACTIVITIES[:10]:

            mol = Molecule(to_rdkit_mol(smiles, to2D=True))
            mol.set_width(corner_molecule.get_width(), stretch=False)
            mol.next_to(mols[-1], DOWN)

            mol_activity = TextMobject('A' if active else 'I',
                                       color=ACTIVE_COLOR if active else INACTIVE_COLOR)
            mol_activity.scale(0.6)
            mol_activity.next_to(mol, RIGHT)

            mol_fingerprint = IntegerMatrix(rng.choice((0, 1), (1, 8)))
            mol_fingerprint.set_height(corner_molecule_fingerprint.get_height())
            mol_fingerprint.next_to(mol_fingerprints[-1], DOWN)

            mols.append(mol)
            mol_activities.append(mol_activity)
            mol_fingerprints.append(mol_fingerprint)
        # noinspection PyTypeChecker

        # Make the whole dataset appear at once
        self.play(FadeIn(VGroup(*(mols[1:] + mol_activities[1:] + mol_fingerprints[1:]))))

        # Grow the fingerprints
        fingerprint_sizes = (16, 32, 128)
        for fingerprint_size in fingerprint_sizes:
            new_fingerprints = []
            for activity in mol_activities:
                mol_fingerprint = IntegerMatrix(rng.choice((0, 1), (1, fingerprint_size)))
                mol_fingerprint.set_height(corner_molecule_fingerprint.get_height())
                new_fingerprints.append(mol_fingerprint)
            new_width = (self.camera_frame.get_width()
                         + new_fingerprints[0].get_width()
                         - mol_fingerprints[0].get_width())
            for activity, mol_fingerprint in zip(mol_activities, new_fingerprints):
                mol_fingerprint.next_to(activity, 2 * RIGHT)
            group = VGroup(*(mols + mol_activities + new_fingerprints))
            group.to_corner(UP + LEFT)
            self.play(
                self.camera_frame.set_width, new_width,
                *[ReplacementTransform(src, dst) for src, dst in zip(mol_fingerprints, new_fingerprints)],
            )


            mol_fingerprints = new_fingerprints


class LogisticRegression(Scene):
    # State what is learned weights, and the strong, simple, linear bias
    # (use swedish example from "How not to be wrong")
    ...


class RepresentationLearning(Scene):
    # Mapping real world objects to Riemannian manifold...
    ...


class Molecule3D(ThreeDScene):
    ...


class Malaria(Scene):

    #
    # Our world in data:
    #   https://ourworldindata.org/malaria
    #   https://ourworldindata.org/causes-of-death
    #

    def construct(self):
        ...


class Eroom(Scene):
    #
    # We could also use things like the Morph transition in powerpoint
    #
    # Data and comments:
    #   - Moore's Law:
    #     https://github.com/wallento/mooreandmore
    #     https://ourworldindata.org/technological-progress
    #     https://towardsdatascience.com/moores-law-is-dying-here-s-how-ai-is-bringing-it-back-to-life-c9a469bc7a5a
    #     https://en.wikipedia.org/wiki/Moore%27s_law
    #     https://www.forbes.com/sites/danwoods/2013/12/12/how-to-create-a-moores-law-for-data/#711eb9eb44ca
    #     https://towardsdatascience.com/moores-law-is-dead-678119754571
    #     https://medium.com/predict/moores-law-is-alive-and-well-adc010ea7a63
    #
    #   - Eroom's Law:
    #     https://blogs.scientificamerican.com/observations/how-to-fight-erooms-law/
    #     https://new.pharmacelera.com/publications/what-is-erooms-law/
    #     https://ourworldindata.org/technological-progress
    #
    def construct(self):
        pass


# --- Entry point

if __name__ == "__main__":

    #
    # We should refactor manim "main" this to accept args
    #   https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
    # In the meantime, we trick it via sys.argv hijacking so we can debug with no problem...
    #

    import sys
    from manimlib import main
    from pathlib import Path

    argv = sys.argv
    try:
        media_dir = Path(__file__).parent.parent / 'media'
        video_dir = media_dir / 'video'
        tex_dir = media_dir / 'tex'
        scenes = (
            # AddingMoreText,
            # MorganFingerprint,
            FeatureMatrix,
        )
        low_quality = True
        for scene in scenes:
            sys.argv = ['manim',
                        '-pl' if low_quality else '-p',
                        '--video_dir', str(video_dir),
                        '--tex_dir', str(tex_dir),
                        __file__,
                        scene.__name__]
            main()
    finally:
        sys.argv = argv


# --- Braindump

# Move camera in scene 3D
# self.move_camera(0.8*np.pi/2, -0.45*np.pi)
# self.begin_ambient_camera_rotation()

#
# Primer, blender and molecules
#   https://github.com/Helpsypoo/primer
#   https://patrickfuller.github.io/molecules-from-smiles-molfiles-in-blender/
#   And there are many molecular visualization tools based on blender
# What about good old POVRay? Missing those days...
#
