from copy import deepcopy

import numpy as np
from manimlib.animation.transform import ApplyMethod
from manimlib.imports import Scene, Dot, Circle, Square, Polygon, ShowCreation, FadeOut, GrowFromCenter, Transform, \
    Mobject, COLOR_MAP, Line, ThreeDScene, DOWN, Sphere
from rdkit.Chem import AllChem

# https://pubchem.ncbi.nlm.nih.gov/compound/68827
# https://en.wikipedia.org/wiki/Artemisinin
ARTEMISININ_SMILES = 'CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C'
ARTEMISININ_ISOMERIC_SMILES = 'C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C'
ARTEMISININ_INCHI = 'InChI=1S/C15H22O5/c1-8-4-5-11-9(2)12(16)17-13-15(11)10(8)6-7-14(3,18-13)19-20-15/' \
                    'h8-11,13H,4-7H2,1-3H3/t8-,9-,10+,11+,13-,14-,15-/m1/s1'


#
# The discovery of artemisinin and Nobel Prize in Physiology or Medicine:
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4966551/
# We could also work with dihydroartemisinin
#
# Tu Youyou: https://www.nobelprize.org/prizes/medicine/2015/tu/facts/
#

# mol = next(AllChem.SDMolSupplier('data/artemisinin/Conformer3D_CID_68827.sdf'))
# print(mol)
# AllChem.AddHs(mol)
# conformer = AllChem.EmbedMolecule(mol)
# print(conformer)
# exit(22)


class Molecule(Mobject):

    def __init__(self,
                 smiles=ARTEMISININ_SMILES,
                 seed=0,
                 **kwargs):

        super().__init__(**kwargs)

        self.seed = seed

        # self.mol = AllChem.MolFromSmiles(smiles, sanitize=True)
        # AllChem.AddHs(self.mol)
        # conformer = AllChem.EmbedMolecule(self.mol, randomSeed=self.seed)
        # AllChem.MMFFOptimizeMolecule(self.mol)
        # if -1 == conformer:
        #     raise Exception('Could not embed molecule in 3D using ETDG ('
        #                     'stereochemistry? try using non-isomeric smiles)')
        # AllChem.RemoveHs(self.mol)
        # conformer = AllChem.Compute2DCoords(self.mol)
        self.mol = next(AllChem.SDMolSupplier('data/artemisinin/Conformer3D_CID_68827.sdf'))
        # self.mol = next(AllChem.SDMolSupplier('data/artemisinin/Structure2D_CID_68827.sdf'))
        conformer = self.mol.GetConformer(0)
        for atom in self.mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            # circle = Circle()
            # circle.set_points([[pos.x, pos.y, pos.z]])
            circle = Sphere()
            circle.set_x(pos.x)
            circle.set_y(pos.y)
            circle.set_z(pos.z)
            circle.set_color(COLOR_MAP['DARK_BLUE'] if atom.GetSymbol() == 'C' else COLOR_MAP['RED_B'])
            circle.set_width(0.6 if atom.GetSymbol() == 'C' else 0.8)
            self.add(circle)
        for bond in self.mol.GetBonds():
            pos1 = conformer.GetAtomPosition(bond.GetBeginAtomIdx())
            pos2 = conformer.GetAtomPosition(bond.GetEndAtomIdx())
            line = Line()
            line.set_points_as_corners([[pos1.x, pos1.y, pos1.z],
                                        [pos2.x, pos2.y, pos2.z]])
            self.add(line)

    @property
    def smiles(self):
        return AllChem.MolToSmiles(self.mol)


class ExampleScene(ThreeDScene):

    def construct(self):
        molecule = Molecule()
        # self.play(ShowCreation(molecule))
        # self.play(FadeOut(molecule))
        self.camera.set_frame_center(molecule.get_center())
        self.play(GrowFromCenter(molecule))
        self.play(ApplyMethod(molecule.shift, 3 * DOWN))
        self.move_camera(0.8*np.pi/2, -0.45*np.pi)
        self.begin_ambient_camera_rotation()
        self.wait(6)


class Shapes(Scene):

    def construct(self):
        circle = Circle()
        square = Square()
        triangle = Polygon(np.array([0, 0, 0]), np.array([1, 1, 0]), np.array([1, -1, 0]))

        self.play(ShowCreation(circle))
        self.play(FadeOut(circle))
        self.play(GrowFromCenter(square))
        self.play(Transform(square, triangle))


if __name__ == "__main__":
    import os
    module_name = os.path.basename(__file__)
    command_A = "manim -pl --video_dir ~/Downloads/  "
    command_B = module_name + " " + "ExampleScene"
    os.system(command_A + command_B)
