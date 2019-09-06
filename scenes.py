import numpy as np
from manimlib.imports import Scene, Dot, Circle, Square, Polygon, ShowCreation, FadeOut, GrowFromCenter, Transform, \
    Mobject
from rdkit import Chem
from rdkit.Chem import AllChem

ARTEMISININ_SMILES = 'O=C3O[C@@H]4O[C@@]1(OO[C@@]42[C@@H](CC1)[C@H](C)CC[C@H]2[C@H]3C)C'


class Molecule(Mobject):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        mol = Chem.MolFromSmiles(ARTEMISININ_SMILES)
        conformer = AllChem.Compute2DCoords(mol)
        conformer = mol.GetConformer(conformer)
        for atom in mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            circle = Circle()
            circle.set_x(pos.x)
            circle.set_y(pos.y)
            self.add(circle)


class ExampleScene(Scene):
    def construct(self):
        molecule = Molecule()
        self.play(ShowCreation(molecule))
        self.play(FadeOut(molecule))


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
