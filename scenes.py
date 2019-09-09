from collections import Counter, defaultdict
from copy import deepcopy
from typing import Dict, Union
import hashlib

import numpy as np
from manimlib.animation.fading import FadeOut, VFadeOut
from manimlib.animation.transform import ApplyFunction, TransformFromCopy
from manimlib.imports import Scene, Circle, ShowCreation, Transform, \
    COLOR_MAP, Line, Write, MovingCameraScene, VMobject, RIGHT, VGroup, LEFT, LEFT_SIDE, RIGHT_SIDE, \
    GrowFromCenter, DOWN, IntegerMatrix, Arrow, ReplacementTransform, UP, ThreeDScene
from manimlib.mobject.svg.tex_mobject import TextMobject
from rdkit.Chem import AllChem

#
# ARTEMISININ
#   https://pubchem.ncbi.nlm.nih.gov/compound/68827
#   https://en.wikipedia.org/wiki/Artemisinin
#
ARTEMISININ_SMILES = 'CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C'
ARTEMISININ_ISOMERIC_SMILES = 'C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C'
ARTEMISININ_INCHI = 'InChI=1S/C15H22O5/c1-8-4-5-11-9(2)12(16)17-13-15(11)10(8)6-7-14(3,18-13)19-20-15/' \
                    'h8-11,13H,4-7H2,1-3H3/t8-,9-,10+,11+,13-,14-,15-/m1/s1'

#
# Primer, blender and molecules
#   https://github.com/Helpsypoo/primer
#   https://patrickfuller.github.io/molecules-from-smiles-molfiles-in-blender/
#   And there are many molecular visualization tools based on blender
# What about POVRay
#


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


# self.mol = AllChem.MolFromSmiles(smiles, sanitize=True)
# AllChem.AddHs(self.mol)
# conformer = AllChem.EmbedMolecule(self.mol, randomSeed=self.seed)
# AllChem.MMFFOptimizeMolecule(self.mol)
# if -1 == conformer:
#     raise Exception('Could not embed molecule in 3D using ETDG ('
#                     'stereochemistry? try using non-isomeric smiles)')
# AllChem.RemoveHs(self.mol)
# conformer = AllChem.Compute2DCoords(self.mol)
# self.mol = next(AllChem.SDMolSupplier('data/artemisinin/Conformer3D_CID_68827.sdf'))
# self.mol = next(AllChem.SDMolSupplier('data/artemisinin/Structure2D_CID_68827.sdf'))


ATOM_COLORS = {
    'C': COLOR_MAP['DARK_GREY'],
    'O': COLOR_MAP['RED_B']
}


def stable_hash(string, hasher=hashlib.md5, fold_to=None):
    shash = int(hasher(string.encode('utf-8')).hexdigest(), 16)
    return shash if not fold_to else shash % fold_to


class Molecule(VMobject):

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
            node.set_color(ATOM_COLORS[atom.GetSymbol()])
            # node.set_width(0.6 if atom.GetSymbol() == 'C' else 0.8)
            node.set_width(0.5)
            if self.atom_in_submol(atom):
                # node.set_width(1.0)
                node.set_opacity(1)
            if atom_index == self.atom_center:
                node.scale(2)
                # node.set_color(COLOR_MAP['YELLOW_E'])
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
                edge.set_stroke(color=COLOR_MAP['YELLOW_E'])
            self._bond_index_to_edge[bond_index] = edge
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
        return VGroup(*nodes+edges)

    @property
    def smiles(self):
        return AllChem.MolToSmiles(self.molecule)


class MorganFingerprintScene(MovingCameraScene):

    def __init__(self, molecule=None, centers=(1, 3, 5, 7), radii=(1, 2, 3), conformer=0, **kwargs):
        if molecule is None:
            molecule = next(AllChem.SDMolSupplier('data/artemisinin/Structure2D_CID_68827.sdf'))
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
        molecule_name = TextMobject('Artemisinin')
        molecule_name.next_to(original_molecule, UP)
        self.play(
            ShowCreation(original_molecule),
            Write(molecule_name)
        )
        self.wait(6)

        # Create fingerprint vectors
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
        matrix = IntegerMatrix(np.zeros((8, 1), dtype=int))
        matrix.next_to(RIGHT_SIDE, 5*LEFT)
        self.play(ReplacementTransform(current_matrix, matrix))
        self.wait(6)

        current_matrix = matrix
        current_molecule = None

        for center in self.centers:

            if current_molecule is None:
                self.play(FadeOut(molecule_name))
            else:
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
                radius_text = TextMobject(f'radius={radius}')
                radius_text.next_to(highlighted_molecule, UP)
                self.play(ReplacementTransform(current_molecule, highlighted_molecule),
                          GrowFromCenter(radius_text))
                current_molecule = highlighted_molecule

                # Animate extracting the submol
                submol_origin = highlighted_molecule.submol_graph(copy=True)
                submol_target = highlighted_molecule.submol_graph(copy=True)
                submol_target.scale(0.8)
                submol_target.next_to(highlighted_molecule, 6 * RIGHT)
                self.play(ReplacementTransform(submol_origin, submol_target))

                # Animate assigning the submol to a fingerprint bucket
                submol_hash = stable_hash(highlighted_molecule.sub_smiles,
                                          fold_to=len(current_matrix.get_entries()))
                bucket_set = seen_substructures[submol_hash]
                entry = current_matrix.get_entries()[submol_hash]
                is_new_collision = len(bucket_set) and (highlighted_molecule.sub_smiles not in bucket_set)
                seen_substructures[submol_hash].add(highlighted_molecule.sub_smiles)
                color = COLOR_MAP['GREEN_E'] if not is_new_collision else COLOR_MAP['RED_E']
                if is_new_collision:
                    collision_text = TextMobject('Collision!', color=COLOR_MAP['RED_E'])
                else:
                    collision_text = TextMobject('No collision', color=COLOR_MAP['GREEN_E'])
                collision_text.next_to(current_matrix, UP)
                arrow = Arrow(submol_target, entry)
                entry.set_value(1)
                entry.set_color(color)
                self.play(
                    entry.set_value, 1,
                    entry.set_color, color,
                    Write(arrow),
                    Write(collision_text),
                    ReplacementTransform(submol_target.copy(), entry)
                )

                # Get rid of all this "bit"
                self.play(FadeOut(arrow),
                          FadeOut(submol_origin),
                          FadeOut(submol_target),
                          FadeOut(radius_text),
                          FadeOut(collision_text))


class FeatureMatrix(Scene):
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


if __name__ == "__main__":

    #
    # We should refactor this to accept args
    # https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
    #
    # In the meantime, we trick it so we can debug with no problem...
    #

    import sys
    from manimlib import main
    from pathlib import Path

    argv = sys.argv
    try:
        media_dir = Path(__file__).parent / 'media'
        video_dir = media_dir / 'video'
        tex_dir = media_dir / 'tex'
        sys.argv = ['manim',
                    '-pl',
                    '--video_dir', str(video_dir),
                    '--tex_dir', str(tex_dir),
                    __file__,
                    MorganFingerprintScene.__name__]
        main()
    finally:
        sys.argv = argv


# TODO: maybe translate the substructure instead of making it appear

# Some other fancy things we can do
# self.play(GrowFromCenter(molecule))
# self.play(ApplyMethod(molecule.shift, 3 * DOWN))
# self.play(FadeOut(molecule))

# Move camera in scene 3D
# self.move_camera(0.8*np.pi/2, -0.45*np.pi)
# self.begin_ambient_camera_rotation()
