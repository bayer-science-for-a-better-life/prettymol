import hashlib
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Dict, Union

import numpy as np
from manimlib.animation.fading import FadeOut
from manimlib.imports import Scene, Circle, ShowCreation, COLOR_MAP, Line, Write, MovingCameraScene, VMobject, RIGHT, \
    VGroup, LEFT, LEFT_SIDE, RIGHT_SIDE, \
    GrowFromCenter, DOWN, IntegerMatrix, Arrow, ReplacementTransform, UP, ThreeDScene, Transform, ApplyMethod
from manimlib.mobject.svg.tex_mobject import TextMobject
from rdkit.Chem import AllChem

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
ARTEMISININ_PUBCHEM_SMILES = 'CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C'
ARTEMISININ_PUBCHEM_ISOMERIC_SMILES = 'C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C'
ARTEMISININ_INCHI = 'InChI=1S/C15H22O5/c1-8-4-5-11-9(2)12(16)17-13-15(11)10(8)6-7-14(3,18-13)19-20-15/' \
                    'h8-11,13H,4-7H2,1-3H3/t8-,9-,10+,11+,13-,14-,15-/m1/s1'
ARTEMISININ_PUBCHEM_2D = Path(__file__).parent / 'data' / 'artemisinin' / 'Structure2D_CID_68827.sdf'
ARTEMISININ_PUBCHEM_3D = Path(__file__).parent / 'data' / 'artemisinin' / 'Conformer3D_CID_68827.sdf'
ARTEMISININ_PUBCHEM_JSON = Path(__file__).parent / 'data' / 'artemisinin' / 'compound_CID_68827.json'

SMILES_ACTIVITIES = (
    ('CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C', True),
    ('FC(F)(F)c1ccc(nc1)N2CCOCC2', False),
    ('CCc1nnc(NC(=O)C(C)Sc2ccccc2)s1', False),
)

# --- Aesthetics

NO_COLLISION_COLOR = COLOR_MAP['GREEN_E']
COLLISION_COLOR = COLOR_MAP['RED_E']
SELECTED_SUBSTRUCTURE_COLOR = COLOR_MAP['YELLOW_E']
ACTIVE_COLOR = COLOR_MAP['GREEN_E']


def rgb2hex(r, g, b, relative=True):
    if relative:
        r = int(round(r * 255))
        g = int(round(g * 255))
        b = int(round(b * 255))
    return "#{0:02x}{1:02x}{2:02x}".format(max(0, min(r, 255)),
                                           max(0, min(g, 255)),
                                           max(0, min(b, 255)))


#
# Ripped from https://github.com/patrickfuller/blender-chemicals/blob/master/blender_chemicals/atoms.json
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

class AddingMoreText(Scene):
    # Good scene from https://talkingphysics.wordpress.com/2018/06/14/creating-text-manim-series-part-4/
    def construct(self):
        quote = TextMobject("Imagination is more important than knowledge")
        quote.set_color(COLOR_MAP['RED_B'])
        quote.to_edge(UP)
        quote2 = TextMobject("A person who never made a mistake never tried anything new")
        quote2.set_color(COLOR_MAP['YELLOW_E'])
        author = TextMobject("-Albert Einstein")
        author.scale(0.75)
        author.next_to(quote.get_corner(DOWN + RIGHT), DOWN)

        self.add(quote)
        self.add(author)
        self.wait(2)
        self.play(Transform(quote, quote2),
                  ApplyMethod(author.move_to, quote2.get_corner(DOWN + RIGHT) + DOWN + 2 * LEFT))
        self.play(ApplyMethod(author.match_color, quote2), Transform(author, author.scale(1)))
        self.play(FadeOut(quote))


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
        self.wait(1)

        # Display fingerprint
        matrix = IntegerMatrix(np.zeros((8, 1), dtype=int))
        matrix.next_to(RIGHT_SIDE, 7 * LEFT)
        matrix_name = TextMobject('Fingerprint (size=8)')
        matrix_name.next_to(matrix, DOWN)
        self.play(Write(matrix), Write(matrix_name))
        self.wait(1)

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


class FeatureMatrix(Scene):
    ...


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
        media_dir = Path(__file__).parent / 'media'
        video_dir = media_dir / 'video'
        tex_dir = media_dir / 'tex'
        scenes = (
            # AddingMoreText,
            MorganFingerprint,
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

# TODO: maybe translate the substructure instead of making it appear

# Some other fancy things we can do
# self.play(GrowFromCenter(molecule))
# self.play(ApplyMethod(molecule.shift, 3 * DOWN))
# self.play(FadeOut(molecule))

# Move camera in scene 3D
# self.move_camera(0.8*np.pi/2, -0.45*np.pi)
# self.begin_ambient_camera_rotation()

#
# Primer, blender and molecules
#   https://github.com/Helpsypoo/primer
#   https://patrickfuller.github.io/molecules-from-smiles-molfiles-in-blender/
#   And there are many molecular visualization tools based on blender
# What about POVRay?
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
#

# Create fingerprint vectors
# def growing_fingerprints():
#     matrices = [IntegerMatrix(np.zeros((num_columns, 1), dtype=int))
#                 for num_columns in (8,)]  # 2, 4, 8, 16, 32, 64, 128, 256
#     for matrix in matrices:
#         matrix.next_to(RIGHT_SIDE, 5*LEFT)
#
#     # Transform vectors to make a point of number of Collisions
#     current_matrix = matrices[0]
#     self.play(Write(current_matrix))
#     for matrix in matrices[1:]:
#         height = (matrix.get_height() if matrix.get_height() > self.camera_frame.get_height()
#                   else self.camera_frame.get_height())
#         self.play(
#             ReplacementTransform(current_matrix, matrix),
#             self.camera_frame.set_height, height
#         )
#         current_matrix = matrix
