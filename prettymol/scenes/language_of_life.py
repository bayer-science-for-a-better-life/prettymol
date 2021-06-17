import string
from functools import partial
from itertools import chain
from pathlib import Path
from typing import Union, Optional, Tuple

import numpy as np

from matplotlib.colors import to_hex
from logomaker.src.colors import get_color_dict

from manim import (
    Scene, ZoomedScene, MovingCameraScene,
    Mobject, VMobject, SVGMobject, VGroup, Text, SurroundingRectangle,
    UP, LEFT, RIGHT, DOWN, ORIGIN,
    RED, BLUE, GREEN, WHITE, YELLOW, TEAL, GOLD, ORANGE, MAROON, PURPLE, PINK,
    BarChart,
    Write, ReplacementTransform, Transform, FadeIn,
    Wiggle, Indicate, Circumscribe, ApplyWave, FocusOn,
    SMALL_BUFF, MED_LARGE_BUFF,
    BOLD, ShrinkToCenter,
)

from prettymol.config import Config
from prettymol.manim_utils import manimce, is_manimce


# --- Utils

def get_aa_hex_color_dict(name='dmslogo_funcgroup', seed=42):
    """A weird overcomplicated function mapping ASCII letters to colors maybe relating to some amino acid property."""

    NOT_AA = 'BJOUXZ'
    AAS = ''.join(c for c in string.ascii_uppercase if c not in NOT_AA)

    aa_color_dict = {char: to_hex(color) for char, color in get_color_dict(name, AAS).items()}

    not_aa_color_dict = {}
    if seed is not None:
        rng = np.random.RandomState(seed=seed)
        for c in NOT_AA:
            not_aa_color_dict[c] = rng.choice(list(aa_color_dict.values()))
    else:
        for c in NOT_AA:
            not_aa_color_dict[c] = to_hex(0)

    aa_color_dict.update(not_aa_color_dict)

    final_color_dict = {}
    for c, color in sorted(aa_color_dict.items()):
        if c not in final_color_dict:
            final_color_dict[c.upper()] = final_color_dict[c.lower()] = color

    return final_color_dict


# --- Media library

XARELTO_SMILES = 'C1COCC(=O)N1C2=CC=C(C=C2)N3CC(OC3=O)CNC(=O)C4=CC=C(S4)Cl'
ADENINE_SMILES = 'C1=NC2=NC=NC(=C2N1)N'


class SVGS:
    DNA = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'dna.svg')
    RNA = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'rna.svg')
    PROTEIN2D = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'protein2D.svg')
    PROTEIN3D = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'protein3D.svg')
    PROTEIN_AND_SUGAR = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'protein-sugar.svg')
    ANTIBODY = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'antibody-simplified.svg')
    PACMAN = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'pacman.svg')
    CYCLIC_PEPTIDE = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'peptide-simplified.svg')
    PILL = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'pill.svg')
    PLANT = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'plant-simplified.svg')
    CELL = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'cell-simplified.svg')
    CELLS = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'cells.svg')
    MUTATION = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'mutation.svg')
    PATIENT = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'patient.svg')
    MODEL = str(Config.DEFAULT_IMAGES_PATH / 'lol' / 'chip-brain.svg')


# --- MObjects

# noinspection PyAbstractClass
class LoLLogo(VGroup):

    def __init__(self, commons_left=SVGS.DNA, commons_right=None, **kwargs):

        super().__init__(**kwargs)

        # Logo level
        lol = Text('Language of Life')
        lol[1].set_color(BLUE)
        if is_manimce(self):
            # Beware: currently semantics seem to have changed in manim ce and this should be 11
            #         are they not counting spaces?
            lol[11].set_color(GREEN)
        else:
            lol[13].set_color(GREEN)

        commons = Text('commons')
        commons.next_to(lol, DOWN)

        self.lol = lol
        self.commons = commons
        self.add(self.lol, self.commons,)

        self.commons_left = self.commons_right = None
        self._set_left_right(commons_left, commons_right)

    def add_to_scene(self, scene, initial_blank_seconds=0., end_blank_seconds=0.):
        scene.wait(initial_blank_seconds)
        scene.play(FadeIn(self.lol))
        scene.wait()
        scene.play(Write(self.commons), FadeIn(self.commons_left), FadeIn(self.commons_right))
        scene.wait(end_blank_seconds)

    def _set_left_right(self,
                        commons_left: Union[VMobject, str, Path],
                        commons_right: Optional[Union[VMobject, str, Path]] = None,
                        scale=0.4,
                        left_colors=(GREEN, BLUE),
                        right_colors=(BLUE, GREEN),
                        stroke_width=0.5):

        if self.commons_left is not None:
            self.remove(self.commons_left, self.commons_right)

        if commons_right is None:
            commons_right = commons_left

        if not isinstance(commons_left, Mobject):
            commons_left = (SVGMobject(commons_left).
                            scale(scale).
                            set_color_by_gradient(*left_colors).
                            set_stroke(width=stroke_width))
        commons_left.next_to(self.commons, LEFT)

        if not isinstance(commons_right, Mobject):
            commons_right = (SVGMobject(commons_right).
                             scale(scale).
                             set_color_by_gradient(*right_colors).
                             set_stroke(width=stroke_width))
        commons_right.next_to(self.commons, RIGHT)

        self.commons_left = commons_left
        self.commons_right = commons_right

        self.add(self.commons_left, self.commons_right)

    def replace_left_right(self, scene, commons_left, commons_right=None, run_time_s: Union[float, int] = 1):
        old_left = self.commons_left
        old_right = self.commons_right
        self._set_left_right(commons_left, commons_right)
        scene.play(
            ReplacementTransform(old_left, self.commons_left),
            ReplacementTransform(old_right, self.commons_right),
            run_time=run_time_s
        )


# --- Scenes

class LoLLogoScene(Scene):

    def __init__(self, initial_blank_seconds=2.5, **kwargs):
        super().__init__(**kwargs)
        self.initial_blank_seconds = initial_blank_seconds

    def construct(self):
        LoLLogo().add_to_scene(self,
                               initial_blank_seconds=self.initial_blank_seconds,
                               end_blank_seconds=self.initial_blank_seconds)


class LoLCommonsIntroScene(Scene):
    DEFAULT_REPLACER_SEQUENCE = (SVGS.RNA,
                                 SVGS.CYCLIC_PEPTIDE,
                                 # SVGS.PROTEIN2D,
                                 SVGS.PACMAN,
                                 SVGS.PROTEIN3D,
                                 SVGS.CELL,
                                 SVGS.PLANT,
                                 # SVGS.PILL,
                                 SVGS.PATIENT,
                                 SVGS.DNA)

    def __init__(self,
                 replacers=DEFAULT_REPLACER_SEQUENCE,
                 add_vector_transform=False,
                 logo_down=False,
                 **kwargs):
        super().__init__(**kwargs)
        self.replacers = replacers
        self.add_vector_transform = add_vector_transform
        self.logo_down = logo_down

    def construct(self):

        # Atom level
        adenine = Text(ADENINE_SMILES)

        # Nucleotide level
        dna = Text('ACTGAATATAGACTATA', t2c={'A': RED, 'C': BLUE, 'T': GREEN, 'G': PURPLE})

        # Amino acid level
        aas = 'VQGGAAVQQEVLA'
        aas = Text(aas, t2c=get_aa_hex_color_dict('skylign_protein'))

        # Logo level
        logo = LoLLogo()

        # Animate... from nucleotide to LoL Commons
        self.wait()
        self.play(Write(adenine))
        self.wait()

        if self.add_vector_transform:
            # FIXME: This is not quite pretty / useful at the moment
            numbers = Text('[0.12 0.32 0.01 0.88 0.99 0.22 0.55 0.42 0.13 0.07]')
            replacements = [adenine, dna, aas, numbers, logo.lol]
        else:
            replacements = [adenine, dna, aas, logo.lol]

        for source, target in zip(replacements, replacements[1:]):
            self.play(ReplacementTransform(source, target))
            self.wait()

        self.play(Write(logo.commons), FadeIn(logo.commons_left), FadeIn(logo.commons_right))
        self.wait()
        for replacer in self.replacers:
            logo.replace_left_right(self, replacer, run_time_s=0.8)
            self.wait(0.3)

        # Move logo to left lower corner
        if self.logo_down:
            mini_logo = logo.copy()
            mini_logo.scale(0.5)
            mini_logo.to_corner()
            self.play(Transform(logo, mini_logo))

        # And keep the scene
        # self.play(Write(Text('Corporate info follows...')))
        self.wait(5)

        # FIXME: also color smiles (see compounds scenes)
        # FIXME: also color nucleotides
        # FIXME: also make dna and protein sequences meaningful
        # FIXME: does it make sense to transform adenine only in these bases?
        # FIXME: does it make sense to transform codons into appropriate amino acids?


class EroomScene(Scene):

    def construct(self):

        a_million_forks = Text('A million forks in the road to our products')  # to drug development
        a_million_forks.to_edge(UP)

        new_drugs_per_billion_RD = [
            ('1950', 100),
            ('1960', 50),
            ('1970', 10),
            ('1980', 5),
            ('1990', 3),
            ('2000', 2),
            ('2010', 1),
            ('2020', 1)
        ]
        x = [x for x, _ in new_drugs_per_billion_RD]
        y = [y for _, y in new_drugs_per_billion_RD]
        colors = ["#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"]
        chart = BarChart(
            values=y,
            max_value=max(y),
            bar_colors=colors,
            bar_names=x,
            bar_label_scale_val=0.5,
        ).scale(0.8).to_edge(LEFT).shift(0.5 * DOWN)

        # text_top = (
        #     Text('EROOM\'s Law: More expensive, slower drug discovery')
        #     .scale(0.9)
        #     .next_to(chart, UP, buff=0.1)
        # )

        # text_left = (
        #     Text('Number of drugs per billion US$ R&D spending')
        #     .rotate(angle=TAU / 4, axis=OUT)
        #     .scale(0.3)
        #     .next_to(chart, LEFT, buff=0.5)
        # )

        text_top = (
            Text('Number of drugs per billion US$ R&D spending (log-scale)')
            .scale(0.3)
            .next_to(chart, UP, buff=0.5)
        )

        text_bottom = (
            Text("Eroom's law: a continuous decline in Pharma R&D productivity")
            .scale(0.3)
            .next_to(chart, DOWN, buff=0.5)
        )

        # --- Examples of forks in the road

        icons_text = [
            (
                SVGMobject(SVGS.ANTIBODY).set_color(BLUE).scale(0.8),
                Text('Is my antibody a potent, functional binder?')
            ),
            (
                SVGMobject(SVGS.CYCLIC_PEPTIDE).set_color_by_gradient(BLUE, GREEN).scale(0.8),
                Text('Can we generate peptides without liabilities?')
            ),
            (
                SVGMobject(SVGS.DNA).set_color_by_gradient(GREEN, BLUE).scale(0.8),
                Text('Can we optimize DNA to better express traits in plants or humans?')
            ),
            (
                SVGMobject(SVGS.PACMAN).set_color(RED).scale(0.8),
                Text('What will an enzyme in the human gut do?')
            ),
            # (
            #     SVGMobject(SVGS.PACMAN).set_color_by_gradient(GREEN, BLUE, RED).scale(0.8),
            #     Text('Can we optimize enzymes to better manufacture bioproducts?')
            # ),
            (
                SVGMobject(SVGS.PROTEIN3D).set_color_by_gradient(GREEN, BLUE, RED).scale(0.8),
                Text('What is the 3D structure of my biomolecule?')
            ),
            (
                SVGMobject(SVGS.PATIENT).set_color_by_gradient(BLUE, RED).scale(0.8),
                Text('Will a patient respond to treatment in a clinical trial?')
            ),
            (
                SVGMobject(SVGS.MUTATION).set_color_by_gradient(BLUE, GREEN, RED).scale(0.8),
                Text('...')
            )
        ]

        # icons_text = icons_text[:1]

        for (icon1, _), (icon2, _) in zip(icons_text, icons_text[1:]):
            icon2.next_to(icon1, DOWN, buff=0.5)

        for icon, text in icons_text:
            text.next_to(icon, RIGHT)

        questions = VGroup(*chain.from_iterable(icons_text))
        questions.scale(0.32)
        questions.next_to(chart, RIGHT, buff=1).shift(0.5 * UP)

        logo = LoLLogo().scale(0.7)
        # logo.next_to(questions, DOWN, buff=0.5)
        logo.move_to(questions)

        # --- Animate

        self.play(FadeIn(text_top))
        self.play(Write(chart), Write(text_bottom), run_time=4)
        self.wait(2)
        self.play(Write(a_million_forks))
        for icon, text in icons_text:
            self.play(FadeIn(icon, text))
            self.wait(3)

        # if we have put the logo under the questions...
        # self.play(FadeIn(logo))

        # if we prefer the logo to morph from the questions
        actionable_insights = (Text('Getting actionable insights from biological sequences')
                               .next_to(logo, DOWN, buff=0.5)
                               .scale(0.4))
        self.play(ReplacementTransform(questions, logo))
        self.play(Write(actionable_insights))
        self.wait(2)


# noinspection PyAbstractClass
class UseCase(VGroup):

    def __init__(self, icon=SVGS.ANTIBODY, name='Antibody binding\nprediction', **kwargs):
        super().__init__(**kwargs)
        self.icon = icon
        self.name = name
        if not isinstance(self.icon, Mobject):
            self.icon = SVGMobject(self.icon).set_color_by_gradient(BLUE).scale(0.6)
        if not isinstance(self.name, Mobject):
            self.name = (VGroup(*(Text(name) for name in self.name.splitlines()))
                         .scale(0.5)
                         .arrange(DOWN, center=True, buff=SMALL_BUFF))
        self.name.next_to(self.icon, DOWN)

        self.icon_name = VGroup(self.icon, self.name)
        self.frame = SurroundingRectangle(self.icon_name).round_corners(0.5).set_color(WHITE)

        self.add(self.frame, self.icon_name)


USE_CASES = (

    UseCase(
        SVGMobject(SVGS.ANTIBODY).set_color(BLUE).scale(0.6),
        'Antibody Binding\nPrediction'
    ),
    UseCase(
        SVGMobject(SVGS.PROTEIN3D).set_color(GREEN).scale(0.6),
        'Toxin\n Design'
    ),
    UseCase(
        SVGMobject(SVGS.PACMAN).set_color_by_gradient(BLUE, GREEN).scale(0.6),
        'Enzyme\nStabilisation'
    ),
    UseCase(
        SVGMobject(SVGS.DNA).set_color_by_gradient(GREEN, BLUE).scale(0.6),
        'Codon\nOptimization'
    ),
    UseCase(
        SVGMobject(SVGS.PATIENT).set_color_by_gradient(RED, BLUE).scale(0.6),
        'Clinical Trial\nSimulation'
    ),

    UseCase(
        SVGMobject(SVGS.CYCLIC_PEPTIDE).set_color_by_gradient(GREEN, BLUE).scale(0.6),
        'Peptides for PPI\nDisruption'
    ),

    UseCase(
        SVGMobject(SVGS.CYCLIC_PEPTIDE).set_color_by_gradient(BLUE, GREEN).scale(0.6),
        'Peptides Solubility\nOptimization'
    ),

    UseCase(
        SVGMobject(SVGS.PROTEIN3D).set_color_by_gradient(GREEN, RED, BLUE).scale(0.6),
        'Microbiome\nScreening'
    ),
    UseCase(
        SVGMobject(SVGS.PROTEIN3D).set_color(GREEN).scale(0.6),
        'Herbicide\nDetoxification'
    ),
    UseCase(
        SVGMobject(SVGS.DNA).set_color_by_gradient(BLUE, GREEN, RED).scale(0.6),
        'Mutation Effect\nPrediction'
    )
)


class UseCasesScene(MovingCameraScene):

    def __init__(self,
                 zoom_each_use_case=True,
                 all_groups_at_same_time=True,
                 **kwargs):
        super().__init__(**kwargs)
        self.zoom_each_use_case = zoom_each_use_case
        self.all_groups_at_same_time = all_groups_at_same_time

    def construct(self):
        super().construct()

        # Group first row and second row
        use_cases = []
        for use_case1, use_case2 in zip(USE_CASES, USE_CASES[5:]):
            use_case2.next_to(use_case1, DOWN, buff=0.5)
            use_cases.append(VGroup(use_case1, use_case2))
        use_cases = VGroup(*use_cases)

        # Make all use cases same width and height
        widest = sorted(use_cases, key=lambda x: x.width)[-1]
        tallest = sorted(chain.from_iterable(use_cases), key=lambda x: x.height)[-1]
        for use_case_top, use_case_bottom in use_cases:
            use_case_top.frame.stretch_to_fit_width(widest.width)
            use_case_top.frame.stretch_to_fit_height(tallest.height)
            use_case_bottom.frame.stretch_to_fit_width(widest.width)
            use_case_bottom.frame.stretch_to_fit_height(tallest.height)

        use_cases.arrange().scale(0.8).center()

        # Human friendly names for use cases

        antibodies = use_cases[0][0]
        peptides_ppi = use_cases[0][1]

        toxins = use_cases[1][0]
        peptides_sol = use_cases[1][1]

        enzyme_stabilisation = use_cases[2][0]
        microbiome = use_cases[2][1]

        codons = use_cases[3][0]
        detoxification = use_cases[3][1]

        clinical = use_cases[4][0]
        mutation = use_cases[4][1]

        # --- Animate!

        # Create the use cases
        self.play(Write(use_cases))
        self.wait(1)

        # Example setting the camera
        old_width = self.camera.frame_width
        if self.zoom_each_use_case:
            for use_case in chain.from_iterable(use_cases):
                self.play(self.camera.frame.animate.move_to(use_case).set(width=use_case.width * 2))
                self.wait(2)

        self.play(self.camera.frame.animate.move_to(ORIGIN).set(width=old_width))
        self.wait(2)

        # Highlight groups with common / synergies
        # Indications to apply to group visually, temporarily, use cases
        # https://docs.manim.community/en/stable/reference/manim.animation.indication.html

        # Similar modalities (some of these are a bit of a stretch to group together)

        if self.all_groups_at_same_time:
            self.play(
                FocusOn(antibodies),
                # Peptides
                Wiggle(peptides_sol), Wiggle(peptides_ppi),
                # DNA + RNA
                Indicate(codons), Indicate(mutation), Indicate(clinical),
                # Proteins
                Circumscribe(toxins.frame), Circumscribe(detoxification.frame),
                # Enzyme
                ApplyWave(enzyme_stabilisation), ApplyWave(microbiome),
                run_time=3
            )
            self.wait(6)

        # Led by divisions

        # Prediction, optimization, preclinical, clinical...


# --- Molecule "datasets"


MOLECULE_COLORS = RED, BLUE, GREEN, YELLOW, TEAL, GOLD, ORANGE, MAROON, PURPLE, PINK


def new_molecule(svs=SVGS.PROTEIN3D, colors=None):
    molecule = SVGMobject(svs)
    if colors is not None:
        try:
            colors = colors()
        except TypeError:
            ...
        try:
            molecule.set_color_by_gradient(*colors)
        except TypeError:
            molecule.set_color_by_gradient(colors)
    return molecule


def molecule_from_molecule(molecule, svs=SVGS.ANTIBODY):
    new_molecule = SVGMobject(svs)
    new_molecule.match_width(molecule)
    new_molecule.match_color(molecule)
    new_molecule.move_to(molecule)
    return new_molecule


def molecule_universe_grid(moleculer=new_molecule, num_molecules: Union[int, Tuple[int, int]] = 5):

    try:
        num_rows, num_cols = num_molecules
    except TypeError:
        num_rows = num_cols = num_molecules

    proteins = []
    prev_row = None
    for _ in range(num_rows):
        protein = moleculer()
        if prev_row is not None:
            protein.next_to(prev_row, RIGHT, buff=SMALL_BUFF)
        prev_row = protein
        proteins.append(protein)
        prev_col = protein
        for _ in range(num_cols):
            protein = moleculer()
            protein.next_to(prev_col, DOWN, buff=SMALL_BUFF)
            prev_col = protein
            proteins.append(protein)

    return VGroup(*proteins)


HIV_ANTIBODY_HEAVY = """
QVQLVQSGAEVKKPGASVKVSCQASGYRFSNFVIHWVRQA
PGQRFEWMGWINPYNGNKEFSAKFQDRVTFTADTSANTAY
MELRSLRSADTAVYYCARVGPYSWDDSPQDNYYMDVWGKG
TTVIVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYF
PEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPS
SSLGTQTYICNVNHKPSNTKVDKKAEPKSCDKTHTCPPCP
APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHED
PEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLH
QDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYT
LPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENN
YKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHE
ALHNHYTQKSLSLSPGK
""".strip()

HIB_ANTIBODY_HCDR3 = 'ARVGPYSWDDSPQDNYYMDV'

HIV_ANTIBODY_LIGHT = """
EIVLTQSPGTLSLSPGERATFSCRSSHSIRSRRVAWYQHK
PGQAPRLVIHGVSNRASGISDRFSGSGSGTDFTLTITRVE
PEDFALYYCQVYGASSYTFGQGTKLERKRTVAAPSVFIFP
PSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNS
QESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQ
GLRSPVTKSFNRGEC
""".strip()


def hiv_antibody():
    """Returns mobjects for the hiv antibody sequence"""

    # Amino acid level - heavy and HCDR3 of the HIV antibody
    # Likely this is too complex for the purpose in this video
    heavy = [
        aa_text(line)
        for line in HIV_ANTIBODY_HEAVY.splitlines(keepends=False)[:-1]
        # N.B. for aesthetics, remove the last line
    ]
    heavy = VGroup(*heavy).scale(0.4).arrange(DOWN, center=False)
    hcdr3 = aa_text(HIB_ANTIBODY_HCDR3).scale(0.4)

    return heavy, hcdr3


MONO_FONTS = (
    'Apercu Mono Pro',
    'Monaco',
)


def aa_text(text, color_dict=None):
    if color_dict is None:
        color_dict = get_aa_hex_color_dict()
    return Text(
        text,
        t2c=color_dict,
        font='Monaco',
        weight=BOLD
    )


def embedding_text(mobject):
    return Text(
        '[0.42, 0.22, 0.33, 0.78, ..., 0.11, 0.92, 0, 0.13]',
        weight=BOLD,
        font='Apercu Mono Pro'
    ).scale(0.4).next_to(mobject, DOWN)


# noinspection PyAbstractClass
class MLM(VGroup):

    def __init__(self,
                 text='PROTEIN',
                 masked_text='PR?TEI?',
                 encoder=None,
                 decoder=None,
                 **kwargs):
        super().__init__(**kwargs)

        self.text = aa_text(text)
        self.masked_text = aa_text(masked_text).move_to(self.text)

        if encoder is None:
            encoder = (
                SVGMobject(SVGS.MODEL)
                .scale(0.5)
                .set_color_by_gradient(RED, GREEN, BLUE)
                .next_to(self.masked_text, DOWN)
            )
        self.encoder = encoder

        self.embedding = embedding_text(encoder)

        if decoder is None:
            decoder = (
                SVGMobject(SVGS.MODEL)
                .scale(0.5)
                .set_color_by_gradient(BLUE, GREEN, RED)
                .next_to(self.embedding, DOWN)
            )
        self.decoder = decoder

        self.masked_decoded = self.masked_text.copy().next_to(decoder, DOWN)
        self.decoded = self.text.copy().next_to(decoder, DOWN)

        self.add(
            self.text,
            self.masked_text,
            self.encoder,
            self.embedding,
            self.decoder,
            self.masked_decoded,
            self.decoded
        )

    def play_encode_decode(self, scene):
        scene.play(FadeIn(self.text))
        scene.wait(1)
        scene.play(ReplacementTransform(self.text, self.masked_text))
        scene.wait(1)
        scene.play(FadeIn(self.encoder), FadeIn(self.decoder))
        scene.wait(1)
        scene.play(
            ReplacementTransform(self.text.copy(), self.encoder),
            ReplacementTransform(self.encoder.copy(), self.embedding),
            ReplacementTransform(self.embedding.copy(), self.decoder),
            ReplacementTransform(self.decoder.copy(), self.masked_decoded)
        )
        scene.wait(1)
        scene.play(ReplacementTransform(self.masked_decoded, self.decoded))
        scene.wait(1)


class PretrainingScene(ZoomedScene):

    def __init__(self, seed=42, **kwargs):
        super().__init__(**kwargs)
        self.seed = seed

    def construct(self):

        rng = np.random.RandomState(self.seed)

        moleculizer = partial(new_molecule, colors=lambda: rng.choice(MOLECULE_COLORS, size=3, replace=True))

        proteins = molecule_universe_grid(moleculizer, num_molecules=(5, 8)).scale(0.3).center()
        proteins_frame = SurroundingRectangle(proteins, color=WHITE).round_corners(0.5)
        universe_title = Text('Known proteins').scale(0.7).next_to(proteins_frame, UP)
        framed_proteins = VGroup(universe_title, proteins_frame, proteins).to_edge(LEFT)

        self.play(FadeIn(framed_proteins))
        self.wait(3)

        # --- Sequence -> Number -> MLM

        heavy, hcdr3 = hiv_antibody()
        heavy.next_to(framed_proteins, RIGHT)
        self.play(Transform(proteins[23], heavy))
        self.wait(5)

        mlm = MLM().next_to(framed_proteins, RIGHT)
        self.play(ShrinkToCenter(heavy))  # , mlm.text
        self.wait(5)

        # FIXME: just to break here
        if mlm.text is not None:
            return

        # --- Subsets
        selected_proteins = [
            proteins[i].copy() for i in
            rng.choice(len(proteins), replace=False, size=len(proteins) // 3)
        ]

        antibodies = VGroup(
            *[molecule_from_molecule(protein) for protein in selected_proteins]
        )
        antibodies_frame = (SurroundingRectangle(antibodies)
                            .match_color(proteins_frame)
                            .match_width(proteins_frame)
                            .match_height(proteins_frame)).round_corners(0.5)
        framed_antibodies = VGroup(antibodies_frame, antibodies)

        framed_antibodies.next_to(proteins, RIGHT, buff=MED_LARGE_BUFF)
        self.play(*[
            ReplacementTransform(protein, antibody) for protein, antibody in zip(selected_proteins, antibodies)],
            Write(antibodies_frame)
        )
        self.wait(4)


class MLMScene(Scene):

    def construct(self):
        super().construct()
        MLM('ANTIBODY', 'A?TIB?DY').play_encode_decode(self)


class FineTuning(Scene):

    # To stress the need for less data in downstream tasks

    def construct(self):
        super().construct()
        raise NotImplemented


if __name__ == '__main__':
    quality = 'l'
    # quality = 'h'
    preview = True
    manimce(
        # LoLLogoScene,
        # LoLCommonsIntroScene,
        # EroomScene,
        # UseCasesScene,
        PretrainingScene,
        # MLMScene,
        quality=quality,
        preview=preview,
        save_last_frame=False,
    )


#
# --- "Company abstractions"
#
# from typing import Union, Iterable
#
# from manim import GREEN, BLUE, RED
#
#
# class Colleague:
#
#     def __init__(self, name):
#         super().__init__()
#         self.name = name
#
#
# class Division:
#
#     def __init__(self, name, short, color):
#         super().__init__()
#         self.name = name
#         self.short = short
#         self.color = color
#
#
# DIVISION_PH = Division(name='Pharma', short='PH', color=BLUE)
# DIVISION_CS = Division(name='Crop Science', short='CS', color=GREEN)
# DIVISION_CH = Division(name='Consumer Health', short='CH', color=RED)
#
#
# class UseCase:
#
#     def __init__(self,
#                  name: str,
#                  divisions: Union[Division, Iterable[Division]],
#                  team):
#         super().__init__()
#         self.name = name
#         self.divisions = divisions
#         self.team = team
#

#
# """Create git repos to visualize with gource."""
#
# from git import Commit, Actor
#
# Actor(name='Olaf')
# Commit()
#
