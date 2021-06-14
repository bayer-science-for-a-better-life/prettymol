from itertools import chain
from pathlib import Path
from typing import Union, Optional

from logomaker.src.colors import get_color_dict
from manim import (Scene,
                   Text, SVGMobject, RED, BLUE, GREEN, PURPLE, LEFT, RIGHT, DOWN,
                   Write, ReplacementTransform, Transform, FadeIn, VGroup, BarChart, UP, Mobject, VMobject)
from matplotlib.colors import to_hex

from prettymol.config import Config
from prettymol.manim_utils import manimce, is_manimce


# --- Utils

def get_hex_color_dict(name, chars):
    return {char: to_hex(color) for char, color in get_color_dict(name, chars).items()}


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


# --- MObjects

# noinspection PyAbstractClass
class LoLLogo(VGroup):

    def __init__(self, commons_left=SVGS.DNA, commons_right=None, **kwargs):
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
        super().__init__(self.lol, self.commons, **kwargs)

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

    def replace_left_right(self, scene, commons_left, commons_right=None, run_time_s=1):
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
                 **kwargs):
        super().__init__(**kwargs)
        self.replacers = replacers
        self.add_vector_transform = add_vector_transform

    def construct(self):

        # Atom level
        adenine = Text(ADENINE_SMILES)

        # Nucleotide level
        dna = Text('ACTGAATATAGACTATA', t2c={'A': RED, 'C': BLUE, 'T': GREEN, 'G': PURPLE})

        # Amino acid level
        aas = 'VQGGAAVQQEVLA'
        aas = Text(aas, t2c=get_hex_color_dict('skylign_protein', aas))

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
            logo.replace_left_right(self, replacer)
            self.wait(1)

        # Move logo to left lower corner
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

        a_million_forks = Text('A million forks in the road to drug development')
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
        ).scale(0.8).to_edge(LEFT).shift(DOWN)

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
                SVGMobject(SVGS.CYCLIC_PEPTIDE).set_color(BLUE, GREEN).scale(0.8),
                Text('Can we generate peptides without liabilities?')
            ),
            (
                SVGMobject(SVGS.DNA).set_color(GREEN).scale(0.8),
                Text('Can we optimize DNA to better express plant traits?')
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
                SVGMobject(SVGS.PATIENT).set_color(BLUE).scale(0.8),
                Text('Will a patient respond to treatment in a clinical trial?')
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
        self.play(ReplacementTransform(questions, logo))
        self.wait(2)


if __name__ == '__main__':

    # quality = 'l'
    quality = 'h'
    preview = True
    manimce(
        # LoLLogoScene,
        # LoLCommonsIntroScene,
        EroomScene,
        quality=quality,
        preview=preview,
        save_last_frame=False,
    )
