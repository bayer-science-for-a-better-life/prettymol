from logomaker.src.colors import get_color_dict
from matplotlib.colors import to_hex

from prettymol.config import Config

if Config.prefer_manimgl():
    from manimlib import (Scene,
                          Text, SVGMobject, VMobject,
                          RED, BLUE, GREEN, PURPLE, LEFT, RIGHT, DOWN,
                          Write, ReplacementTransform, Transform, FadeIn)
else:
    from manim import (Scene,
                       Text, SVGMobject, VMobject,
                       RED, BLUE, GREEN, PURPLE, LEFT, RIGHT, DOWN,
                       Write, ReplacementTransform, Transform, FadeIn)


def get_hex_color_dict(name, chars):
    return {char: to_hex(color) for char, color in get_color_dict(name, chars).items()}


XARELTO_SMILES = 'C1COCC(=O)N1C2=CC=C(C=C2)N3CC(OC3=O)CNC(=O)C4=CC=C(S4)Cl'
ADENINE_SMILES = 'C1=NC2=NC=NC(=C2N1)N'
DNA_SVG = str(Config.DEFAULT_IMAGES_PATH / 'dna.svg')


class LoLCommonsIntro(Scene):

    def construct(self):

        # Atom level
        adenine = Text(ADENINE_SMILES)

        # Nucleotide level
        dna = Text('ACTGAATATAGACTATA', t2c={'A': RED, 'C': BLUE, 'T': GREEN, 'G': PURPLE})

        # Amino acid level
        aas = 'VQGGAAVQQEVLA'
        aas = Text(aas, t2c=get_hex_color_dict('skylign_protein', aas))

        # Logo level
        lol = Text('Language of Life')
        lol[1].set_color(BLUE)
        lol[13].set_color(GREEN)
        # Beware: currently semantics seem to have changed in manim ce and this should be 11
        #         are they not counting spaces?

        commons = Text('commons')
        commons.next_to(lol, DOWN)

        double_helix_left = (SVGMobject(DNA_SVG).
                             scale(0.4).
                             set_color_by_gradient(GREEN, BLUE).
                             set_stroke(width=0.5))
        double_helix_left.next_to(commons, LEFT)

        double_helix_right = double_helix_left.copy().set_color_by_gradient(BLUE, GREEN)
        double_helix_right.next_to(commons, RIGHT)

        # Animate... from nucleotide to LoL Commons
        self.play(Write(adenine))
        self.wait()
        self.play(ReplacementTransform(adenine, dna))
        self.wait()
        self.play(ReplacementTransform(dna, aas))
        self.wait()
        self.play(ReplacementTransform(aas, lol))
        self.wait()
        self.play(Write(commons), FadeIn(double_helix_left), FadeIn(double_helix_right))
        self.wait()

        # Move logo to left lower corner
        logo = VMobject()
        logo.add(lol, double_helix_left, double_helix_right, commons)
        mini_logo = logo.copy()
        mini_logo.scale(0.5)
        mini_logo.to_corner()
        self.play(Transform(logo, mini_logo))

        # And keep the scene
        self.play(Write(Text('Corporate info follows...')))
        self.wait(5)

        # FIXME: also color smiles (see compounds scenes)
        # FIXME: also color nucleotides
        # FIXME: also make dna and protein sequences meaningful
        # FIXME: does it make sense to transform adenine only in these bases?
        # FIXME: does it make sense to transform codons into appropriate amino acids?


if __name__ == '__main__':
    if Config.prefer_manimgl():
        from prettymol.manim_utils import manimgl
        manimgl(LoLCommonsIntro, write=True)
    else:
        from prettymol.manim_utils import manimce
        manimce(LoLCommonsIntro)
    # TODO: convenience to route depending on the actual type of the scene in manim_utils
