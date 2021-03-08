from manimlib import *

from prettymol.config import Config

IMAGES_PATH = Config.DEFAULT_IMAGES_PATH

# TODO: likely we should provide a "media library"

DNA_SVG = str(IMAGES_PATH / 'dna.svg')
DNASTRUCT_SVG = str(IMAGES_PATH / 'DNA_pymol.svg')
RNASTRUCT_SVG = str(IMAGES_PATH / 'RNA.svg')
PEPSTRUCT_SVG = str(IMAGES_PATH / 'Peptide.svg')
BCIRCLE_SVG = str(IMAGES_PATH / 'CircleLogo.svg')
MACHINE_SVG = str(IMAGES_PATH / 'machine3.svg')
BLOGO_SVG = str(IMAGES_PATH / 'Logo.svg')
dnastr = (SVGMobject(DNASTRUCT_SVG).
          scale(3).
          set_color_by_gradient(GREEN, BLUE).
          set_stroke(width=0.5))


def lol_commons_minilogo(scene):
    # TODO: make a little library of logos, use everywhere
    lol = Text('Language of Life')
    lol[1].set_color(BLUE)
    lol[13].set_color(GREEN)

    commons = Text('commons')
    commons.next_to(lol, DOWN)
    double_helix_left = (SVGMobject(DNA_SVG).
                         scale(0.4).
                         set_color_by_gradient(GREEN, BLUE).
                         set_stroke(width=0.5))
    double_helix_left.next_to(commons, LEFT)
    double_helix_right = double_helix_left.copy().set_color_by_gradient(BLUE, GREEN)
    double_helix_right.next_to(commons, RIGHT)

    logo = VMobject()
    logo.add(lol, double_helix_left, double_helix_right, commons)
    mini_logo = logo.copy()
    mini_logo.scale(0.5)
    mini_logo.to_corner()
    scene.play(Transform(mini_logo, mini_logo))

    return logo, mini_logo


class Info1(Scene):

    def construct(self):

        _, _ = lol_commons_minilogo(self)

        source = Text("imagine...", height=1)
        source[2].set_color(BLUE)
        source[4].set_color(GREEN)
        target = Text("a toolbox..", height=1)
        target2 = Text("for the...", height=1)
        target3 = Text("language of life", height=1)
        target3[1].set_color(BLUE)
        target3[13].set_color(GREEN)

        self.play(Write(source))
        self.wait()
        kw = {"run_time": 3, "path_arc": PI / 2}
        self.play(TransformMatchingShapes(source, target, **kw))
        self.wait()
        self.play(TransformMatchingShapes(target, target2, **kw))
        self.play(TransformMatchingShapes(target2, target3, **kw))
        self.play(FadeOut(target3))


class Info2(Scene):
    def construct(self):
        # Alternatively, if you don't want to think about breaking up
        # the tex strings deliberately, you can TransformMatchingShapes,
        # which will try to line up all pieces of a source mobject with
        # those of a target, regardless of the submobject hierarchy in
        # each one, according to whether those pieces have the same
        # shape (as best it can).

        _, _ = lol_commons_minilogo(self)

        rnatext = Text("RNA", font="Consolas", font_size=90)
        rnastr = (SVGMobject(RNASTRUCT_SVG).
                  scale(3).
                  set_color_by_gradient(GREEN, BLUE).
                  set_stroke(width=0.5))
        VGroup(rnatext, rnastr).arrange(DOWN, buff=1)
        self.play(Write(rnatext))
        self.play(FadeIn(rnastr))
        self.play(FadeOut(rnastr))
        self.play(FadeOut(rnatext))

        dnatext = Text("DNA", font="Consolas", font_size=90)

        VGroup(dnatext, dnastr).arrange(DOWN, buff=1)
        # self.play(ReplacementTransform(rnastr, dnastr))
        self.play(FadeIn(dnastr))
        self.play(Write(dnatext))
        self.play(FadeOut(dnastr))
        self.play(FadeOut(dnatext))

        peptext = Text("Peptides", font="Consolas", font_size=90)
        pepstr = (SVGMobject(PEPSTRUCT_SVG).
                  scale(2).
                  set_color_by_gradient(GREEN, BLUE).
                  set_stroke(width=0.1))
        VGroup(peptext, pepstr).arrange(DOWN, buff=1)
        # self.play(ReplacementTransform(dnastr,pepstr))
        self.play(FadeIn(pepstr))
        self.play(Write(peptext))
        self.play(FadeOut(peptext))
        self.play(FadeOut(pepstr))

        # TODO: Peptide doesn't look right - make new figure

        prottext = Text("Proteins", font="Consolas", font_size=90)
        #        VGroup(peptext, pepstr).arrange(DOWN, buff=1)
        # self.play(ReplacementTransform(dnastr,pepstr))
        #  self.play(FadeIn(pepstr))
        self.play(Write(prottext))
        self.play(FadeOut(prottext))


# ---  Add in protein picture


class Info3(Scene):
    def construct(self):

        _, _ = lol_commons_minilogo(self)

        self.play(FadeIn(dnastr))
        grid = Tex("01").get_grid(10, 10, height=7)
        # self.add(grid)
        self.play(ReplacementTransform(dnastr, grid))

        self.play(grid.animate.apply_complex_function(np.exp), run_time=5)
        self.wait()

        machine = (SVGMobject(MACHINE_SVG).
                   scale(0.5).
                   set_color_by_gradient(GREEN, BLUE).
                   set_stroke(width=0.5))
        self.play(FadeIn(machine))
        self.play(machine.scale, 6.0)
        self.play(FadeOut(grid))

        # TODO: Add in products for Crop Science and Pharam (maybe call them disruptive?)


if __name__ == '__main__':
    from prettymol.manim_utils import manimgl
    manimgl(Info1, Info2, Info3, write=True)
