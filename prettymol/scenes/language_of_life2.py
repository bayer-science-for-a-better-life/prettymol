from manimlib import *
import numpy as np
from pathlib import Path

DNA_SVG = "..\..\media\images\dna.svg"
#TOOLS = "toolbox.png"
DNASTRUCT_SVG = "..\..\media\images\DNA_pymol.svg"
RNASTRUCT_SVG = "..\..\media\images\RNA.svg"
PEPSTRUCT_SVG = "..\..\media\images\Peptide.svg"

BLOGO_SVG = "..\..\media\images\Logo.svg"


class TexTransformExample(Scene):
    def construct(self):

        # Alternatively, if you don't want to think about breaking up
        # the tex strings deliberately, you can TransformMatchingShapes,
        # which will try to line up all pieces of a source mobject with
        # those of a target, regardless of the submobject hierarchy in
        # each one, according to whether those pieces have the same
        # shape (as best it can).
        
        
        lol = Text('Language of Life')
        
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
        self.play(Transform(mini_logo, mini_logo))
        
        
        source = Text("imagine...", height=1)
        target = Text("ai", height=1)

        self.play(Write(source))
        self.wait()
        kw = {"run_time": 3, "path_arc": PI / 2}
        self.play(TransformMatchingShapes(source, target, **kw))
        self.wait()

        rnastr= (SVGMobject(RNASTRUCT_SVG).
                             scale(3).
                             set_color_by_gradient(GREEN, BLUE).
                             set_stroke(width=0.5))
        
        self.play(ReplacementTransform(target, rnastr))
        
        dnastr= (SVGMobject(DNASTRUCT_SVG).
                             scale(3).
                             set_color_by_gradient(GREEN, BLUE).
                             set_stroke(width=0.5))
        
        self.play(ReplacementTransform(rnastr, dnastr))
        
        pepstr= (SVGMobject(PEPSTRUCT_SVG).
                             scale(2).
                             set_color_by_gradient(GREEN, BLUE).
                             set_stroke(width=0.1))
        
        self.play(ReplacementTransform(dnastr,pepstr))
        
        blogo = (SVGMobject(BLOGO_SVG).
                             scale(5).
                             set_color_by_gradient(GREEN, BLUE).
                             set_stroke(width=0.1))
        
        self.play(ReplacementTransform(pepstr,blogo))