from manimlib import *
import numpy as np
from pathlib import Path

DNA_SVG = "..\..\media\images\dna.svg"
#TOOLS = "toolbox.png"
DNASTRUCT_SVG = "..\..\media\images\DNA_pymol.svg"
RNASTRUCT_SVG = "..\..\media\images\RNA.svg"
PEPSTRUCT_SVG = "..\..\media\images\Peptide.svg"
BCIRCLE_SVG = "..\..\media\images\CircleLogo.svg"

BLOGO_SVG = "..\..\media\images\Logo.svg"
dnastr= (SVGMobject(DNASTRUCT_SVG).
                             scale(3).
                             set_color_by_gradient(GREEN, BLUE).
                             set_stroke(width=0.5))

class Info1(Scene):
    def construct(self):

        # Alternatively, if you don't want to think about breaking up
        # the tex strings deliberately, you can TransformMatchingShapes,
        # which will try to line up all pieces of a source mobject with
        # those of a target, regardless of the submobject hierarchy in
        # each one, according to whether those pieces have the same
        # shape (as best it can).
        
        
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
        self.play(Transform(mini_logo, mini_logo))
        
        
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
        self.play(Transform(mini_logo, mini_logo))
        
        rnatext = Text("RNA", font="Consolas", font_size=90)      
        rnastr= (SVGMobject(RNASTRUCT_SVG).
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
        #self.play(ReplacementTransform(rnastr, dnastr))
        self.play(FadeIn(dnastr))
        self.play(Write(dnatext))
        self.play(FadeOut(dnastr))
        self.play(FadeOut(dnatext))
        
        peptext = Text("Peptides", font="Consolas", font_size=90)
        pepstr= (SVGMobject(PEPSTRUCT_SVG).
                             scale(2).
                             set_color_by_gradient(GREEN, BLUE).
                             set_stroke(width=0.1))
        VGroup(peptext, pepstr).arrange(DOWN, buff=1)
        #self.play(ReplacementTransform(dnastr,pepstr))
        self.play(FadeIn(pepstr))
        self.play(Write(peptext)) 
        self.play(FadeOut(peptext))
        self.play(FadeOut(pepstr))
        
        prottext = Text("Proteins", font="Consolas", font_size=90)
#        VGroup(peptext, pepstr).arrange(DOWN, buff=1)
        #self.play(ReplacementTransform(dnastr,pepstr))
      #  self.play(FadeIn(pepstr))
        self.play(Write(prottext)) 
        self.play(FadeOut(prottext))
    #    self.play(FadeOut(pepstr))
    
    
class Info3(Scene):
    def construct(self):

        # Alternatively, if you don't want to think about breaking up
        # the tex strings deliberately, you can TransformMatchingShapes,
        # which will try to line up all pieces of a source mobject with
        # those of a target, regardless of the submobject hierarchy in
        # each one, according to whether those pieces have the same
        # shape (as best it can).
        
        
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
        self.play(Transform(mini_logo, mini_logo))
        
        self.play(FadeIn(dnastr))
        grid = Tex("01").get_grid(10, 10, height=4)
        #self.add(grid)
        self.play(ReplacementTransform(dnastr,grid))
        # You can animate the application of mobject methods with the
        # ".animate" syntax:
        #self.play(grid.animate.shift(LEFT))

 