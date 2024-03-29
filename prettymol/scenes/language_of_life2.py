from manimlib import *

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

IMAGES_PATH = Config.DEFAULT_IMAGES_PATH

# TODO: likely we should provide a "media library"

DNA_SVG = str(IMAGES_PATH / 'dna.svg')
DNASTRUCT_SVG = str(IMAGES_PATH / 'DNA_pymol.svg')
RNASTRUCT_SVG = str(IMAGES_PATH / 'RNA.svg')
PEPSTRUCT_SVG = str(IMAGES_PATH / 'Peptide.svg')
BCIRCLE_SVG = str(IMAGES_PATH / 'CircleLogo.svg')
MACHINE_SVG = str(IMAGES_PATH / 'machine3.svg')
PILLS_SVG = str(IMAGES_PATH / 'Pills.svg')
CORN_SVG = str(IMAGES_PATH / 'Corn.svg')
TOXIN_PNG = str(IMAGES_PATH / 'toxin.png')
MACHINE_PNG = str(IMAGES_PATH / 'machine.png')
CROSS_PNG = str(IMAGES_PATH / 'Corp-Logo_BG_Bayer-Cross_Basic_150dpi_on-screen_RGB.png')
AI_PNG = str(IMAGES_PATH / 'Ai-image-deep-net-v3-690px.jpg')
VIP_PNG = str(IMAGES_PATH / 'VIP2.png')
CRY_PNG = str(IMAGES_PATH / 'Cry2.png')

WCR = str(IMAGES_PATH / 'WCR.jpg')
STINK = str(IMAGES_PATH / 'StinkBug.jpg')
FAW = str(IMAGES_PATH / 'FAW.jpg')

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
        cornstr = (SVGMobject(CORN_SVG).
                  scale(1).
                  set_color_by_gradient(GREEN, BLUE).
                  set_stroke(width=0.5))
        pillsstr = (SVGMobject(PILLS_SVG).
                  scale(1).
                  set_color_by_gradient(GREEN, BLUE).
                  set_stroke(width=0.5))
        
        self.play(FadeIn(machine))
        self.play(FadeOut(grid))
        self.play(machine.scale, 4.0)
        
        
        cornstr_left=cornstr.next_to(machine, LEFT)
        cornstr_left_up=cornstr_left.shift(UP)
        pillsstr.next_to(machine, RIGHT)
        self.play(FadeIn(cornstr_left_up))
        self.play(FadeIn(pillsstr))
       

class CS_Use(Scene):
    def construct(self):

        _, _ = lol_commons_minilogo(self)

        source = Text("imagine...", height=1)
        source[2].set_color(BLUE)
        source[4].set_color(GREEN)
        target = Text("a toolbox that...", height=0.5)
        target[2:10].set_color(BLUE)
        target2a = Text("turns sequence", height=0.5)
        target2a[6:14].set_color(GREEN)

        aas = Text("VQGGAAVQQEVLA", height=0.5)
        aas.next_to(target2a, DOWN)

        target2b= Text("and structure", height=0.5) 
        target2b[4:13].set_color(BLUE)
        ### Protein structure
        
        target2c= Text("into a form for", height=0.5) 
        target2d= Text("machine learning", height=0.5) 
        target2d[1].set_color(BLUE)
        target2d[4].set_color(GREEN)
        target2d[10].set_color(BLUE)
        target2d[13].set_color(GREEN)
        
        target2 = Text("to discover", height=0.5)
        target2[3:15].set_color(GREEN)
        target3 = Text("and design", height=0.5)
        target3[4:10].set_color(GREEN)
        target4 = Text("new insect control", height=0.5)    
        target5 = Text("toxins", height=1)
        target6 = Text("to combat insect pests", height=0.5)
        
        corona= ImageMobject(VIP_PNG)
        corona.scale(1.2)
        
        ai = ImageMobject(AI_PNG)
        ai.scale(1.2)
        


        
        self.play(Write(source))
        self.wait()
        kw = {"run_time": 3, "path_arc": PI / 2}

        self.play(TransformMatchingShapes(source, target, **kw))
        self.wait()
        self.play(TransformMatchingShapes(target, target2a, **kw))
        self.wait()
        self.add(aas)
        self.wait()
        self.play(FadeOut(aas))
        self.play(TransformMatchingShapes(target2a, target2b, **kw))
        self.wait()
        self.add(corona)
        self.bring_to_back(corona)
        self.wait()
        self.play(FadeOut(corona))
        self.wait()
        self.play(TransformMatchingShapes(target2b, target2c, **kw))
        self.wait()

        self.play(TransformMatchingShapes(target2c, target2d, **kw))
        self.wait()
        self.play(FadeOut(target2d))

        grid = Tex("01").get_grid(10, 10, height=7)

        self.add(aas)
        self.play(aas.move_to, UP*3, run_time=2)
        self.add(corona)
        self.wait()
        self.play(FadeOut(corona))
        self.play(TransformMatchingShapes(aas, grid, **kw))
        self.play(grid.animate.apply_complex_function(np.exp), run_time=5)

        self.play(FadeOut(grid))
        self.play(FadeIn(ai))
        self.wait()
        self.play(FadeOut(ai)) 
               
        self.play(TransformMatchingShapes(target2, target3, **kw))
        self.wait()
        self.play(TransformMatchingShapes(target3, target4, **kw))
        self.wait()
        self.play(TransformMatchingShapes(target4, target5, **kw))
        self.play(TransformMatchingShapes(target5, target6, **kw))
        self.play(target6.move_to, UP*3, run_time=2)
        faw = ImageMobject(FAW)
        self.play(faw.scale, 0.6)
        self.add(faw)
        self.play(faw.move_to, LEFT*3, run_time=2)
        stink = ImageMobject(STINK)
        self.play(stink.scale, 0.6)
        self.add(stink)
        self.play(stink.move_to, RIGHT*3, run_time=2)
        wcr = ImageMobject(WCR)
        self.play(wcr.scale, 0.6)
        self.add(wcr)
        
        
        
        
        
class Info4(Scene):
    def construct(self):

        _, _ = lol_commons_minilogo(self)
        
        corona= ImageMobject(TOXIN_PNG)
        corona.scale(1.2)
        corona.to_edge(RIGHT, buff=1)
        self.add(corona)        
        
if __name__ == '__main__': 
    from prettymol.manim_utils import manimgl
    manimgl(Info1, Info2, Info3, Info4, write=True)
