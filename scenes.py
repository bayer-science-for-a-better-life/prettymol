import numpy as np
from manimlib.imports import Scene, Dot, Circle, Square, Polygon, ShowCreation, FadeOut, GrowFromCenter, Transform


ARTEMISININ_SMILES = 'O=C3O[C@@H]4O[C@@]1(OO[C@@]42[C@@H](CC1)[C@H](C)CC[C@H]2[C@H]3C)C'


class ExampleScene(Scene):
    def construct(self):
        dot = Dot()
        self.add(dot)
        self.wait(2)


class Shapes(Scene):

    def construct(self):
        circle = Circle()
        square = Square()
        triangle = Polygon(np.array([0, 0, 0]), np.array([1, 1, 0]), np.array([1, -1, 0]))

        self.play(ShowCreation(circle, lag_ratio=1000))
        self.play(FadeOut(circle, lag_ratio=200))
        self.play(GrowFromCenter(square))
        self.play(Transform(square, triangle))


if __name__ == "__main__":
    import os
    module_name = os.path.basename(__file__)
    command_A = "manim -pl --video_dir ~/Downloads/  "
    command_B = module_name + " " + "Shapes"
    os.system(command_A + command_B)
