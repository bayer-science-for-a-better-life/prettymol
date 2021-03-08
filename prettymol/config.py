import os
from enum import Enum


class ManimImplementation(Enum):
    manimgl = 'manimgl'
    manimce = 'manimce'


class Config:

    @staticmethod
    def prefer_manimgl() -> bool:
        preferred_manim = os.environ.get('PRETTYMOL_PREFERRED_MANIM', ManimImplementation.manimgl.name).lower()
        possible_manims = [item.value for item in ManimImplementation]
        if preferred_manim not in possible_manims:
            raise ValueError(f'Unknown manim {preferred_manim}, must be one of {possible_manims}')
        return preferred_manim == ManimImplementation.manimgl.value
