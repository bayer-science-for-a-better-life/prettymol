import os
from enum import Enum
from pathlib import Path


class ManimImplementation(Enum):
    manimgl = 'manimgl'
    manimce = 'manimce'


class Config:

    DEFAULT_MEDIA_PATH = Path(__file__).parent.parent / 'media'
    DEFAULT_VIDEOS_PATH = DEFAULT_MEDIA_PATH / 'videos'
    DEFAULT_RENDERED_SCENES_PATH = DEFAULT_VIDEOS_PATH
    DEFAULT_IMAGES_PATH = DEFAULT_MEDIA_PATH / 'images'

    @staticmethod
    def prefer_manimgl(preferred_manim='manimgl') -> bool:
        if preferred_manim is None:
            preferred_manim = os.environ.get('PRETTYMOL_PREFERRED_MANIM', ManimImplementation.manimgl.name).lower()
        possible_manims = [item.value for item in ManimImplementation]
        if preferred_manim not in possible_manims:
            raise ValueError(f'Unknown manim {preferred_manim}, must be one of {possible_manims}')
        return preferred_manim == ManimImplementation.manimgl.value
