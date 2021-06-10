"""Helpers to deal with manim."""
import inspect
import sys
from contextlib import contextmanager
from itertools import groupby
from pathlib import Path

from typing import Union, Optional

from prettymol.config import Config


# --- utils

@contextmanager
def sys_argv(program_name: str, *args: str):
    """
    A context manager that temporarily changes the value of sys.argv.
    """
    # We should refactor manim "main" this to accept args
    #   https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
    # In the meantime, we trick it via sys.argv hijacking so we can debug with no problem...
    old_argv = sys.argv
    sys.argv = [program_name] + list(args)
    yield
    sys.argv = old_argv


# --- manimgl

def _patch_osx_window_screeninfo_get_monitors():

    # This is so that screeninfo does not choke in OSX when XRandr is also installed
    # Should instead open a PR for upstream

    def get_monitors_osx():
        """Returns a list of :class:`Monitor` objects based on active monitors."""
        # noinspection PyUnresolvedReferences,PyProtectedMember
        from screeninfo.screeninfo import Enumerator, _get_monitors, ScreenInfoError
        # noinspection TryExceptPass,PyBroadException
        try:
            return _get_monitors(Enumerator.OSX)
        except Exception:
            pass
        raise ScreenInfoError("No enumerators available")

    if sys.platform == 'darwin':
        import manimlib
        manimlib.get_monitors = get_monitors_osx
        manimlib.window.get_monitors = get_monitors_osx
        # noinspection PyUnresolvedReferences
        manimlib.config.get_monitors = get_monitors_osx


def manimgl(*scenes,
            write_file=False,
            skip_animations=False,
            quality=None,
            full_screen=False,
            save_pngs=False,
            gif=False,
            transparent=False,
            quiet=False,
            write_all=False,
            open_file=False,
            finder=False,
            config=None,
            file_name=None,
            start_at_animation_number=None,
            resolution=None,
            frame_rate=None,
            color=None,
            leave_progress_bars=False,
            video_dir: Optional[Union[str, Path]] = None,
            show_help=False):
    """
    Python friendly manimgl caller.

    This currently has been tested with manimgl version 1.0.

    From manimgl --help
    -----------------------------------
    usage: manimgl [-h] [-w] [-s] [-l] [-m] [--hd] [--uhd] [-f] [-g] [-i] [-t] [-q] [-a] [-o] [--finder] [--config]
                   [--file_name FILE_NAME] [-n START_AT_ANIMATION_NUMBER] [-r RESOLUTION]
                   [--frame_rate FRAME_RATE] [-c COLOR] [--leave_progress_bars] [--video_dir VIDEO_DIR]
                   [file] [scene_names [scene_names ...]]

    positional arguments:
      file                  path to file holding the python code for the scene
      scene_names           Name of the Scene class you want to see

    optional arguments:
      -h, --help            show this help message and exit
      -w, --write_file      Render the scene as a movie file
      -s, --skip_animations
                            Save the last frame
      -l, --low_quality     Render at a low quality (for faster rendering)
      -m, --medium_quality  Render at a medium quality
      --hd                  Render at a 1080p
      --uhd                 Render at a 4k
      -f, --full_screen     Show window in full screen
      -g, --save_pngs       Save each frame as a png
      -i, --gif             Save the video as gif
      -t, --transparent     Render to a movie file with an alpha channel
      -q, --quiet
      -a, --write_all       Write all the scenes from a file
      -o, --open            Automatically open the saved file once its done
      --finder              Show the output file in finder
      --config              Guide for automatic configuration
      --file_name FILE_NAME
                            Name for the movie or image file
      -n START_AT_ANIMATION_NUMBER, --start_at_animation_number START_AT_ANIMATION_NUMBER
                            Start rendering not from the first animation, butfrom another, specified by its index.
                            If you passin two comma separated values, e.g. "3,6", it will end
                            the rendering at the second value
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution, passed as "WxH", e.g. "1920x1080"
      --frame_rate FRAME_RATE
                            Frame rate, as an integer
      -c COLOR, --color COLOR
                            Background color
      --leave_progress_bars
                            Leave progress bars displayed in terminal
      --video_dir VIDEO_DIR
                            directory to write video
    -----------------------------------
    """

    _patch_osx_window_screeninfo_get_monitors()

    from manimlib.__main__ import main as manimgl_main

    arguments = []

    # --- Options

    if write_file:
        arguments += ['--write_file']

    if skip_animations:
        arguments += ['--skip_animations']

    if quality is not None:
        if quality in ('l', 'low'):
            arguments += ['--low_quality']
        elif quality in ('m', 'medium'):
            arguments += ['--medium_quality']
        elif quality in ('h', 'hd', 'high'):
            arguments += ['--hd']
        elif quality in ('u', 'uhd', 'ultra'):
            arguments += ['--uhd']
        else:
            raise ValueError(f'quality must be one of ["l", "m", "h", "u"]')

    if full_screen:
        arguments += ['--full_screen']

    if save_pngs:
        arguments += ['--save_pngs']

    if gif:
        arguments += ['--gif']

    if transparent:
        arguments += ['--transparent']

    if quiet:
        arguments += ['--quiet']

    if write_all:
        arguments += ['--write_all']

    if open_file:
        arguments += ['--open']

    if finder:
        arguments += ['--finder']

    if config:
        arguments += ['--config']

    if file_name is not None:
        arguments += ['--file_name', str(file_name)]

    if start_at_animation_number is not None:
        arguments += ['--start_at_animation_number', str(start_at_animation_number)]

    if resolution is not None:
        arguments += ['--resolution', str(resolution)]

    if frame_rate is not None:
        arguments += ['--frame_rate', str(frame_rate)]

    if color is not None:
        arguments += ['--color', str(color)]

    if leave_progress_bars:
        arguments += ['--leave_progress_bars']

    if video_dir is None:
        video_dir = Config.DEFAULT_MEDIA_PATH
    arguments += ['--video_dir', str(video_dir)]

    if show_help:
        arguments += ['--help']

    # --- Run!

    path2scenes = sorted((py_file_path, sorted(set(scene.__name__ for scene in scenes_in_file)))
                         for py_file_path, scenes_in_file in
                         groupby(sorted(scenes, key=inspect.getfile), inspect.getfile))

    for py_file_path, scene_names in path2scenes:
        with sys_argv('manimgl', *(arguments + [py_file_path] + scene_names)):
            manimgl_main()


# --- manim.community

def manimce(*scenes,
            # Global options
            config_file=None,
            custom_folders=False,
            disable_caching=False,
            flush_cache=False,
            tex_template=None,
            verbosity='INFO',
            notify_outdated_version=False,
            # Output options
            output_file=None,
            write_to_movie=False,
            media_dir = None,
            log_dir=None,
            log_to_file=False,
            # Render options
            from_animation_number=None,
            write_all=False,
            output_format='mp4',
            save_last_frame=False,
            quality=None,
            resolution=None,
            frame_rate=None,
            renderer='cairo',
            webgl_renderer_path=None,
            transparent=False,
            # Ease of access options
            progress_bar='display',
            preview=False,
            show_in_file_browser=False,
            jupyter=False,
            # Other options
            show_help=False):
    """
    Python friendly manim community caller.

    From manim render --help
    -----------------------------------
    Manim Community v0.7.0

    Usage: manim render [OPTIONS] FILE [SCENE_NAMES]...

      Render SCENE(S) from the input FILE.

      FILE is the file path of the script.

      SCENES is an optional list of scenes in the file.

    Global options:
      -c, --config_file TEXT          Specify the configuration file to use for
                                      render settings.

      --custom_folders                Use the folders defined in the
                                      [custom_folders] section of the config file
                                      to define the output folder structure.

      --disable_caching               Disable the use of the cache (still
                                      generates cache files).

      --flush_cache                   Remove cached partial movie files.
      --tex_template TEXT             Specify a custom TeX template file.
      -v, --verbosity [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                      Verbosity of CLI output. Changes ffmpeg log
                                      level unless 5+.

      --notify_outdated_version / --silent
                                      Display warnings for outdated installation.

    Output options:
      -o, --output_file TEXT          Specify the filename(s) of the rendered
                                      scene(s).

      --write_to_movie                Write to a file.
      --media_dir PATH                Path to store rendered videos and latex.
      --log_dir PATH                  Path to store render logs.
      --log_to_file                   Log terminal output to file.

    Render Options:
      -n, --from_animation_number TEXT
                                      Start rendering from n_0 until n_1. If n_1
                                      is left unspecified, renders all scenes
                                      after n_0.

      -a, --write_all                 Render all scenes in the input file.
      --format [png|gif|mp4]
      -s, --save_last_frame
      -q, --quality [l|m|h|p|k]       Render quality at the follow resolution
                                      framerates, respectively: 854x480 30FPS,
                                      1280x720 30FPS, 1920x1080 60FPS, 2560x1440
                                      60FPS, 3840x2160 60FPS

      -r, --resolution TEXT           Resolution in (W,H) for when 16:9 aspect
                                      ratio isn't possible.

      --fps, --frame_rate FLOAT       Render at this frame rate.
      --renderer [cairo|opengl|webgl]
                                      Select a renderer for your Scene.
      --use_opengl_renderer           Render scenes using OpenGL (Deprecated).
      --use_webgl_renderer            Render scenes using the WebGL frontend
                                      (Deprecated).

      --webgl_renderer_path PATH      The path to the WebGL frontend.
      -g, --save_pngs                 Save each frame as png (Deprecated).
      -i, --save_as_gif               Save as a gif (Deprecated).
      -s, --save_last_frame           Save last frame as png (Deprecated).
      -t, --transparent               Render scenes with alpha channel.

    Ease of access options:
      --progress_bar [display|leave|none]
                                      Display progress bars and/or keep them
                                      displayed.  [default: display]

      -p, --preview                   Preview the Scene's animation. OpenGL does a
                                      live preview in a popup window. Cairo opens
                                      the rendered video file in the system
                                      default media player.

      -f, --show_in_file_browser      Show the output file in the file browser.
      --jupyter                       Using jupyter notebook magic.

    Other options:
      --help                          Show this message and exit.

      Made with <3 by Manim Community developers.
    -----------------------------------
    """

    from manim.__main__ import main as manimce_main

    arguments = ['render']

    # --- Global options

    if config_file is not None:
        arguments += ['--config_file', str(config_file)]

    if custom_folders:
        arguments += ['--custom_folders']

    if disable_caching:
        arguments += ['--disable_caching']

    if flush_cache:
        arguments += ['--flush_cache']

    if tex_template is not None:
        arguments += ['--tex_template', str(tex_template)]

    if verbosity is not None:
        arguments += ['--verbosity', verbosity]

    if notify_outdated_version:
        arguments += ['--notify_outdated_version']

    # --- Output options

    if output_file is not None:
        arguments += ['--output_file', str(output_file)]

    if write_to_movie:
        arguments += ['write_to_movie']

    if media_dir is None:
        media_dir = Config.DEFAULT_MEDIA_PATH
    arguments += ['--media_dir', str(media_dir)]

    if log_dir is not None:
        arguments += ['--log_dir', str(log_dir)]

    if log_to_file:
        arguments += ['--log_to_file']

    # --- Render options

    if from_animation_number is not None:
        arguments += ['--from_animation_number', str(from_animation_number)]

    if write_all:
        arguments += ['--write_all']

    if output_format is not None:
        arguments += ['--format', str(output_format)]

    if save_last_frame:
        arguments += ['--save_last_frame']

    if quality is not None:
        arguments += ['--quality', str(quality)]

    if resolution is not None:
        arguments += ['--resolution', str(resolution)]

    if frame_rate is not None:
        arguments += ['--frame_rate', str(frame_rate)]

    if renderer is not None:
        arguments += ['--renderer', str(renderer)]

    if webgl_renderer_path is not None:
        arguments += ['--webgl_renderer_path', str(webgl_renderer_path)]

    if transparent:
        arguments += ['--transparent']

    # --- Ease of access options

    if progress_bar is not None:
        arguments += ['--progress_bar', str(progress_bar)]

    if preview:
        arguments += ['--preview']

    if show_in_file_browser:
        arguments += ['--show_in_file_browser']

    if jupyter:
        arguments += ['--jupyter']

    # --- Other options

    if show_help:
        arguments += ['--help']

    # --- Run!

    path2scenes = sorted((py_file_path, sorted(set(scene.__name__ for scene in scenes_in_file)))
                         for py_file_path, scenes_in_file in
                         groupby(sorted(scenes, key=inspect.getfile), inspect.getfile))

    for py_file_path, scene_names in path2scenes:
        with sys_argv('manim', *(arguments + [py_file_path] + scene_names)):
            manimce_main()


# --- Conveniences

def is_manimgl(scene):
    try:
        from manimlib import Scene
        return issubclass(scene, Scene)
    except ImportError:
        return False


def is_manimce(scene):
    try:
        from manim import Scene
        return issubclass(scene, Scene)
    except ImportError:
        return False


def manim(*scenes, **kwargs):

    manim_ce_scenes = [scene for scene in scenes if is_manimce(scene)]
    if manim_ce_scenes:
        manimce(*manim_ce_scenes, **kwargs)

    manim_gl_scenes = [scene for scene in scenes if is_manimgl(scene)]
    if manim_gl_scenes:
        manimgl(*manim_gl_scenes, **kwargs)
