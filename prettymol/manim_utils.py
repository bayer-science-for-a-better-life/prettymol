"""Helpers to deal with manim."""
import inspect
import sys
from contextlib import contextmanager
from itertools import groupby
from pathlib import Path

from manimlib.__main__ import main as manimgl_main
from manim.__main__ import main as manimce_main


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


def manimgl(*scenes,
            write: bool = False):
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

    arguments = []

    if write:
        arguments += ['--write_file']

    media_dir = Path(__file__).parent.parent / 'media'

    path2scenes = sorted((py_file_path, sorted(set(scene.__name__ for scene in scenes_in_file)))
                         for py_file_path, scenes_in_file in
                         groupby(sorted(scenes, key=inspect.getfile), inspect.getfile))

    for py_file_path, scene_names in path2scenes:
        with sys_argv('manimgl', *(arguments + [py_file_path] + scene_names)):
            manimgl_main()


def manim_community(*scenes):
    """
    Python friendly manim caller.

    This currently has been tested with manim.community version 1.0.

    From manim --help
    -----------------------------------
    Manim Community v0.4.0
    usage: manim file [flags] [scene [scene ...]]
           manim {cfg,init,plugins} [opts]

    Animation engine for explanatory math videos

    positional arguments:
      file                  Path to file holding the python code for the scene
      scene_names           Name of the Scene class you want to see

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_FILE, --output_file OUTPUT_FILE
                            Specify the name of the output file, if it should be different from the scene class name
      -p, --preview         Automatically open the saved file once its done
      -f, --show_in_file_browser
                            Show the output file in the File Browser
      --sound               Play a success/failure sound
      --leave_progress_bars
                            Leave progress bars displayed in terminal
      -a, --write_all       Write all the scenes from a file
      -w, --write_to_movie  Render the scene as a movie file (this is on by default)
      -s, --save_last_frame
                            Save the last frame only (no movie file is generated)
      -g, --save_pngs       Save each frame as a png
      -i, --save_as_gif     Save the video as gif
      --disable_caching     Disable caching (will generate partial-movie-files anyway)
      --flush_cache         Remove all cached partial-movie-files
      --log_to_file         Log terminal output to file
      -c BACKGROUND_COLOR, --background_color BACKGROUND_COLOR
                            Specify background color
      --media_dir MEDIA_DIR
                            Directory to store media (including video files)
      --log_dir LOG_DIR     Directory to store log files
      --tex_template TEX_TEMPLATE
                            Specify a custom TeX template file
      --dry_run             Do a dry run (render scenes but generate no output files)
      -t, --transparent     Render a scene with an alpha channel
      -q {k,p,h,m,l}, --quality {k,p,h,m,l}
                            Render at specific quality, short form of the --*_quality flags
      --low_quality         Render at low quality
      --medium_quality      Render at medium quality
      --high_quality        Render at high quality
      --production_quality  Render at default production quality
      --fourk_quality       Render at 4K quality
      -l                    DEPRECATED: USE -ql or --quality l
      -m                    DEPRECATED: USE -qm or --quality m
      -e                    DEPRECATED: USE -qh or --quality h
      -k                    DEPRECATED: USE -qk or --quality k
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution, passed as "height,width". Overrides the -l, -m, -e, and -k flags, if present
      -n FROM_ANIMATION_NUMBER, --from_animation_number FROM_ANIMATION_NUMBER
                            Start rendering at the specified animation index, instead of the first animation.
                            If you pass in two comma separated values, e.g. '3,6', it will end the rendering
                            at the second value
      --use_webgl_renderer  Render animations using the WebGL frontend
      --webgl_renderer_path WEBGL_RENDERER_PATH
                            Path to the WebGL frontend
      --webgl_updater_fps WEBGL_UPDATER_FPS
                            Frame rate to use when generating keyframe data for animations that use updaters
                            while using the WebGL frontend
      --config_file CONFIG_FILE
                            Specify the configuration file
      --custom_folders      Use the folders defined in the [custom_folders] section of the config file to define
                            the output folder structure
      -v {DEBUG,INFO,WARNING,ERROR,CRITICAL}, --verbosity {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                            Verbosity level. Also changes the ffmpeg log level unless the latter is specified
                            in the config
      --version             Print the current version of Manim you are using
      --progress_bar True/False
                            Display the progress bar

    Made with <3 by the ManimCommunity devs

    """

    arguments = []

    path2scenes = sorted((py_file_path, sorted(set(scene.__name__ for scene in scenes_in_file)))
                         for py_file_path, scenes_in_file in
                         groupby(sorted(scenes, key=inspect.getfile), inspect.getfile))

    for py_file_path, scene_names in path2scenes:
        with sys_argv('manimgl', *(arguments + [py_file_path] + scene_names)):
            manimce_main()


if __name__ == '__main__':
    from prettymol.scenes.language_of_life import LoLCommonsIntro
    manim_community(LoLCommonsIntro)
