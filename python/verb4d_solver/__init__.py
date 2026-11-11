# SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre for Geosciences
#
# SPDX-License-Identifier: BSD-3-Clause

"""The VERB-4D solver, installed with pip: the solver executable, its build options and the Python bindings.

``executable_path()`` gives the installed ``VERB4D_Solver``, ``build_options`` the CMake options it was built with.
The functions of the bindings (``Diffusion_2D``, ``Convection_2D``, ...) are available directly from this package.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from ._build_options import BUILD_OPTIONS as build_options

try:
    from ._core import *  # noqa: F403 (the bindings: Diffusion_2D, Convection_2D, ...)
except ModuleNotFoundError:  # built with PYTHON_BINDINGS=OFF
    pass

try:
    __version__ = version("verb4d-solver")
except PackageNotFoundError:  # not installed as a distribution
    __version__ = "unknown"

_EXECUTABLE_NAME = "VERB4D_Solver"


def executable_path() -> Path:
    """Return the path to the installed solver executable ``VERB4D_Solver``.

    The package folders (``__path__``) are searched in order: an editable install (``pip install -e``) has the
    Python files in the source tree and the built solver in site-packages.

    Returns
    -------
    Path
        Absolute path to the executable (``VERB4D_Solver.exe`` on Windows). If none is found, the path in the first
        package folder.
    """
    folders = [Path(folder).resolve() for folder in __path__]  # noqa: F405 (set for every package)
    for folder in folders:
        for name in (_EXECUTABLE_NAME, f"{_EXECUTABLE_NAME}.exe"):
            if (folder / name).is_file():
                return folder / name
    return folders[0] / _EXECUTABLE_NAME


__all__ = ["__version__", "build_options", "executable_path"]
