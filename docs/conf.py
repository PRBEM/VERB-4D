# SPDX-FileCopyrightText: 2026 GFZ Helmholtz Centre for Geosciences
# SPDX-FileContributor: Sahil Jhawar
#
# SPDX-License-Identifier: BSD-3-Clause

"""Sphinx configuration for the VERB-4D C++ API documentation (Breathe + Exhale)."""

from pathlib import Path

docs_dir = Path(__file__).resolve().parent

project = "VERB-4D"
copyright = "2026, GFZ Helmholtz Centre for Geosciences"
author = "GFZ Helmholtz Centre for Geosciences"

master_doc = "index"
exclude_patterns = ["_build", "html", "latex", "api/**/.doctrees"]

extensions = [
    "breathe",
    "exhale",
    "myst_parser",
]

breathe_projects = {"VERB-4D": str(docs_dir / "_doxygen" / "xml")}
breathe_default_project = "VERB-4D"

# Exhale runs doxygen itself, so the existing Doxyfile stays usable on its own.
exhale_args = {
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "doxygenStripFromPath": "..",
    "rootFileTitle": "C++ API",
    "createTreeView": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": "INPUT = ../src\nRECURSIVE = YES\nEXCLUDE = ../src/python_bindings\nEXTRACT_ALL = YES",
}

primary_domain = "cpp"
highlight_language = "cpp"

html_theme = "pydata_sphinx_theme"
