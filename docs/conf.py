"""Sphinx configuration for the immune documentation website."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import immune

project = "immune"
author = "Tim Holy"
copyright = "2026, immune contributors"
version = release = immune.__version__

needs_sphinx = "9.1"
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_typehints_format = "short"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = False

myst_enable_extensions = [
    "attrs_inline",
    "colon_fence",
    "deflist",
    "fieldlist",
    "substitution",
]
myst_heading_anchors = 3
myst_substitutions = {"version": release}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "skbio": ("https://scikit.bio/docs/latest", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable", None),
    "scirpy": ("https://scirpy.scverse.org/en/latest", None),
}

html_theme = "pydata_sphinx_theme"
html_title = f"immune {release}"
html_logo = "_static/immune-logo.svg"
html_favicon = "_static/favicon.svg"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_last_updated_fmt = "%Y-%m-%d"
html_last_updated_use_utc = True

html_theme_options = {
    "navbar_align": "left",
    "header_links_before_dropdown": 6,
    "show_nav_level": 2,
    "show_toc_level": 2,
    "navigation_depth": 4,
    "navigation_with_keys": True,
    "search_as_you_type": True,
    "back_to_top_button": True,
    "pygments_light_style": "tango",
    "pygments_dark_style": "github-dark",
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_persistent": ["search-button"],
    "icon_links": [
        {
            "name": "GitHub repository",
            "url": "https://github.com/ScienceAdvances/immune",
            "icon": "fa-brands fa-github",
        }
    ],
    "secondary_sidebar_items": ["page-toc", "sourcelink"],
    "footer_start": ["copyright", "sphinx-version"],
    "footer_end": ["theme-version"],
}

html_sidebars = {"index": []}

copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.?: "
copybutton_prompt_is_regexp = True


def setup(app):
    """Register website metadata consumed by templates and extensions."""

    app.add_config_value("immune_release", release, "html")
