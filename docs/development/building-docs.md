# Building the documentation

## Install

Documentation builds require Python 3.12 or newer because Sphinx 9.1 requires
it.

```bash
python -m pip install -e '.[docs]'
```

## Strict HTML build

```bash
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

Or use the documentation Makefile:

```bash
make -C docs html
```

Warnings fail the build. Resolve missing references, malformed directives and
autodoc import failures before publishing.

## Link checking

```bash
make -C docs linkcheck
```

External services can be transient, so investigate failures before changing or
disabling a reference.

## Read the Docs

The repository includes `.readthedocs.yaml`. Read the Docs installs the `docs`
extra, builds `docs/conf.py` with Python 3.12 and fails on Sphinx warnings.

## Theme and footer

The site uses Sphinx 9.1 and PyData Sphinx Theme 0.20. The footer includes the
actual Sphinx and theme versions through the theme's built-in `sphinx-version`
and `theme-version` templates.

Brand styling is isolated in `docs/_static/custom.css`; content and API pages
remain standard MyST Markdown/reStructuredText so future theme upgrades do not
require rewriting the documentation.
