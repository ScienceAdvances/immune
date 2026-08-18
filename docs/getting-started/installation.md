# Installation

## Base package

`immune` requires Python 3.10 or newer. The base installation contains readers,
the canonical data model, clonotype construction, longitudinal statistics and
dependency-light repertoire summaries.

```bash
python -m pip install immune
```

For local development:

```bash
git clone <your-repository-url>
cd immune
python -m pip install -e .
```

## Optional analysis families

Install only the mature backends required by a study.

| Extra | Install command | Main capabilities |
|---|---|---|
| `bulk` | `pip install 'immune[bulk]'` | scikit-bio diversity and downsampling |
| `differential` | `pip install 'immune[differential]'` | PyDESeq2 clonotype differential abundance |
| `singlecell` | `pip install 'immune[singlecell]'` | AnnData, Scanpy, MuData and Scirpy |
| `advanced` | `pip install 'immune[advanced]'` | scvi-tools, Pertpy Milo and decoupler |
| `plot` | `pip install 'immune[plot]'` | Matplotlib and seaborn bulk plots |
| `all` | `pip install 'immune[all]'` | All analysis families |

Most bulk studies will begin with:

```bash
python -m pip install 'immune[bulk,differential,plot]'
```

A joint scRNA-seq/TCR study can use:

```bash
python -m pip install 'immune[bulk,singlecell,advanced,plot]'
```

## Check optional backends

```python
import immune as iu

iu.backend_status()
```

Optional packages are imported only when their functions are called. A missing
backend therefore does not prevent the base package from importing.

## Build the documentation website

The documentation toolchain requires Python 3.12 or newer. The current website
targets Sphinx 9.1 and PyData Sphinx Theme 0.20.

```bash
python -m pip install -e '.[docs]'
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` in a browser to inspect the local website.
