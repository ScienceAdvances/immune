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
git clone https://github.com/ScienceAdvances/immune.git
cd immune
python -m pip install -e .
```

## Optional analysis families

Install only the mature backends required by a study.

| Extra | Install command | Main capabilities |
|---|---|---|
| `bulk` | `pip install 'immune[bulk]'` | scikit-bio diversity and downsampling |
| `differential` | `pip install 'immune[differential]'` | PyDESeq2 clonotype differential abundance |
| `singlecell` | `pip install 'immune[singlecell]'` | AnnData, MuData and Scirpy |
| `plot` | `pip install 'immune[plot]'` | Matplotlib and seaborn bulk plots |
| `all` | `pip install 'immune[all]'` | All analysis families |

Most bulk studies will begin with:

```bash
python -m pip install 'immune[bulk,differential,plot]'
```

A joint scRNA-seq/TCR study can use:

```bash
python -m pip install 'immune[bulk,singlecell,plot]'
```

## Check optional backends

```python
import immune as iu

iu.backend_status()
```

Optional packages are imported only when their functions are called. A missing
backend therefore does not prevent the base package from importing.

RNA analysis is installed separately through cellscope. immune has no RNA or R
analysis extra and no dependency on cellscope.
