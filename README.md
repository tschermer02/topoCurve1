# TopoCurve

TopoCurve is a Python library for processing digital elevation models (DEM) stored in GeoTIFF format. This library provides functionalities to extract metadata, calculate principal curvatures and curvature features, as well as plot elevation values. It also includes spectral filtering capabilities for advanced DEM processing.

## Installation

To install TopoCurve, simply clone the repository and install the dependencies listed in `requirements.txt`:

```
git clone https://github.com/username/topo_curve.git
cd topo_curve
pip install -r requirements.txt
```

## Usage

```
from topo_curve import TopoCurve
from topo_curve.spectralfiltering import SpectralFiltering

# Instantiate TopoCurve object with a GeoTIFF file
dem = TopoCurve('path/to/your/file.tif')

# Calculate principal curvatures and curvature features
K1, K2, KM, KG = dem.CurveCalc(ZFilt, dx, dy, kt)

# Plot and save elevation values
dem.plot(input_array, 'output_image.png')

# Instantiate SpectralFiltering object with a GeoTIFF file
spectral_filter = SpectralFiltering('path/to/your/file.tif')

# Apply FFT filtering
dx, dy, filtered_elevation = spectral_filter.FFT(filter, filterType, alphaIn)
```

## Project Organization

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io

---

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
```
