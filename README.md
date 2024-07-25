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
from TopoCurve import TopoCurve
from SpectralFiltering import SpectralFiltering

# Instantiate TopoCurve object with a GeoTIFF file
dem = TopoCurve('path/to/your/file.tif')

# Instantiate SpectralFiltering object with a GeoTIFF file
spectral_filter = SpectralFiltering('path/to/your/file.tif')

# Apply FFT filtering
dx, dy, filtered_elevation = spectral_filter.FFT(filter, 'filtertype', alphaIn)

# Calculate principal curvatures and curvature features
K1, K2, KM, KG = dem.CurveCalc(ZFilt, dx, dy, kt)

# Plot and save elevation values
dem.plot(input_array, 'output_image.png')
```

## Example

```
from TopoCurve import TopoCurve
from SpectralFiltering import SpectralFiltering

# Instantiate TopoCurve object with a GeoTIFF file
dem = TopoCurve('references\DEM_files\Purgatory.tif')

# Instantiate SpectralFiltering object with a GeoTIFF file
spectral_filter = SpectralFiltering('references\DEM_files\Purgatory.tif')

# Apply FFT filtering
dx, dy, filtered_elevation = spectral_filter.FFT([190, 200], 'lowpass', 0)

# Calculate principal curvatures and curvature features
K1, K2, KM, KG = dem.CurveCalc(filtered_elevation, dx, dy, 0)

# Plot and save elevation values
dem.plot(filtered_elevation, 'output_image.png')
```

## API Documentation

### TopoCurve Class

#### `TopoCurve(tiff_file)`

Initialize the TopoCurve object with a GeoTIFF file.

- **Parameters:**
  - `tiff_file` (str): Path to the GeoTIFF file.

---

#### `CurveCalc(ZFilt, dx, dy, kt)`

Calculate principal curvatures and curvature features.

- **Parameters:**

  - `ZFilt` (numpy.ndarray): Filtered surface data.
  - `dx` (float): Grid spacing in the x direction.
  - `dy` (float): Grid spacing in the y direction.
  - `kt` (float): Threshold value.

- **Returns:**
  - `K1, K2` (tuple): Principal curvatures.
  - `KM` (tuple): Mean curvature.
  - `KG` (tuple): Gaussian curvature.

---

#### `plot(input, filename)`

Plot the elevation values and save the plot as an image file to reports/figures.

- **Parameters:**
  - `input` (numpy.ndarray): Input elevation values.
  - `filename` (str): Name of the output image file.

---

### Spectral Filtering Class

#### `SpectralFiltering(tiff_file)`

Initialize the SpectralFiltering object with a GeoTIFF file.

- **Parameters:**
  - `tiff_file` (str): Path to the GeoTIFF file.

#### `detrend()`

Detrend the elevation values using least squares plane fitting.

- **Returns:**
  - `Z_detrended` (numpy.ndarray): Detrended elevation values.
  - `plane` (numpy.ndarray): Trend component of the elevation values.

#### `mirror_dem()`

Mirror the elevation values.

- **Returns:**
  - `mirrored_array` (numpy.ndarray): Mirrored elevation values.

#### `tukeyWindow(alphaIn)`

Apply a Tukey window to the elevation values.

- **Arguments:**

  - `alphaIn` (float): Parameter controlling the shape of the Tukey window.

- **Returns:**
  - `tukey_array` (numpy.ndarray): Elevation values after applying the Tukey window.

#### `padding(alphaIn)`

Pad the elevation values.

- **Arguments:**

  - `alphaIn` (float): Parameter controlling the shape of the Tukey window.

- **Returns:**
  - `padded_window_array` (numpy.ndarray): Padded elevation values.

#### `FFT(filter, filterType, alphaIn)`

Apply FFT filtering to the elevation values.

- **Arguments:**

  - `filter` (float): Filter parameter.
  - `filterType` (str): Type of filter ('lowpass' or 'highpass').
  - `alphaIn` (float): Parameter controlling the shape of the Tukey window.

- **Returns:**
  - `dx` (float): Grid spacing in the x direction.
  - `dy` (float): Grid spacing in the y direction.
  - `ZFilt` (numpy.ndarray): Filtered elevation values.

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

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

## Credits

- **Principal Contributor:** [Sonie Taylor Schermer](https://github.com/tschermer02)
- **Contributions:** [Joel Nash](https://github.com/tschermer02](https://github.com/jxnash), [Nate Klema](https://github.com/ntklema)



<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
```
