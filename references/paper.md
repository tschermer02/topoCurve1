---
title: "TopoCurve: A Python Package for Advanced Digital Elevation Model Analysis"
tags:
  - Python
  - Digital Elevation Models
  - Curvature Analysis
  - Geographic Information Systems (GIS)
  - Topographic Modeling
  - Geospatial Analysis
  - Environmental Modeling
authors:
- name: Sonie T. Schermer
  orcid: 0009-0005-9439-4470
  equal-contrib: false
  affiliation: 1
- name: Joel Nash
  equal-contrib: false
  affiliation: 1
- name: Nate Klema
- equal-contrib: false
  affiliation: 1
affiliations:
  - name: Fort Lewis College, United States of America
    index: 1
date: 25 July 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The analysis of digital elevation models (DEMs) is crucial for understanding 
and interpreting landscape features and their formations. The curvature and 
surface attributes of surface features provide insights into geological processes and 
environmental conditions. `TopoCurve` is a Python package designed to analyze 
and process DEMs by calculating various curvature attributes, including principal,
Gaussian, and mean curvatures. This tool facilitates the extraction and 
classification of surface attributes, aiding in the interpretation of topographic data.
The package leverages efficient numerical methods and advanced processing techniques 
to handle large datasets, making it a valuable asset for researchers and practitioners 
in geomatics, geography, and earth sciences.

# Statement of need

`TopoCurve` addresses a significant need in the field of topographic surface analysis by providing a specialized Python package for calculating and visualizing curvature attributes. Python's versatility allows for the integration of advanced numerical methods and algorithms while maintaining an accessible and user-friendly interface. `TopoCurve` provides a comprehensive set of tools for analyzing terrain attributes, such as curvature and surface classification, through a well-structured API. This includes functionalities for calculating principal curvatures, Gaussian and mean curvatures, and generating curvature maps based on elevation data. This tool is particularly valuable for researchers and professionals in geomorphology, volcanology, and related fields who require precise and computationally efficient methods for analyzing surface topography.

The mathematical foundation of `TopoCurve` is rooted in the classical work of Gauss (1827), which laid the groundwork for differential geometry and Riemannian manifolds. By using Gauss's method for curvature derivation, TopoCurve offers a computationally faster alternative to modern matrix-based approaches, making it suitable for large datasets and real-time applications. The package focuses on the derivation of intrinsic and extrinsic curvatures of topographic surfaces and includes functionalities for spectral filtering and scale selection.

Beyond theoretical analysis, TopoCurve is applied to practical problems, such as mapping groundwater in the Oregon Cascades, and has attracted interest from high-level geomorphologists and volcanologists eager to leverage its capabilities. The package builds on existing work in the spectral domain, including the approach of Perron et al. (2008), and incorporates rigorous methods for filtering and characterizing topography. Additionally, it draws inspiration from the curvature classification methods of Bergbauer and Pollard (2001), adapted to suit the needs of contemporary topographic analysis.

By providing a robust, user-friendly interface for these advanced calculations, TopoCurve enables researchers to explore and interpret topographic features with greater accuracy and efficiency, facilitating new insights and applications in the study of Earth's surface processes.


# Mathematics

\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int\_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:

- `@author:2001` -> "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
