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
and interpreting landscape features and their formations. The curvature and attributes of surface features provide insights into geological processes and 
environmental conditions. `TopoCurve` is a Python package designed to analyze 
and process DEMs by calculating various curvature attributes, including principal,
Gaussian, and mean curvatures. This tool facilitates the extraction and 
classification of surface attributes, aiding in the interpretation of topographic data.
The package leverages efficient numerical methods and advanced processing techniques 
to handle large datasets, making it a valuable asset for researchers and practitioners 
in geomatics, geography, and earth sciences.

# Statement of need

`TopoCurve` addresses a significant need in the field of topographic surface analysis by providing a specialized Python package for calculating and visualizing curvature attributes. Python's versatility allows for the integration of advanced numerical methods and algorithms while maintaining an accessible and user-friendly interface. `TopoCurve` provides a comprehensive set of tools for analyzing terrain attributes, such as curvature and surface classification, through a well-structured API. This includes functionalities for calculating principal curvatures, Gaussian and mean curvatures, and generating curvature maps based on elevation data. This tool is particularly valuable for researchers and professionals in geomorphology, volcanology, and related fields who require precise and computationally efficient methods for analyzing surface topography.

The mathematical foundation of `TopoCurve` is rooted in the classical work of Gauss (1827), which laid the groundwork for differential geometry and Riemannian manifolds. By using Gauss's method for curvature derivation, `TopoCurve` offers a computationally faster alternative to modern matrix-based approaches, making it suitable for large datasets and real-time applications. The package focuses on the derivation of intrinsic and extrinsic curvatures of topographic surfaces and includes functionalities for spectral filtering and scale selection.

Beyond theoretical analysis, `TopoCurve` is applied to practical problems, such as mapping groundwater in the Oregon Cascades, and has attracted interest from high-level geomorphologists and volcanologists eager to leverage its capabilities. The package builds on existing work in the spectral domain, including the approach of Perron et al. (2008), and incorporates rigorous methods for filtering and characterizing topography. Additionally, it draws inspiration from the curvature classification methods of Bergbauer and Pollard (2001), adapted to suit the needs of topographic analysis.

By providing a robust, user-friendly interface for these advanced calculations, `TopoCurve` enables researchers to explore and interpret topographic features with greater accuracy and efficiency, facilitating new insights and applications in the study of Earth's surface processes.


# Mathematics

### Intrinsic Curvature

The intrinsic curvature of a surface can be computed using the following formulas:

1. **Principal Curvatures $K_1$ and $K_2$**:

   Principal curvatures $K_1$ and $K_2$ are determined from the coefficients of the second fundamental form and the first fundamental form. The principal curvatures are given by:

   $$K_1,K_2=\frac{-b\pm\sqrt{b^2-4ac}}{2a}$$

   with:
   
   $$a=E\cdot G-F^2$$
   
   $$b=-(g\cdot E-2\cdot f\cdot F+e\cdot G)$$
   
   $$c=e\cdot g-f^2$$

3. **Gaussian Curvature $K_G$**:

   The Gaussian curvature $K_G$ is the product of the principal curvatures:

   $$K_G=K_1\cdot K_2$$

4. **Mean Curvature $K_M$**:

   The mean curvature $K_M$ is the average of the principal curvatures:

   $$K_M=\frac{1}{2}\left(K_1+K_2\right)$$

### Extrinsic Curvature

The extrinsic curvature of a surface involves normal vectors and their derivatives. The calculations include:

1. **Normal Vector Derivatives**:

   The normal vectors $\mathbf{N}_x$ and $\mathbf{N}_y$ are derived from the surface normal components:

   $$\mathbf{N}=\frac{\mathbf{N}_x}{\|\mathbf{N}\|},\quad \mathbf{N}_y=\frac{\mathbf{N}_y}{\|\mathbf{N}\|}$$

   where $\mathbf{N}_x$ and $\mathbf{N}_y$ are the partial derivatives of the normal vector with respect to $x$ and $y$, respectively.

2. **Second Fundamental Form Coefficients**:

   The coefficients $e, f, g$ of the second fundamental form are calculated from:

   $$e=-\left(\mathbf{N}_x\cdot\mathbf{S}_U\right)$$

   $$f=-0.5\left(\mathbf{N}_x\cdot\mathbf{S}_V+\mathbf{N}_y\cdot\mathbf{S}_U\right)$$

   $$g=-\left(\mathbf{N}_y\cdot\mathbf{S}_V\right)$$

   where $\mathbf{S}_U$ and $\mathbf{S}_V$ are surface derivatives in the $u$ and $v$ directions, respectively.

3. **Mean Curvature $H$**:

   The mean curvature $H$ is derived from:

   $$H=\frac{1}{2}\left(\frac{E\cdot N_x+F\cdot (N_x+N_y)+G\cdot N_y}{E\cdot G-F^2}\right)$$

   where $N_x$ and $N_y$ are the partial derivatives of the normal vector with respect to $x$ and $y$.

These equations and coefficients are fundamental to the analysis of surface geometry and attributes in digital elevation models.


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
