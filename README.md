[![DOI](https://zenodo.org/badge/754718683.svg)](https://zenodo.org/doi/10.5281/zenodo.11612296)

# Code and data for "Vector data cubes for features evolving in space and time"

## Short paper at AGILE 2024:

Abad, L., Sudmanns, M., and Hölbling, D.: Vector data cubes for features evolving in space and time, AGILE GIScience Ser., 5, 16, <https://doi.org/10.5194/agile-giss-5-16-2024>, 2024.

## Software documentation:

### [Reproduce manuscript figures](https://loreabad6.github.io/vdc-space-time-feats/code/manuscript-figures.html)

The directory [`/code`](/code) contains a Quarto file [manuscript-figures.qmd](code/manuscript-figures.qmd) to recreate the figures from the manuscript.
The data for the figures comes from Zenodo and is downloaded with the code above.
This file has been rendered to a Quarto output ([manuscript-figures.html](https://loreabad6.github.io/vdc-space-time-feats/code/manuscript-figures.html)) or [manuscript-figures.md](code/manuscript-figures.md) to easily explore the rendered code. 

### [VDC showcase notebook](https://loreabad6.github.io/vdc-space-time-feats/notebook/vdc-showcase.html)

The directory [`/notebook`](/notebook) contains a computational notebook written in Quarto [vdc-showcase.qmd](notebook/vdc-showcase.qmd) that documents the creation of vector data cubes (VDC) for features that evolve in space and time, namely, lava flows and landslides, as presented in the manuscript. The notebook contains further examples of the usage of VDCs, beyond the ones included in the manuscript. See [vdc-showcase.html](https://loreabad6.github.io/vdc-space-time-feats/notebook/vdc-showcase.html) or  [vdc-showcase.md](notebook/vdc-showcase.md)for a rendered version.

At the end of both files, the output of `sessioninfo::sessioninfo()` is included to document set-up information such as software and packages used.

## [Manuscript](https://loreabad6.github.io/vdc-space-time-feats/manuscript/manuscript.html)

## [Presentation](https://loreabad6.github.io/vdc-space-time-feats/presentation/agile24.html#/title-slide)

## Data sources:

- Hölbling, Daniel, Lorena Abad, Zahra Dabiri, Prasicek, and Anne-Laure Argentin. 2024. “Butangbunasi Landslide and Landslide-Dammed Lake Outlines Based on Landsat Time Series with Respect to Typhoons.” Zenodo. https://doi.org/10.5281/zenodo.10635102.

- Pedersen, Gro B. M., Joaquin M. C. Belart, Birgir Vilhelm Óskarsson, Magnús Tumi Gudmundsson, Nils Gies, Thórdís Högnadóttir, Ásta Rut Hjartardóttir, et al. 2023. “Digital Elevation Models, Orthoimages and Lava Outlines of the 2021 Fagradalsfjall Eruption: Results from Near Real-Time Photogrammetric Monitoring.” Zenodo. https://doi.org/10.5281/zenodo.7866738.
