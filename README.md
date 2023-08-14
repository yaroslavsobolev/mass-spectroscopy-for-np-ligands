# Assignments for mass spectra of surface of Au nanoparticles

Data and code used in a research article (citation info and DOI will be added later).

## Dependencies:
- Python 3.7
- Numpy 1.17.3
- Scipy 1.5.2
- Matplotlib 3.1.0
- pybel (OpenBabel Python bindings)
- OpenBabel 2


## Reproducing the figures from the article
For reproducing the figures from the article's Supplementary Information file, run the
following scripts.

To reproduce Supplementary Fig. 1, run Python script `tma-esi_2cd.py`

To reproduce Supplementary Fig. 2, run Python scripts:
- `tma-mixes-esi.py` for bottom (purple) panel. Uncomment appropriate lines to get middle (orange) panel.
- `tma-esi_2cd.py` fir the main (blue) panel
- `tma-esi-isotopic-patterns.py` for insets with isotopic patterns

To reproduce Supplementary Fig. 3, run Python scripts `mua-dart-neg.py` 

To reproduce Supplementary Fig. 4, run Python script `mus-esi-neg.py`

To reproduce Supplementary Fig. 5, run Python script `mup-esi-neg.py`

To reproduce Supplementary Fig. 6, run Python script `mup-esi-neg-highMZ.py`
