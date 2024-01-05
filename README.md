
# fastsum-gpu

This repository is for the "FAST" FASTSUM Analysis Software Tools prepared during the 2024 DiRAC GPU Workshop.


---
# Pre-commit notes
`pre-commit` will be installed. Activate the conda environment (below) and run
`pre-commit install`

# Conda Notes
Install your favourite conda solution, such as https://docs.conda.io/en/latest/miniconda.html
## Switch to a faster environment solver
This is optional, but likely will solve the dependencies much much faster. See https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community
```conda update -n base conda```
```conda install -n base conda-libmamba-solver```
```conda config --set solver libmamba```

## Install Environment
```conda env create -f linux_environment.yml```
## Activate/Use
```conda activate fastsum```
## Update (w. new packages)

1. Edit `environment.yml`
2. Deactivate conda environment with `conda deactivate`
3. Update conda environment with `conda env update -f=linux_environment.yml`
