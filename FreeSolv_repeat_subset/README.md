# FreeSolv hydration free energy calculations (repeats on a subset of molecules)
Code/tools relating to FreeSolv hydration free energ calculations as reported in the SMIRNOFF paper.
This directory repeats benchmark FreeSolv calculations on a subset of molecules (those with updated cyclobutane parameters and/or updated hydroxyl parameters).

## Manifest
* `results/`: contains a summary of the analysis.
* `scripts/`: the python and YANK scripts to reproduce the calculations and run the analysis.

## Installation
On Dec. 26, 2017.
```
conda config --add channels omnia --add channels conda-forge
conda create --name smirnoff yank openforcefield packmol 'icu=58.*'
source activate smirnoff
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
```
The final installation looks like
```
(smirnoff) bash$ conda list
alabaster                 0.7.10                   py35_0    conda-forge
ambermini                 16.16.0                       7    omnia
babel                     2.3.4                    py35_0    conda-forge
ca-certificates           2017.4.17                     0    conda-forge
certifi                   2017.4.17                py35_0    conda-forge
clusterutils              0.3.0                    py35_0    omnia
curl                      7.52.1                        0    conda-forge
cycler                    0.10.0                   py35_0    conda-forge
cython                    0.25.2                   py35_1    conda-forge
decorator                 4.0.11                   py35_0    conda-forge
docopt                    0.6.2                    py35_0    conda-forge
docutils                  0.13.1                   py35_0    conda-forge
fftw3f                    3.3.4                         2    omnia
fontconfig                2.12.1                        1    conda-forge
freetype                  2.6.3                         1    conda-forge
hdf4                      4.2.12                        0    conda-forge
hdf5                      1.8.17                       11    conda-forge
icu                       58.1                          1    conda-forge
imagesize                 0.7.1                    py35_0    conda-forge
jinja2                    2.9.5                    py35_0    conda-forge
jpeg                      9b                            0    conda-forge
latexcodec                1.0.4                    py35_0    conda-forge
libgcc                    5.2.0                         0
libgfortran               1.0                           0
libiconv                  1.14                          4    conda-forge
libnetcdf                 4.4.1.1                       4    conda-forge
libpng                    1.6.28                        0    conda-forge
libxml2                   2.9.4                         4    conda-forge
libxslt                   1.1.29                        4    conda-forge
lxml                      3.7.3                    py35_0    conda-forge
markupsafe                1.0                      py35_0    conda-forge
matplotlib                2.0.0b3             np110py35_3    conda-forge
mdtraj                    1.8.0               np110py35_1    conda-forge
mkl                       11.3.1                        0
mpi4py                    2.0.0                    py35_2    conda-forge
mpich                     3.2                           4    conda-forge
ncurses                   5.9                          10    conda-forge
netcdf4                   1.2.7               np110py35_0    conda-forge
networkx                  1.11                     py35_0    conda-forge
numexpr                   2.6.1               np110py35_1    conda-forge
numpy                     1.10.4                   py35_1
numpydoc                  0.6.0                    py35_0    conda-forge
OpenEye-toolkits          2017.2.1                  <pip>
OpenEye-toolkits-python3-linux-x64 2017.2.1                  <pip>
openforcefield            0.0.2                    py35_0    omnia
openmm                    7.1.1                    py35_0    omnia
openmmtools               0.10.0                   py35_0    omnia
openmoltools              0.8.1                    py35_0    omnia
openssl                   1.0.2k                        0    conda-forge
oset                      0.1.3                    py35_1    omnia
packmol                   2016.06.09                    1    omnia
pandas                    0.19.2              np110py35_1    conda-forge
parmed                    2.7.3                    py35_1    omnia
pip                       9.0.1                    py35_0    conda-forge
pybtex                    0.18                     py35_0    omnia
pybtex-docutils           0.2.1                    py35_1    omnia
pygments                  2.2.0                    py35_0    conda-forge
pymbar                    3.0.1.beta0         np110py35_0    omnia
pyparsing                 2.2.0                    py35_0    conda-forge
pyqt                      4.11.4                   py35_2    conda-forge
pytables                  3.3.0               np110py35_0    conda-forge
python                    3.5.3                         3    conda-forge
python-dateutil           2.6.0                    py35_0    conda-forge
pytz                      2017.2                   py35_0    conda-forge
pyyaml                    3.12                     py35_1    conda-forge
qt                        4.8.7                         3
readline                  6.2                           0    conda-forge
requests                  2.13.0                   py35_0    conda-forge
schema                    0.6.2                    py35_0    omnia
scipy                     0.17.0              np110py35_1
setuptools                33.1.1                   py35_0    conda-forge
sip                       4.18                     py35_1    conda-forge
six                       1.10.0                   py35_1    conda-forge
snowballstemmer           1.2.1                    py35_0    conda-forge
sphinx                    1.6.2                    py35_0    conda-forge
sphinxcontrib-bibtex      0.3.2                    py35_0    omnia
sphinxcontrib-websupport  1.0.1                    py35_0    conda-forge
sqlite                    3.13.0                        1    conda-forge
tk                        8.5.19                        1    conda-forge
typing                    3.6.1                    py35_0    conda-forge
wheel                     0.29.0                   py35_0    conda-forge
xz                        5.2.2                         0    conda-forge
yaml                      0.1.6                         0    conda-forge
yank                      0.16.0                   py35_0    omnia
zlib                      1.2.11                        0    conda-forge
```

## Create setup files and run YANK
```
cd scripts/
python create_input_files.py > create_input_files.log
qsub -t 1-40 -v n_jobs=40 run-torque.sh  # or bash run-all-torque.sh
```

## Author
Credit for these scripts/tools goes to Andrea Rizzi (Chodera lab) and David Mobley (UCI)
