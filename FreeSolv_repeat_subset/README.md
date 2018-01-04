# FreeSolv hydration free energy calculations (repeats on a subset of molecules)
Code/tools relating to FreeSolv hydration free energ calculations as reported in the SMIRNOFF paper.
This directory repeats benchmark FreeSolv calculations on a subset of molecules (those with updated cyclobutane parameters and/or updated hydroxyl parameters).

## Manifest
* `results/`: contains a summary of the analysis.
* `scripts/`: the python and YANK scripts to reproduce the calculations and run the analysis.
* `openmmfiles.tar.gz`: Exact OpenMM files used running calculations.
* `pdbfiles.tar.gz`: Exact PDB files used running calculations.

## Installation
On Dec. 27, 2017; using Yank 0.16.0 for consistency with earlier calculations.
```
conda config --add channels omnia --add channels conda-forge
conda create --name smirnoff yank=0.16.0 openforcefield packmol 'icu=58.*'
source activate smirnoff
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
```
The final installation looks like
```
(smirnoff) conda list
alabaster                 0.7.10                   py36_1    conda-forge
ambermini                 16.16.0                       7    omnia
asn1crypto                0.22.0                   py36_0    conda-forge
babel                     2.5.1                    py36_0    conda-forge
backports                 1.0                      py36_1    conda-forge
backports.functools_lru_cache 1.4                      py36_1    conda-forge
bleach                    2.0.0                    py36_0    conda-forge
bzip2                     1.0.6                         1    conda-forge
ca-certificates           2017.11.5                     0    conda-forge
cerberus                  1.1                      py36_0    conda-forge
certifi                   2017.11.5                py36_0    conda-forge
cffi                      1.11.2                   py36_0    conda-forge
chardet                   3.0.4                    py36_0    conda-forge
clusterutils              0.3.1                    py36_0    omnia
cryptography              2.1.4                    py36_0    conda-forge
curl                      7.55.1                        0    conda-forge
cycler                    0.10.0                   py36_0    conda-forge
cython                    0.27.3                   py36_0    conda-forge
dbus                      1.10.22                       0    conda-forge
decorator                 4.1.2                    py36_0    conda-forge
docopt                    0.6.2                    py36_0    conda-forge
docutils                  0.14                     py36_0    conda-forge
entrypoints               0.2.3                    py36_1    conda-forge
expat                     2.2.5                         0    conda-forge
fftw3f                    3.3.4                         2    omnia
fontconfig                2.12.6                        0    conda-forge
freetype                  2.8.1                         0    conda-forge
gettext                   0.19.7                        1    conda-forge
glib                      2.55.0                        0    conda-forge
gmp                       6.1.2                         0    conda-forge
gst-plugins-base          1.8.0                         0    conda-forge
gstreamer                 1.8.0                         1    conda-forge
hdf4                      4.2.13                        0    conda-forge
hdf5                      1.10.1                        1    conda-forge
html5lib                  1.0.1                      py_0    conda-forge
icu                       58.2                          0    conda-forge
idna                      2.6                      py36_1    conda-forge
imagesize                 0.7.1                    py36_0    conda-forge
intel-openmp              2018.0.0             hc7b2577_8
ipykernel                 4.7.0                    py36_0    conda-forge
ipython                   6.2.1                    py36_0    conda-forge
ipython_genutils          0.2.0                    py36_0    conda-forge
ipywidgets                7.0.5                    py36_0    conda-forge
jedi                      0.10.2                   py36_0    conda-forge
jinja2                    2.10                     py36_0    conda-forge
jpeg                      9b                            2    conda-forge
jsonschema                2.6.0                    py36_0    conda-forge
jupyter                   1.0.0                    py36_0    conda-forge
jupyter_client            5.2.0                    py36_0    conda-forge
jupyter_console           5.2.0                    py36_0    conda-forge
jupyter_core              4.4.0                      py_0    conda-forge
krb5                      1.14.2                        0    conda-forge
latexcodec                1.0.4                    py36_0    conda-forge
libffi                    3.2.1                         3    conda-forge
libgcc                    7.2.0                h69d50b8_2
libgcc-ng                 7.2.0                h7cc24e2_2
libgfortran               1.0                           0
libgfortran-ng            7.2.0                h9f7466a_2
libiconv                  1.15                          0    conda-forge
libnetcdf                 4.5.0                         3    conda-forge
libpng                    1.6.34                        0    conda-forge
libsodium                 1.0.15                        1    conda-forge
libssh2                   1.8.0                         2    conda-forge
libstdcxx-ng              7.2.0                h7a57d05_2
libxcb                    1.12                          1    conda-forge
libxml2                   2.9.5                         2    conda-forge
libxslt                   1.1.32                        0    conda-forge
lxml                      4.1.1                    py36_0    conda-forge
lzo                       2.10                          0    conda-forge
markupsafe                1.0                      py36_0    conda-forge
matplotlib                2.1.1                    py36_2    conda-forge
mdtraj                    1.9.1                    py36_1    conda-forge
mistune                   0.8.3                      py_0    conda-forge
mkl                       2018.0.1             h19d6760_4
mpi4py                    3.0.0                    py36_0    conda-forge
mpich                     3.2                           5    conda-forge
nbconvert                 5.3.1                      py_1    conda-forge
nbformat                  4.4.0                    py36_0    conda-forge
ncurses                   5.9                          10    conda-forge
netcdf4                   1.3.1                    py36_2    conda-forge
networkx                  2.0                      py36_1    conda-forge
notebook                  5.2.2                    py36_1    conda-forge
numexpr                   2.6.4                    py36_0    conda-forge
numpy                     1.13.3           py36ha12f23b_0
numpydoc                  0.7.0                    py36_0    conda-forge
OpenEye-toolkits          2017.10.1                 <pip>
OpenEye-toolkits-python3-linux-x64 2017.10.1                 <pip>
openforcefield            0.0.3                    py36_1    omnia
openmm                    7.1.1                    py36_0    omnia
openmmtools               0.13.4                   py36_0    omnia
openmoltools              0.8.1                    py36_0    omnia
openssl                   1.0.2n                        0    conda-forge
oset                      0.1.3                    py36_0    conda-forge
packmol                   1!17.221                      1    omnia
pandas                    0.21.1                   py36_0    conda-forge
pandoc                    2.0.5                         0    conda-forge
pandocfilters             1.4.1                    py36_0    conda-forge
parmed                    2.7.3                    py36_1    omnia
pcre                      8.39                          0    conda-forge
pexpect                   4.3.1                    py36_0    conda-forge
pickleshare               0.7.4                    py36_0    conda-forge
pip                       9.0.1                    py36_0    conda-forge
prompt_toolkit            1.0.15                   py36_0    conda-forge
ptyprocess                0.5.2                    py36_0    conda-forge
pybtex                    0.21                     py36_0    conda-forge
pybtex-docutils           0.2.1                    py36_0    conda-forge
pycparser                 2.18                     py36_0    conda-forge
pygments                  2.2.0                    py36_0    conda-forge
pymbar                    3.0.3                    py36_1    conda-forge
pyopenssl                 17.4.0                   py36_0    conda-forge
pyparsing                 2.2.0                    py36_0    conda-forge
pyqt                      5.6.0                    py36_4    conda-forge
pysocks                   1.6.7                    py36_0    conda-forge
pytables                  3.4.2                    py36_7    conda-forge
python                    3.6.4                         0    conda-forge
python-dateutil           2.6.1                    py36_0    conda-forge
pytz                      2017.3                     py_2    conda-forge
pyyaml                    3.12                     py36_1    conda-forge
pyzmq                     16.0.2                   py36_2    conda-forge
qt                        5.6.2                         7    conda-forge
qtconsole                 4.3.1                    py36_0    conda-forge
readline                  7.0                           0    conda-forge
requests                  2.18.4                   py36_1    conda-forge
schema                    0.6.2                    py36_0    omnia
scipy                     1.0.0            py36hbf646e7_0
setuptools                38.2.4                   py36_0    conda-forge
simplegeneric             0.8.1                    py36_0    conda-forge
sip                       4.18                     py36_1    conda-forge
six                       1.11.0                   py36_1    conda-forge
snowballstemmer           1.2.1                    py36_0    conda-forge
sphinx                    1.6.5                    py36_0    conda-forge
sphinxcontrib-bibtex      0.3.6                    py36_0    conda-forge
sphinxcontrib-websupport  1.0.1                    py36_0    conda-forge
sqlite                    3.20.1                        2    conda-forge
terminado                 0.8.1                    py36_0    conda-forge
testpath                  0.3.1                    py36_0    conda-forge
tk                        8.6.7                         0    conda-forge
tornado                   4.5.2                    py36_0    conda-forge
traitlets                 4.3.2                    py36_0    conda-forge
typing                    3.6.2                    py36_0    conda-forge
urllib3                   1.22                     py36_0    conda-forge
wcwidth                   0.1.7                    py36_0    conda-forge
webencodings              0.5                      py36_0    conda-forge
wheel                     0.30.0                     py_1    conda-forge
widgetsnbextension        3.1.0                    py36_0    conda-forge
xorg-libxau               1.0.8                         3    conda-forge
xorg-libxdmcp             1.1.2                         3    conda-forge
xz                        5.2.3                         0    conda-forge
yaml                      0.1.6                         0    conda-forge
yank                      0.16.0                   py36_0    omnia
zeromq                    4.2.1                         1    conda-forge
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
