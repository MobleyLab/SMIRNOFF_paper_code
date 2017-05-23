# SMIRNOFF paper code
Code/tools relating to the initial SMIRNOFF format paper.

## Installation
On May 20, 2017, after adding the `http://conda.anaconda.org/omnia` channel.
```
conda create --name smirnoff yank openforcefield packmol
```
I've removed and reinstalled `yank` (hash `8f0798e`), `openmmtools` (from PR#200, hash `3baac91`),
and `openforcefield` (hash `14a9021`) from source code. Then
```
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
```
The final installation looks like
```
(smirnoff) bash$ conda list
alabaster                 0.7.10                   py35_0
ambermini                 16.16.0                       7    omnia
babel                     2.4.0                    py35_0
cairo                     1.14.8                        0
clusterutils              0.3.0                    py35_0    omnia
curl                      7.45.0                        0
cycler                    0.10.0                   py35_0
cython                    0.25.2                   py35_0
decorator                 4.0.11                   py35_0
docopt                    0.6.2                    py35_1    omnia
docutils                  0.13.1                   py35_0
fftw3f                    3.3.4                         2    omnia
fontconfig                2.12.1                        3
freetype                  2.5.5                         2
glib                      2.50.2                        1
harfbuzz                  0.9.39                        2
hdf4                      4.2.12                        0
hdf5                      1.8.15.1                      3
imagesize                 0.7.1                    py35_0
jinja2                    2.9.6                    py35_0
jpeg                      8d                            2
latexcodec                1.0.1                    py35_0    omnia
libffi                    3.2.1                         1
libgcc                    5.2.0                         0
libgfortran               1.0                           0
libiconv                  1.14                          0
libnetcdf                 4.3.3.1                       3
libpng                    1.6.27                        0
libxml2                   2.9.4                         0
libxslt                   1.1.29                        0
lxml                      3.7.3                    py35_0
markupsafe                0.23                     py35_2
matplotlib                1.5.1               np110py35_0
mdtraj                    1.8.0               np110py35_1    omnia
mkl                       11.3.1                        0
mpi4py                    2.0.0                    py35_2    omnia
mpich                     3.2                           4    omnia
netcdf4                   1.2.2               np110py35_0
networkx                  1.11                     py35_0
nose                      1.3.7                    py35_1
numexpr                   2.5.2               np110py35_0
numpy                     1.10.4                   py35_1
numpydoc                  0.6.0                    py35_0
OpenEye-toolkits          2017.2.1                  <pip>
OpenEye-toolkits-python3-linux-x64 2017.2.1                  <pip>
openmm                    7.1.0rc1                 py35_0    omnia
openmmtools               0.9.5                     <pip>
openmoltools              0.7.5                    py35_0    omnia
openssl                   1.0.2k                        2
oset                      0.1.3                    py35_1    omnia
packmol                   2016.06.09                    1    omnia
pandas                    0.18.1              np110py35_0
pango                     1.40.3                        1
parmed                    2.7.3                    py35_1    omnia
pcre                      8.39                          1
pip                       9.0.1                    py35_1
pixman                    0.34.0                        0
pybtex                    0.18                     py35_0    omnia
pybtex-docutils           0.2.1                    py35_1    omnia
pygments                  2.2.0                    py35_0
pymbar                    3.0.1.beta0         np110py35_0    omnia
pyparsing                 2.1.4                    py35_0
pyqt                      4.11.4                   py35_4
pytables                  3.2.2               np110py35_1
python                    3.5.3                         1
python-dateutil           2.6.0                    py35_0
pytz                      2017.2                   py35_0
pyyaml                    3.12                     py35_0
qt                        4.8.7                         4
readline                  6.2                           2
requests                  2.14.2                   py35_0
schema                    0.6.2                    py35_0    omnia
scipy                     0.17.0              np110py35_1
setuptools                27.2.0                   py35_0
sip                       4.18                     py35_0
six                       1.10.0                   py35_0
snowballstemmer           1.2.1                    py35_0
sphinx                    1.5.6                    py35_0
sphinxcontrib-bibtex      0.3.2                    py35_0    omnia
sqlite                    3.13.0                        0
tk                        8.5.18                        0
wheel                     0.29.0                   py35_0
xz                        5.2.2                         1
yaml                      0.1.6                         0
yank                      0.16.0                    <pip>
zlib                      1.2.8                         3
```

## Create setup files
```
cd scripts/
python create_input_files.py > create_input_files.log
bash run-all-torque.sh
```