# Photomap Tool

Photomap is an open source tool for extraction of nick and backbone signals from provided Bionano optical map tiff images. Photomap is able to transform nick signals into nicking site coordinates and write in BNX format.


## Requirements

**`Python >3.5`** and 
**`sqlite`**

Required Python libraries:

* [Numpy and Scipy](http://www.numpy.org)
* [Numba](http://numba.pydata.org)
* [Sqlalchemy](https://www.sqlalchemy.org)
* [intervaltree](https://pypi.org/project/intervaltree)
* [scikit-image](https://scikit-image.org)


## Installation

```bash
pip install numpy numba scipy sqlalchemy pylab scikit-image intervaltree
git clone https://gitlab.com/akdel/OptiScan.git
cd photomap
python setup.py install
cd pipelines
chmod +x ./extract_molecules
chmod +x ./write_bnx
```

## Usage

`optiscan/pipelines/` folder contains two pipeline scripts for running OptiScan.


1. `extract_molecules` is for extracting and recording OM signals into disk.

    arguments: 
    ```bash
    <list of runs: use commas to separate> 
    <chip dimension: x,y> 
    <database name: str>
    <number of threads: int> 
    <organism name: str>
    ```

    example usage:

    ```bash
    ./extract_molecules '/path/to/run1,/path/to/run/2' '12,95' apple.db 10 apple
    ```

2. `write_bnx` is for writing extracted molecules in `.bnx` format.
    
    arguments:
    ```bash 
    <database name> 
    <signal to noise ratio> 
    <bnx file template>
    ```

    BNX file template can be found in photomap folder as 'bnx_header.txt'.
    
    example usage:
    ```bash
    ./write_bnx /path/to/apple.db 3 ../bnx_header.txt
    ```