# OptiScan Tool

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
* [dash](https://github.com/plotly/dash) (optional for user interface)

## Installation

### Library
```bash
git clone https://gitlab.com/akdel/OptiScan.git
cd OptiScan
pip install .
chmod +x pipelines/extract_molecules
chmod +x pipelines/write_bnx
```

### Dashboard
1. Run the http server with:
    ```bash
    cd dashboard
    python optiscan_app.py localhost 8080 > log.txt &
    ```
2. Access from browser with http://localhost:8080/optiscan.

## Usage

### Running from OptiScan dashboard

OptiScan dashboard can be run as a web application built with [`Dash`](https://github.com/plotly/dash) python library. In this interface you can execute molecule detection/extraction and inspect the raw molecules and molecule length distributions. These inspections aid the choice of SNR, maximum DNA backbone intensity and minimum molecule length thresholds prior to exporting the data into `bnx` format.

#### Using dashboard with test data

You can find two identical scans in `dashboard/test_data/test_run/` directory in `test_scans.tar.gz`. After extracting the scans with `tar -xvf scans.tar.gz`, the path to `test_run` directory can be given in dashboard as the **Run Folders**.

Providing all the test parameters in **Molecule detection** section and running OptiScan is required to extract molecules and create the test database. Below are a list of parameters needed to detect molecules from the test data:

|: Parameter :|: Input :|
|:---------:|:-----|
|  Platform |  Irys   |
| Database Name | test.db |
| Num. of threads | 1 or 2 (should not be more than number of scans)|
| Run folders | dashboard/test_data/test_run |
| Chip dimensions | test dimensions |
| Organism | test |



#### Dashboard screenshots

![](screenshot1.png)
Molecule detection

![](screenshot2.png)
Molecule inspection

### Commandline scripts
`optiscan/pipelines/` folder contains two pipeline scripts for running OptiScan.


1. `extract_molecules` is for extracting and recording OM signals into disk.

    arguments: 
    ```bash
    <list of runs: use commas to separate> 
    <chip dimension: x,y> 
    <database name: str>
    <number of threads: int> 
    <organism name: str>
    <platform name: 'irys' or 'saphyr'>
    ```

    example usage:

    ```bash
    ./extract_molecules '/path/to/run1,/path/to/run/2' '12,95' apple.db 10 apple irys
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