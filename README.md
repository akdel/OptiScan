# OptiScan Tool

OptiScan is an open source tool for extraction of nick and backbone signals
from provided Bionano Genomics optical map TIFF images. OptiScan is able to
transform nick signals into nicking site coordinates and output these in
the BNX format for further processing.


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
pip3 install . --user
chmod +x pipelines/extract_molecules
chmod +x pipelines/write_bnx
```

### Dashboard
1. Run the http server with:
    ```bash
    cd dashboard
    python3 optiscan_app.py localhost 8080 > log.txt &
    ```
2. Access it from a browser at http://localhost:8080/optiscan.

## Usage

### Running from OptiScan dashboard

While ordinarily OptiScan is best used from the command line, a dashboard is available. 
 It is built as a web application with 
[`Dash`](https://github.com/plotly/dash).  In this interface, you can
execute molecule detection/extraction and inspect the raw molecules and
molecule length distributions.  Inspecting these values supports the choice
of SNR, maximum DNA backbone intensity and minimum molecule length
thresholds prior to exporting the data in the BNX format.

#### Using dashboard with test data

You can find two scans in the directory `dashboard/test_data/test_run/`
in `test_scans.tar.gz`.  After extracting the scans with `tar -xvf
test_scans.tar.gz`, the path to the `test_run` directory can be specified in
the dashboard under **Run Folders**.

Providing all the test parameters in **Molecule detection** section and
running OptiScan is required to extract molecules and create the test
database.  Test parameters are already filled (except the database name) when
the dashboard is launched.  Molecules in the test data can be simply
detected by entering a database name and clicking **RUN OPTISCAN**.  This
process can take up to 5-10 minutes to complete, depending on your 
workstation speed.

You can find the live demonstration of OptiScan dashboard at [http://bioinformatics.nl/optiscan/](http://bioinformatics.nl/optiscan/). (note: the demo interface lacks the molecule detection function.)

#### Dashboard screenshots

![](screenshot1.png)
Molecule detection

![](screenshot2.png)
Molecule inspection

### Command line scripts
The `optiscan/pipelines/` folder contains two pipeline scripts for running OptiScan.

1. `extract_molecules` is for extracting and recording optical map signals.

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

2. `write_bnx` is for writing extracted molecules in the BNX format.
    
    arguments:
    ```bash 
    <database name> 
    <signal to noise ratio> 
    <bnx file template>
    ```

    The BNX file template can be found in the OptiScan directory as 'bnx_head.txt'.
    
    example usage:
    ```bash
    ./write_bnx /path/to/apple.db 3 ../bnx_head.txt
    ```