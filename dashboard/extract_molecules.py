#!/usr/bin/env python

from OptiScan.database import MoleculeDB
import multiprocessing as mp
from sys import argv
from os import path


def run_for_extraction(run_scan_pairs: [tuple]):
    runs = list(set([run_id for run_id, _ in run_scan_pairs]))
    moldb_ = MoleculeDB(database_name, runs, organism, chip_dimension=chip_dimension)
    for run_id, scan_id in run_scan_pairs:
        moldb_.connect_db()
        moldb_.extract_molecules_from_scan_in_run(run_id, scan_id)
        moldb_.insert_columns_for_scan(run_id, scan_id)


def create_tables():
    moldb_ = MoleculeDB(database_name, list(set(list_of_run_folders)), organism, chip_dimension=chip_dimension)
    moldb_.create_tables_first_time()


if __name__ == "__main__":
    if argv[1] == "-h" or argv[1] == "--help" or argv[1] == "help":
        print("arguments: <list of runs: use commas to separate> <chip dimension: x,y> <database name: str> "
              "<number of threads: int> <organism name: str> <data type: 'irys' or 'saphyr'>\n\n"
              "Example command: extract_molecules '/path/to/run1,/path/to/run/2' '12,95' apple.db 10 apple")
    else:
        try:
            assert len(argv) == 7
        except AssertionError:
            print("All parameters should be provided. use -h or --help for help.")
            exit()
        if argv[6] == "saphyr":
            saphyr_ = True
        elif argv[6] == "irys":
            saphyr_ = False
        else:
            raise TypeError("data type should be set to irys or saphyr")
        list_of_run_folders = argv[1].split(",")
        for folder in list_of_run_folders:
            assert path.isdir(folder)
        chip_dimension = tuple([int(x.strip()) for x in argv[2].split(',')])
        assert len(chip_dimension) == 2
        database_name = argv[3]
        number_of_threads = int(argv[4])
        organism = argv[5]

        create_tables()
        print("Database initiated!")

        moldb = MoleculeDB(database_name, list_of_run_folders, organism, chip_dimension=chip_dimension, saphyr=saphyr_)
        scan_jobs = moldb.split_scans_to_threads(number_of_threads)
        number_of_jobs = len(scan_jobs)

        with mp.Pool(number_of_jobs) as job:
            job.map(run_for_extraction, scan_jobs)
    exit()

