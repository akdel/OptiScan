#!/usr/bin/env python

from OptiScan.scan import SaphyrExtract, SaphyrFileTree
import multiprocessing as mp
from sys import argv
from os import path


def run_for_extraction(run_scan_pairs: [tuple]):
    run = SaphyrExtract(run_folder)
    print(True)
    run.save_molecules_for_bank_ids(run_scan_pairs, abstraction_threshold=contrast)


def split_banks_to_threads(bank_ids: str, number_of_threads: int):
    scan_jobs = bank_ids
    assert len(scan_jobs) >= number_of_threads
    number_of_jobs_in_thread = int(len(scan_jobs)/number_of_threads)
    left_overs = int(len(scan_jobs) % number_of_threads)
    if number_of_jobs_in_thread < 1:
        return [[x] for x in scan_jobs]
    else:
        current = [scan_jobs[i:i+number_of_jobs_in_thread] for i in range(0, len(scan_jobs) - left_overs,
                                                                          number_of_jobs_in_thread)]
        [current[i].append(scan_jobs[-i]) for i in range(1, left_overs + 1)]
        return current


if __name__ == "__main__":
    if argv[1] == "-h" or argv[1] == "--help" or argv[1] == "help": # Input will be the run folder and number of threads..
        print("arguments: <run folder>  <number of threads> <contrast threshold>")
    else:
        try:
            assert len(argv) == 4
        except AssertionError:
            print("All parameters should be provided. use -h or --help for help.")
            exit()
        run_folder = argv[1]
        assert path.isdir(run_folder)
        number_of_threads = int(argv[2])
        print(number_of_threads)
        contrast = int(argv[3])
        run_tree = SaphyrFileTree.from_run_folder(run_folder)
        scan_jobs = split_banks_to_threads(run_tree.get_bank_level_ids(), number_of_threads)
        number_of_jobs = len(scan_jobs)

        with mp.Pool(number_of_jobs) as job:
            job.map(run_for_extraction, scan_jobs)
    exit()

