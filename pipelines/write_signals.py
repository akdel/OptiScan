#!/usr/bin/env python

from sys import argv
import numpy as np

def create_bnx_file():
    for run_id in connection.db_runs:
        mols = connection.yield_molecule_signals_in_run(run_id)
        np.save(str(run_id)+"OptiScan_nick.npy", np.array([x[0] for x in mols], dtype="O"))
        mols = connection.yield_molecule_signals_in_run(run_id)
        np.save(str(run_id)+"OptiScan_back.npy", np.array([x[1] for x in mols], dtype="O"))
        #molecules_to_bnxv2(mols, 10., 485., str(run_id) + "OptiScan.bnx", bnx_template_path=bnx_file_template,
        #                 signal_to_noise_ratio=snr)


def get_number_of_molecules_in_run(run):
    n = 0
    for scan_id in connection.db_runs[run]:
        scan_info = connection.db_runs[run][scan_id]
        for col_id in scan_info["columns"]:
            col_info = scan_info["columns"][col_id]
            n += int(col_info["molecule_count"])
    return n


if __name__ == "__main__":
    if argv[1] == "-h" or argv[1] == "--help" or argv[1] == "help":
        print("parameters: <database name> <signal to noise ratio> <bnx file template>\n\n"
              "bnx file template can be found in OptiScan folder as 'bnx_header.txt'.\n"
              "\nExample Usage:\n"
              "write_bnx /path/to/apple.db 3 ../bnx_head.txt\n\n")
    else:
        try:
            assert len(argv) == 4
        except AssertionError:
            print("All parameters should be provided. use -h or --help for help.")
            exit(3)
        from OptiScan.database import MoleculeConnector, molecules_to_bnxv2
        from sys import argv
        database_name = argv[1]
        snr = float(argv[2])
        bnx_file_template = argv[3]
        connection = MoleculeConnector(database_name)
        connection.connect_db()
        print("Run folders in database:")
        for run_id in connection.db_runs:
            print(run_id)
            print("\tNumber of scans in run -> %s" % len(connection.db_runs[run_id]))
            print("\tNumber of molecules in run -> %s" % get_number_of_molecules_in_run(run_id))
        create_bnx_file()
    exit()