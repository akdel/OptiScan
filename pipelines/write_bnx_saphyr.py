#!/usr/bin/env python

from sys import argv
import numpy as np
from scipy import ndimage


def create_bnx_file():
    for bank_id in run_tree.bank_info:
        bank_path = run_tree.saphyr_tree.root + bank_id
        number_of_columns = len(run_tree.bank_info[bank_id])//2
        molecule_labels = list()
        molecule_backbones = list()
        for i in range(number_of_columns):
            molecule_labels += list(np.load(f"{bank_path}/{i}_label.npy", allow_pickle=True))
            molecule_backbones += list(np.load(f"{bank_path}/{i}_backbone.npy", allow_pickle=True))

    mols = [(ndimage.white_tophat(molecule_labels[x], structure=np.ones(20)),
             molecule_backbones[x]) for x in range(len(molecule_labels))]
    molecules_to_bnxv2(mols, 10., 485., str(run_tree.saphyr_tree.root) + "OptiScan.bnx", bnx_template_path=bnx_file_template,
                     signal_to_noise_ratio=snr)


if __name__ == "__main__":
    if argv[1] == "-h" or argv[1] == "--help" or argv[1] == "help": # parameters will be the run path and snr
        print("parameters: <run path> <signal to noise ratio> <bnx file template>\n\n"
              "bnx file template can be found in OptiScan folder as 'bnx_header.txt'.\n"
              "\nExample Usage:\n"
              "write_bnx /path/to/apple.db 3 ../bnx_head.txt\n\n")
    else:
        try:
            assert len(argv) == 4
        except AssertionError:
            print("All parameters should be provided. use -h or --help for help.")
            exit(3)
        from OptiScan.database import molecules_to_bnxv2
        from OptiScan.scan import SaphyrFileTree, SaphyrExtract
        from sys import argv
        run_path = argv[1]
        snr = float(argv[2])
        bnx_file_template = argv[3]
        run_tree = SaphyrExtract(run_path)
        create_bnx_file()
    exit()
