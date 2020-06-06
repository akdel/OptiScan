#!/usr/bin/env python

from sys import argv
import numpy as np
from scipy import ndimage


def create_bnx_file():
    molecule_labels = list()
    molecule_backbones = list()
    second_labels = list()
    for bank_id in run_tree.bank_info:
        bank_path = run_tree.saphyr_tree.root + bank_id
        if dual:
            number_of_columns = len(run_tree.bank_info[bank_id])//3
        else:
            number_of_columns = len(run_tree.bank_info[bank_id]) // 2
        for i in range(number_of_columns):
            try:
                molecule_labels += list(np.load(f"{bank_path}/{i}_label.npy", allow_pickle=True))
                molecule_backbones += list(np.load(f"{bank_path}/{i}_backbone.npy", allow_pickle=True))
                if dual:
                    second_labels += list(np.load(f"{bank_path}/{i}_label2.npy", allow_pickle=True))
            except FileNotFoundError:
                continue
    if not dual:
        mols = [ndimage.white_tophat(molecule_labels[x], structure=np.ones(10)) for x in range(len(molecule_labels))]

        np.save("labels.npy", mols)
        np.save("backbones.npy", molecule_backbones)
    else:
        mols1 = [ndimage.white_tophat(molecule_labels[x], structure=np.ones(10)) for x in range(len(molecule_labels))]
        mols2 = [ndimage.white_tophat(second_labels[x], structure=np.ones(10)) for x in range(len(second_labels))]
        np.save("labels.npy", mols1)
        np.save("labels2.npy", mols2)
        np.save("backbones.npy", molecule_backbones)

if __name__ == "__main__":
    if argv[1] == "-h" or argv[1] == "--help" or argv[1] == "help": # parameters will be the run path and snr
        print("parameters: <run path> <signal to noise ratio> <bnx file template>\n\n"
              "bnx file template can be found in OptiScan folder as 'bnx_header.txt'.\n"
              "\nExample Usage:\n"
              "write_bnx /path/to/apple.db 3 ../bnx_head.txt\n\n")
    else:
        try:
            assert len(argv) == 5
        except AssertionError:
            print("All parameters should be provided. use -h or --help for help.")
            exit(3)
        from OptiScan.database import molecules_to_bnxv2
        from OptiScan.scan import SaphyrFileTree, SaphyrExtract
        from sys import argv
        if argv[4] == "True":
            dual = True
        else:
            dual = False
        run_path = argv[1]
        snr = float(argv[2])
        bnx_file_template = argv[3]
        run_tree = SaphyrExtract(run_path)
        create_bnx_file()
    exit()
