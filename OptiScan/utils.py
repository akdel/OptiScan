from OptiScan import *
from OptiScan import numba_funcs
from scipy import signal
from scipy import ndimage
import numpy as np
from numba import njit, jit, vectorize, prange
from intervaltree import IntervalTree


CMAP_HEADER = """# CMAP File Version:    0.1
# Label Channels:       1
# Nickase Recognition Site 1:   %s
# Enzyme1:      %s
# Number of Consensus Nanomaps: 13
#h CMapId       ContigLength    NumSites        SiteID  LabelChannel    Position        StdDev  Coverage        Occurrence
#f int  float   int     int     int     float   float   int     int
"""

__author__ = "MAkdel"


"""
Functions for signal pairs. These should include bnx maker.
"""


DTYPE = np.dtype({'names': ["overlap_score", "long_id", "short_id", "long_start", "long_end",
                            "short_start", "short_end", "long_len", "short_len", "contained", "reversed"],
                  'formats': [np.float64, np.int64, np.int64, np.float64, np.float64, np.float64,
                              np.float64, np.float64, np.float64, np.bool, np.bool]})


def get_bnx_info(nick_signal: np.ndarray, snr: float) -> dict:
    bnx_dict = dict()
    median = np.median(nick_signal)
    if median < 1:
        median = 1
    bnx_dict["length"] = nick_signal.shape[0]
    bnx_dict["nick_distances"] = get_peaks(nick_signal, snr, median)
    bnx_dict["nick_snrs"] = [nick_signal[x] / median for x in bnx_dict["nick_distances"]]
    bnx_dict["nick_intensities"] = [nick_signal[x] for x in bnx_dict["nick_distances"]]
    bnx_dict["median"] = median
    return bnx_dict


def get_bnx_info_force_median(nick_signal: np.ndarray, snr: float, median: float) -> dict:
    bnx_dict = dict()
    bnx_dict["length"] = nick_signal.shape[0]
    bnx_dict["nick_distances"] = get_peaks(nick_signal, snr, median)
    bnx_dict["nick_snrs"] = [nick_signal[x] / median for x in bnx_dict["nick_distances"]]
    bnx_dict["nick_intensities"] = [nick_signal[x] for x in bnx_dict["nick_distances"]]
    bnx_dict["median"] = median
    return bnx_dict


def get_peaks(nick_signal: np.ndarray, snr: float, median: float) -> list:
    return [x for x in signal.argrelextrema(nick_signal, np.greater)[0] if nick_signal[x] >= median*snr]


# def get_peaks(nick_signal: np.ndarray, snr: float, median: float) -> list:
#     return list(signal.find_peaks(nick_signal, width=20, height=snr*median)[0])


def get_bnx_info_with_upsampling(nick_signal: np.ndarray, snr: float, zoom_ratio: float, final_ratio: float) -> dict:
    final_ratio /= zoom_ratio
    nick_signal = ndimage.zoom(nick_signal, zoom_ratio)
    bnx_info = get_bnx_info(nick_signal, snr)
    bnx_info["nick_distances"] = [x*final_ratio for x in bnx_info["nick_distances"]]
    bnx_info["length"] = bnx_info["length"] * final_ratio
    return bnx_info


def get_bnx_info_with_upsampling_with_forced_median(nick_signal: np.ndarray, snr: float, median: float,
                                                    zoom_ratio: float, final_ratio: float) -> dict:
    final_ratio /= zoom_ratio
    nick_signal = ndimage.zoom(nick_signal, zoom_ratio)
    bnx_info = get_bnx_info_force_median(nick_signal, snr, median)
    bnx_info["nick_distances"] = [x*final_ratio for x in bnx_info["nick_distances"]]
    return bnx_info


class MolStats:
    """
    Parse and perform stats on mol files in bionano pipeline.
    """
    def __init__(self, mol_file_destionation):
        self.lines = list()
        self.mol_dict = dict()
        self.filename = mol_file_destionation
        self.headers = list()
        self.labels = list()
        self.label_snr = list()
        self.sorted_mols = list()

    def get_sorted_mols(self):
        mols = [(int(self.mol_dict[x]["PixelCount"]), x) for x in self.mol_dict]
        self.sorted_mols = sorted(mols)

    def get_labels_from_bnx(self, bnx_file):
        bnx = BnxParser(bnx_file)
        bnx.read_bnx_file()
        self.labels = bnx.bnx_labels
        self.label_snr = bnx.bnx_label_snr

    def parse_mols(self):
        f = open(self.filename, 'r')
        self.lines = f.readlines()[3:]
        f.close()
        self.headers = [x.strip() for x in self.lines[0].split()]
        self.lines = self.lines[1:]
        mol_id = 1
        for line in self.lines:
            line = [x.strip() for x in line.split()]
            self.mol_dict[mol_id] = {}
            for i in range(len(self.headers)):
                self.mol_dict[mol_id][self.headers[i]] = line[i]
            mol_id +=1

    def show_length_histogram(self, bin_width=50, min_mol_len=100):
        lengths = [float(self.mol_dict[x]["Length"]) for x in self.mol_dict if float(self.mol_dict[x]["Length"]) >= min_mol_len]
        max_len = int(max(lengths)) + bin_width
        histogram = dict()
        [histogram.update({(r, r + bin_width): 0}) for r in range(0, max_len - bin_width, bin_width)]
        for mol_len in lengths:
            for interval in histogram.keys():
                if interval[0] <= mol_len < interval[1]:
                    histogram[interval] += 1
                    break
                else:
                    continue
        average = sum(lengths)
        current = int()
        i = int()
        while not current >= (average/2):
            current += lengths[i]
            i += 1
        print("N50: %s" % lengths[i+1])
        import matplotlib.pyplot as plt
        plt.hist(lengths, bins=50, range=(0, 1100))
        plt.show()

    def return_label_count_according_to_length(self, min_mol_len=400, snr=15):
        lengths = [float(self.mol_dict[x]["Length"]) for x in self.mol_dict]
        total_label_count = int()
        mol_count = int()
        if len(lengths) != len(self.label_snr):
            print("Molecule counts between bnx and mol files did not match.")
            return None
        for i in range(len(lengths)):
            if lengths[i] >= min_mol_len:
                total_label_count += len([x for x in self.label_snr[i] if float(x) >= snr])
                mol_count += 1
        print("number of molecules counted: %s" % mol_count)
        print("number of lables which fits criteria: %s" % total_label_count)
        return total_label_count


class BnxParser:
    """
    A basic parses for bnx files. Keeps the complete information from the file. Functions will be added as needed..
    """
    def __init__(self, bnx_file_destination):
        """
        Several arrays are initiated. these arrays will store the raw text lines of the bnx and also the parsed info
        from the lines.
        :param bnx_file_destination: Path to bnx file; str
        """
        self.bnx_filename = bnx_file_destination
        self.bnx_lines = list()
        self.bnx_arrays = list()
        self.bnx_labels = list()
        self.bnx_label_snr = list()
        self.bnx_headers = list()
        self.bnx_head_text = str()

    def read_bnx_file(self):
        """
        Parses the information on the BNX file, while filling the bnx containers. (see init)
        :return:
        """
        f = open(self.bnx_filename, 'r')
        self.bnx_lines = f.readlines()
        f.close()
        header_line = ""
        for line in self.bnx_lines:
            if line.startswith("#0h"): header_line = line; break
            if line.startswith("#"): self.bnx_head_text += line
        self.bnx_lines = [x for x in self.bnx_lines if not x.startswith('#') and x != "\n"]
        self.bnx_headers = [x.strip() for x in header_line.split()][1:]
        for r in range(0, len(self.bnx_lines), 4):
            self.bnx_label_snr.append(self.bnx_lines[r+2].split()[1:])
            self.bnx_labels.append(self.bnx_lines[r+1].split()[1:])
            self.bnx_arrays.append({"info": [x.strip() for x in self.bnx_lines[r].split("\t")],
                                    "labels": [float(x) for x in self.bnx_lines[r+1].split()[1:]],
                                    "label_snr": [float(x) for x in self.bnx_lines[r+2].split()[1:]],
                                    "raw_intensities": [float(x) for x in self.bnx_lines[r+3].split()[1:]],
                                    "signal": None})
        sort_guide = list()
        for i in range(len(self.bnx_arrays)):
            merged_id = int(self.bnx_arrays[i]["info"][1])
            sort_guide.append((merged_id, i))
        sort_guide = sorted(sort_guide)
        sorted_bnx = list()
        for entry in sort_guide:
            index = entry[1]
            sorted_bnx.append(self.bnx_arrays[index])
        self.bnx_arrays = sorted_bnx


    def write_arrays_as_bnx(self, outname):
        """
        writes the bnx_arrays into file
        ie// this is useful if bnx_arrays are modified.
        """
        f = open(outname, "w")
        lines = [str(self.bnx_head_text)]
        for arr in self.bnx_arrays:
            lines.append("\t".join(arr["info"]))
            lines.append("1\t" +  "\t".join([str(x) for x in arr["labels"]]))
            lines.append("QX11\t" +  "\t".join([str(x) for x in arr["label_snr"]]))
            lines.append("QX12\t" +  "\t".join([str(x) for x in arr["raw_intensities"]]))
        f.write("\n".join(lines))
        f.close()


    def get_molecule(self, mol_id, start_from=0):
        """
        :param mol_id:
        :param start_from:
        :return:
        """
        i = 0
        for entry in self.bnx_arrays[start_from:]:
            if int(entry["info"][1]) == int(mol_id):
                return int(entry["info"][6]), int(entry["info"][7]), start_from + i
            i += 1

    def get_molecule_complete(self, mol_id):
        """
        Returns the complete info dict for a molecule.
        :param mol_id: molecule id; str or int
        :return: The molecule info dictionary; dict
        """
        return [entry for entry in self.bnx_arrays if int(entry["info"][1]) == int(mol_id)][0]

    def get_all_molecule_signals(self, bp=560.):
        for mol in self.bnx_arrays:
            label_distances = [int(x/bp) for x in mol["labels"]]
            label_intensities = [int(x) for x in mol["label_snr"]]
            total_len = label_distances[-1]
            signal_array = np.zeros(total_len, dtype=int)
            for i in range(len(label_intensities)):
                signal_array[label_distances[i] - 1] = label_intensities[i]
            mol["signal"] = signal_array

    def bin_molecules_according_to_label_densities(self, number_of_bins):
        """
        Obtains label density per molecule and bins them.
        :param number_of_bins: number of bins to be obtained at the end of the process.
        :return:
        """
        molecule_densities = dict()
        for mol in self.bnx_arrays:
            mol_length = float(mol["info"][2])
            number_of_labels = len(mol["labels"])
            molecule_densities[mol_length/number_of_labels] = int(mol["info"][1])
        top_density = int(max(molecule_densities))
        least_density = int(min(molecule_densities))
        density_space = top_density - least_density
        single_bin_range = int(density_space / number_of_bins)
        ordered_densities = sorted(molecule_densities.keys(), reverse=False)
        bin_ranges = [(x, x + single_bin_range) for x in
                      range(least_density, top_density - single_bin_range, single_bin_range)]
        bins = [list() for x in bin_ranges]
        mol_counter = 0
        for i in range(len(bin_ranges)):
            for density in ordered_densities[mol_counter:]:
                if bin_ranges[i][0] < density <= bin_ranges[i][1]:
                    bins[i].append(molecule_densities[density])
                else:
                    break
                mol_counter += 1
        return bin_ranges, bins

    def _molecule_id_to_bnx_lines(self, molecule_id):
        line_multiplier = 4
        line_number = [x for x in range(len(self.bnx_arrays))
                       if int(self.bnx_arrays[x]["info"][1]) == molecule_id][0] * line_multiplier
        return "".join(self.bnx_lines[line_number:line_number+4])

    def _molecule_ids_to_bnx_lines(self, list_of_molecule_ids):
        list_of_lines = list()
        for mol_no in list_of_molecule_ids:
            line = self._molecule_id_to_bnx_lines(mol_no)
            list_of_lines.append(line)
        return "".join(list_of_lines)

    def write_bnx_for_list_of_molecules(self, list_of_molecule_ids, output_file_name):
        """
        Filters bnx file by given list of molecule ids.
        :param list_of_molecule_ids: List of bnx molecule ids
        :param output_file_name:
        :param bnx_head_file:
        :return:
        """
        from OptiScan import bnx_head
        bnx_head_text = str(bnx_head.HEADER)
        bnx_body = self._molecule_ids_to_bnx_lines(list_of_molecule_ids)
        output_file = open(output_file_name, "w")
        output_file.write(bnx_head_text + bnx_body)
        output_file.close()

    def filter_with_minimum_length(self, min_len):
        """
        Filters out molecules which are lower than a given minimum molecule length.
        :param min_len: minimum molecule length
        :return:
        """
        to_del = [x for x in range(len(self.bnx_arrays)) if float(self.bnx_arrays[x]["info"][2]) < min_len]
        for idx in to_del[::-1]:
            del self.bnx_arrays[idx]
            del self.bnx_lines[idx*4:idx*4+4]


class MqrMapParser:
    def __init__(self, map_file):
        self.map_filename = map_file
        self.lines = list()
        self.headers = list()
        self._get_lines_and_header()

    def _get_lines_and_header(self):
        map_file = open(self.map_filename, "r")
        all_lines = map_file.readlines()
        map_file.close()
        all_lines = [x for x in all_lines if not x.startswith("#")]
        self.headers = [x.strip() for x in all_lines[0].split('\t')]
        self.lines = [[y.strip() for y in x.split()] for x in all_lines[1:]]

    def get_total_map_length(self):
        pass


class CmapParser:
    def __init__(self, cmap_file_location, ini=False):
        self.cmap_file_location = cmap_file_location
        self.cmap_lines = list()
        self.cmap_headers = list()
        self.intervals = dict()  # tuples as keys and molecule ids as values
        self.position_index = dict() # cmap ids as keys and list of matching positions as lists
        self.contig_match_stats = dict()
        self.mqr = None
        self.xmap = None
        if ini:
            self._initiate2()

    def _initiate2(self):
        self.read_and_load_cmap_file()
        self.get_position_indexes()
        self.positions_to_intervals()

    def read_and_load_cmap_file(self):
        f = open(self.cmap_file_location, "r")
        self.cmap_lines = f.readlines()
        f.close()
        header_line = ""
        for line in self.cmap_lines:
            if line.startswith("#h"): header_line = line; break
        self.cmap_headers = [x.strip() for x in header_line.split()][1:]
        self.cmap_lines = [x.split() for x in self.cmap_lines if not x.startswith("#")]
        self.mqr = None

    def get_position_indexes(self):
        for line in self.cmap_lines:
            _id = int(line[0])
            _position = float(line[5])
            if _id not in self.position_index:
                self.position_index[_id] = [_position]
            else:
                self.position_index[_id].append(_position)

    def positions_to_intervals(self):
        for _id in self.position_index:
            start = self.position_index[_id][0]; end = self.position_index[_id][-1]
            interval = (start, end)
            if interval not in self.intervals:
                self.intervals[interval] = [_id]
            else:
                self.intervals[interval].append(_id)

    def parse_map_file_stats(self, mqr_path):
        self.mqr = MqrMapParser(mqr_path)
        for line in self.mqr.lines:
            contig_id = int(line[3])
            print(contig_id)
            if contig_id in self.contig_match_stats:
                self.contig_match_stats[contig_id].append((int(line[1]), int(line[9]), int(line[10])))
            else:
                self.contig_match_stats[contig_id] = [(int(line[1]), int(line[9]), int(line[10]))]

    def getN50(self):
        all_matches = list()
        if self.mqr:
            for line in self.mqr.lines:
                all_matches.append(int(line[10]) - int(line[9]))
        sorted_matches = sorted(all_matches)
        return (sorted_matches[-1], sorted(all_matches)[len(all_matches)/2], sorted_matches[0])

    def obtain_contig_coverage(self):
        results = dict()
        for contig_id in self.contig_match_stats:
            total_contig_length = float(self.position_index[contig_id][-1]) - float(self.position_index[contig_id][0])
            mat = np.zeros((int(self.position_index[contig_id][-1]), len(self.contig_match_stats)))
            mat[:,:] = np.NaN
            for i, stat_tuple in enumerate(self.contig_match_stats[contig_id]):
                mat[int(stat_tuple[1]):int(stat_tuple[2]),i] = 1
            mat = np.nanmedian(mat, axis=1)
            results[contig_id] = (total_contig_length, np.nansum(mat), (np.nansum(mat)/total_contig_length)*100,
                                  len(self.contig_match_stats[contig_id]))
        return results

    def map_and_show_all(self,z):
        import matplotlib.pyplot as plt
        ydim = len(self.intervals) * 7
        xdim = int(max([x[1] for x in self.intervals])/z) + 1
        all_om = np.zeros((ydim, xdim), dtype=np.bool)
        _id = 0
        intervals = sorted([(x[0]-x[1],x) for x in self.intervals.keys()])
        for interval in [x[1] for x in intervals]:
            array_index = _id*6
            print(_id)
            all_om[array_index:array_index + 4, int(interval[0]/z):int(interval[1]/z)] = True
            _id += 1
        plt.imsave("test.png", all_om)
        return all_om


class XmapParser:
    def __init__(self, xmap_file_location, ini=False):
        self.xmap_file_location = xmap_file_location
        self.xmap_lines = list()
        self.xmap_headers = list()
        self.intervals = dict()  # tuples as keys and molecule ids as values
        self.label_intervals = dict()
        self.position_index = dict() # cmap ids as keys and list of matching positions as lists
        self.contig_match_stats = dict()
        self.mqr = None
        self.last_contig_label = 0
        self.maps = None
        self.molecule_mathced_labels = dict()
        if ini:
            self._initiate2()
        self.moa_tree = IntervalTree()
        self.pairwise_overlaps = dict()
        self.total_overlaps = 0

    def _initiate2(self):
        self.read_and_load_xmap_file()
        self.get_intervals()
        self.obtain_molecule_match_info()

    def read_and_load_xmap_file(self):
        f = open(self.xmap_file_location, "r")
        self.xmap_lines = f.readlines()
        f.close()
        header_line = ""
        for line in self.xmap_lines:
            if line.startswith("#h"): header_line = line; break
        self.xmap_headers = [x.strip() for x in header_line.split()][1:]
        self.xmap_lines = [x.split() for x in self.xmap_lines if not x.startswith("#")]
        self.mqr = None

    def obtain_molecule_match_info(self):
        for line in self.xmap_lines:
            mol_id = int(line[1])
            ori = line[7]
            matches = list()
            for match in line[-1].split(')('):
                ids = tuple([int(x) for x in match.replace(')', "").replace("(", "").split(',')])
                matches.append(ids)
            max_label = max([x[0] for x in matches])
            min_label = min([x[0] for x in matches])
            if max_label > self.last_contig_label:
                self.last_contig_label = max_label
            self.molecule_mathced_labels[mol_id] = {"ori": ori, "matches": matches, "label_interval": (min_label,max_label)}
            self.moa_tree.addi(min_label, max_label, mol_id)

    def _mol_vs_tree(self, mol_id, space):
        mol_interval = self.molecule_mathced_labels[mol_id]["label_interval"]
        if mol_interval[-1] - mol_interval[0] < space*2:
            return None
        return self.moa_tree[mol_interval[0]-6:mol_interval[0]+6]

    def find_pairwise_overlaps(self, space=6):
        self.total_overlaps = 0
        for mol_id in self.molecule_mathced_labels:
            overlaps = self._mol_vs_tree(mol_id, space)
            self.pairwise_overlaps[mol_id] = overlaps
            if overlaps:
                self.total_overlaps += len(list(overlaps))


    def show_label_matrix(self, space=15, highlight=[]):
        empty = np.zeros((len(self.xmap_lines), self.last_contig_label))
        for i in range(len(list(self.molecule_mathced_labels.keys()))):
            mol_id = self.molecule_mathced_labels.keys()[i]
            interval = self.molecule_mathced_labels[mol_id]["label_interval"]
            if self._mol_vs_tree(mol_id, space):
                empty[i,interval[0]:interval[1]] = 1
            else:
                empty[i, interval[0]:interval[1]] = 0
            if mol_id in highlight:
                empty[i, interval[0]:interval[1]] = 2


        import matplotlib.pyplot as plt
        plt.imshow(empty)
        plt.show()


    def get_intervals(self):
        for line in self.xmap_lines:
            _id = int(line[1])
            start = float(line[5])
            end = float(line[6])
            self.intervals[(start, end)] = _id

    def parse_map_file_stats(self, mqr_path):
        self.mqr = MqrMapParser(mqr_path)
        for line in self.mqr.lines:
            contig_id = int(line[3])
            print(contig_id)
            if contig_id in self.contig_match_stats:
                self.contig_match_stats[contig_id].append((int(line[1]), int(line[9]), int(line[10])))
            else:
                self.contig_match_stats[contig_id] = [(int(line[1]), int(line[9]), int(line[10]))]

    def getN50(self):
        all_matches = list()
        if self.mqr:
            for line in self.mqr.lines:
                all_matches.append(int(line[10]) - int(line[9]))
        sorted_matches = sorted(all_matches)
        return (sorted_matches[-1], sorted(all_matches)[len(all_matches)/2], sorted_matches[0])

    def obtain_contig_coverage(self):
        results = dict()
        for contig_id in self.contig_match_stats:
            total_contig_length = float(self.position_index[contig_id][-1]) - float(self.position_index[contig_id][0])
            mat = np.zeros((int(self.position_index[contig_id][-1]), len(self.contig_match_stats)))
            mat[:,:] = np.NaN
            for i, stat_tuple in enumerate(self.contig_match_stats[contig_id]):
                mat[int(stat_tuple[1]):int(stat_tuple[2]),i] = 1
            mat = np.nanmedian(mat, axis=1)
            results[contig_id] = (total_contig_length, np.nansum(mat), (np.nansum(mat)/total_contig_length)*100,
                                  len(self.contig_match_stats[contig_id]))
        return results

    def map_and_show_all(self,z):
        import matplotlib.pyplot as plt
        ydim = len(self.intervals) * 7
        xdim = int(max([x[1] for x in self.intervals])/z) + 1
        all_om = np.zeros((ydim, xdim), dtype=np.bool)
        _id = 0
        intervals = sorted([(x[0]-x[1],x) for x in self.intervals.keys()])
        for interval in [x[1] for x in intervals]:
            array_index = _id*6
            print(_id)
            all_om[array_index:array_index + 4, int(interval[0]/z):int(interval[1]/z)] = True
            _id += 1
        plt.imsave("test.png", all_om)
        return all_om

    def xmap_line_details(self, xmap_line:list):
        res = {'chr_no': None, 'chr_start': None, 'chr_end': None, 'ori': None, 'mol_id': None}
        res['ori'] = xmap_line[-7]
        res['chr_no'] = int(xmap_line[2]) - 1
        res['mol_id'] = int(xmap_line[1]) - 1
        if res['ori'] == '+':
            res['chr_start'] = int(float(xmap_line[5])) - int(float(xmap_line[3]))
            res['chr_end'] = res['chr_start'] + int(float(xmap_line[-4]))
        elif res['ori'] == '-':
            res['chr_start'] = int(float(xmap_line[5])) - (int(float(xmap_line[-4])) - int(float(xmap_line[3])))
            res['chr_end'] = res['chr_start'] + int(float(xmap_line[-4]))
        return res


class LinkXmapToBnx(BnxParser, XmapParser):
    def __init__(self, bnx_file_path, xmap_file_path):
        BnxParser.__init__(self, bnx_file_path)
        XmapParser.__init__(self, xmap_file_path)
        self.mapped_mols = None
        self._initiate2()
        self.read_bnx_file()
        self._initiate_array()
        self.ref_starts = list()

    def _initiate_array(self):
        y_dim = len(self.xmap_lines)
        x_dim = int(sorted([x[1] for x in self.intervals.keys()])[-1])
        self.mapped_mols = np.zeros((y_dim, x_dim), dtype=float)
        self.mapped_mols[:] = np.nan

    def place_molecules(self):
        mol_rank = 0
        for mol in self.xmap_lines:
            mol_id = int(mol[1])
            ref_start = int(float(mol[5]))
            # print ref_start
            q_start = float(mol[3])
            q_end = float(mol[4])
            if int(q_start) > int(q_end):
                flip = True
            else:
                flip = False
            bnx_mol = self.get_molecule_complete(mol_id)
            bnx_start = bnx_mol["labels"].index(q_start)
            bnx_end = bnx_mol["labels"].index(q_end)
            if flip:
                bnx_labels = bnx_mol["labels"][bnx_end:bnx_start]
                bnx_snr = bnx_mol["label_snr"][bnx_end:bnx_start][::-1]
                bnx_labels = [int(bnx_labels[-1] - x) for x in bnx_labels[::-1]]
            else:
                bnx_labels = bnx_mol["labels"][bnx_start:bnx_end]
                bnx_snr = bnx_mol["label_snr"][bnx_start:bnx_end]
                bnx_labels = [int(x - bnx_labels[0]) for x in bnx_labels]
            ref_map_index = [x + ref_start for x in bnx_labels]
            self.ref_starts.append(ref_map_index[0])
            print("\t".join([str(x) for x in ref_map_index]))
            print("\t".join([str(x) for x in bnx_snr]))
            print("\n")
            self.mapped_mols[mol_rank, ref_map_index[0]:ref_map_index[-1]] = 1
            for i in range(len(ref_map_index)):
                self.mapped_mols[mol_rank, ref_map_index[i] - 5:ref_map_index[i] + 5] = bnx_snr[i]
            mol_rank += 1


class FastaObject:
    """
    Fasta file is transformed into numpy arrays and numba functions are used on it to compute string matching in parallel.
    """
    def __init__(self, fasta_file_path: str):
        """
        Opens the fasta/multifasta file and creates a dictionary from it's content. Keys being the fasta title and
        values being the corresponding sequence.
        :param fasta_file_path: filte path to fasta file
        :return:

        NOTE example usage:
        >>> digestion_sequence = "ACCTGACCA"
        >>> fasta = FastaObject("/path/to/file.fasta")
        ... fasta.fill_complete_fasta_array() # fills in array with sequence
        ... fasta.digest_array(digestion_sequence) # in-silico digestion of the sequence
        ... fasta.clean_fasta_dict() # removes the arbitary sequence from memory
        >>> np.where(np.sum(fasta.fasta_digestion_array, axis=0) > 0) # This returns the digestion indices
        """
        self.file_path = fasta_file_path
        self.fasta = dict()
        self.pattern = str()
        self.info_txt = 'This is a fasta sequence object.'
        self.contig_order_in_fasta = list()
        self.lengths = list()
        if fasta_file_path:
            fasta_file = open(self.file_path, 'r')
            fasta_lines = fasta_file.readlines()
            fasta_file.close()
        else:
            fasta_lines = None

        if fasta_lines:
            for line in fasta_lines:
                if line.startswith('>'):
                    current_header = line[1:].strip()
                    self.contig_order_in_fasta.append(current_header)
                    self.fasta[current_header] = []
                else:
                    current_sequence = line.strip()
                    try:
                        self.fasta[current_header] += [current_sequence]
                    except KeyError:
                        print("Incorrect file format. (file not fasta)")
            for entry in self.fasta:
                self.fasta[entry] = "".join(self.fasta[entry])
            self._get_lengths()
        self.fasta_array = None
        self.fasta_digestion_array = None

    def initiate_fasta_array(self):
        """
        Initiates the fasta array in both directions.
        """
        self.fasta_array = np.zeros((2, np.sum(self.lengths)), dtype=np.dtype("|S1"))

    def fill_fasta_chrom(self, chromosome_name, chromosome_index):
        """
        Fills in a chromosome within the array with the corresponding sequence.
        """
        string = self.fasta[chromosome_name]
        numba_funcs.fill_chromosomes_with_fasta(string, chromosome_index, self.fasta_array[0])

    def fill_reverse_complement(self):
        """
        Fills in the reverse complement.
        """
        self.fasta_array[1] = numba_funcs.reverse_complement_v2(self.fasta_array[0].view("int8")).view("S1")

    def clean_fasta_dict(self):
        """
        Deletes the fasta sequence to save memory.
        """
        del self.fasta

    def fill_complete_fasta_array(self):
        """
        Fills in the fasta array with sequences.
        """
        for i in range(len(self.contig_order_in_fasta)):
            chromosome_name = self.contig_order_in_fasta[i]
            print(chromosome_name)
            chromosome_index = np.sum(self.lengths[:i])
            self.fill_fasta_chrom(chromosome_name, chromosome_index)
        self.fill_reverse_complement()

    def _get_lengths(self):
        """
        Obtains the chromosome lengths.
        """
        for header in self.contig_order_in_fasta:
            contig_sequence_length = len(self.fasta[header])
            self.lengths.append(contig_sequence_length)

    def digest_fasta_array(self, digestion_sequence:str):
        """
        Digests the fasta array with a given digestion sequence.
        """
        self.fasta_digestion_array = np.zeros((2, self.fasta_array[0].shape[0]), dtype=np.bool)
        digestion_array = np.array(tuple(digestion_sequence), dtype="|S1").view("int8")
        self.fasta_digestion_array[0] = numba_funcs.digest_fasta(self.fasta_array[0].view("int8"), digestion_array, self.fasta_digestion_array[0])
        self.fasta_digestion_array[1] = numba_funcs.digest_fasta(self.fasta_array[1].view("int8"), digestion_array, self.fasta_digestion_array[1])
        self.fasta_digestion_array[1] = self.fasta_digestion_array[1][::-1]

    def write_fasta_to_cmap(self, digestion_sequence:str, output_file_name:str, enzyme_name="BSP1Q", channel=1):
        """
        Writes the sequences into a CMAP file with an in silico digestion by the given digestion sequence.
        """
        def gen_line():
            return '%s\t%s\t%s\t%s\t%s\t%s\t1.0\t1\t1\n'

        self.digest_fasta_array(digestion_sequence)
        digested_array = np.sum(self.fasta_digestion_array , axis=0)

        f = open(output_file_name, 'w')
        f.write(CMAP_HEADER % (digestion_sequence, enzyme_name))
        
        for i in range(len(self.lengths)-1):
            start = sum(self.lengths[:i])
            end = sum(self.lengths[:i+1])
            _id = i + 1
            nicking_sites = np.where(digested_array[start:end] > 0.0)[0]
            if nicking_sites.shape[0] == 0:
                continue
            else:
                length = end-start
                for j in range(nicking_sites.shape[0]):
                    line = gen_line() % (_id, length, nicking_sites.shape[0], j+1, channel, nicking_sites[j])
                    f.write(line)
                line = gen_line()
                line = line % (_id, length, nicking_sites.shape[0], nicking_sites.shape[0]+1, 0, length)
                f.write(line)
        f.close()


class MQR:
    """
    A simple wrapper for aligning molecules by RefAligner.
    """
    def __init__(self, output_folder, ref_align_path, score="1e-10"):
        """
        Sets and stores the output folder and refaligner part. Also sets the BNG command.

        output_folder: Output folder path where refaligner will use.
        ref_align_path: Path to the refaligner software.
        score: Score for filtering results. default is 1e-10

        NOTE usage example of this class:
        >>> mqr = MQR("/path/to/folder/", "/path/to/RefAlign", score="1e-11")
        ... mqr.run_ref_align("reference.cmap", "molecules.bnx", 10000)
        ... mqr.load_results()
        ... xmap = mqr.xmap                     # This is an XmapParser instance which can 
        ... xmap.read_and_load_xmap_file()      # be used to work on the resulting xmap file.
        ... print(xmap.xmap_lines[0])    # This prints the first line of results.
        """
        self.score = score
        self.ref_align = ref_align_path
        self.output_dir = output_folder
        self.command = """%s -f -ref %s -i %s -o %s -nosplit 2 -BestRef 1 -biaswt 0 -Mfast 0 -FP 1.5 -FN 0.05 \
            -sf 0.2 -sd 0.0 -A 5 -outlier 1e-3 -outlierMax 40 -endoutlier 1e-4 -S -1000 -sr 0.03 -se 0.2 -MaxSF \
             0.3 -MaxSE 0.5 -MaxSD 0.12 -resbias 4 64 -maxmem 64 -M 3 3 -minlen 50  -T %s -maxthreads 12 \
             -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 3 -hash -hashdelta 10 -hashoffset 1 -hashmaxmem 64 -insertThreads 4 \
             -maptype 0 -PVres 2 -PVendoutlier -AlignRes 2.0 -rres 0.9 -resEstimate -ScanScaling 2 -RepeatMask 5 0.01 \
             -RepeatRec 0.7 0.6 1.4 -maxEnd 50 -usecolor %s -stdout -stderr -subset 1 %s"""
        self.pairwise_command = """%s -i %s -o %s -usecolor 1 -FP 1.5 -FN 0.15 -sd 0.0 -sf 0.2 -sr 0.03 -res 3.3 -T %s -maxmem 7.5 \
            -minlen 150 -minsites 8 -MaxIntensity 20000 -usecolor 1 -maxsites 200 -mres 0.9 -usecolor 1 -A 5 \
            -S 1 -MaxSE 0.5 -outlier 0.0001 -outlierMax 40. -endoutlier 0 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 1.4 -PVres 2 \
            -alignscore -maptype 0 -HSDrange 1.0 -hashoffset 1 -f -hashgen 5 3 2.2 1.2 0.05 3.0 1 1 1 -hash \
            -nosplit 2 -align_format 0 -stdout -stderr -maxthreads 24"""
        self.denovo_command = """python2 %s -l %s -t %s -b %s -R -y 1 -i 1 -a %s"""
        self.xmap = None

    def run_ref_align(self, reference_cmap_path, bnx_file, number_of_molecules, color=1):
        from subprocess import check_call as ck
        ck(self.command % (self.ref_align, reference_cmap_path, bnx_file, self.output_dir + "bnx_quality", self.score,
                           color, number_of_molecules), shell=True)

    def run_ref_align_pairwise(self, bnx_file):
        from subprocess import check_call as ck
        ck(self.pairwise_command % (self.ref_align, bnx_file, self.output_dir + "bnx_quality", self.score), shell=True)

    def rough_denovo_assembly(self, pipeline_location, output_folder, refaligner_directory, bnx_file, xml_file):
        from subprocess import check_call as ck
        ck(self.denovo_command % (pipeline_location, output_folder, refaligner_directory, bnx_file, xml_file))

    def load_results(self):
        xmap_file = self.output_dir + 'bnx_quality.xmap'
        self.xmap = XmapParser(xmap_file)

class MoleculeToMolecule(MQR):
    def __init__(self, output_folder, ref_align_path, score="1e-10"):
        MQR.__init__(self, output_folder=output_folder, ref_align_path=ref_align_path, score=score)
        self.some_to_all = """%s -i %s -i %s -o %s -usecolor 1 -FP 1.5 -FN 0.15 -sd 0.0 -sf 0.2 -sr 0.03 -res 3.3 -T %s -maxmem 7.5 \
            -minlen 150 -minsites 8 -MaxIntensity 20000 -usecolor 1 -maxsites 200 -mres 0.9 -usecolor 1 -A 5 \
            -S 1 -MaxSE 0.5 -outlier 0.0001 -outlierMax 40. -endoutlier 0 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 1.4 -PVres 2 \
            -alignscore -maptype 0 -HSDrange 1.0 -hashoffset 1 -f -hashgen 5 3 2.2 1.2 0.05 3.0 1 1 1 -hash \
            -nosplit 2 -align_format 0 -stdout -stderr -maxthreads 24"""
    
    def run_some_vs_all(self, bnx1, bnx2):
        from subprocess import check_call as ck
        ck(self.some_to_all % (self.ref_align, bnx1, bnx2, self.output_dir + "bnx_quality", self.score), shell=True)
    

class AlignParser:
    def __init__(self, alignment_file):
        self.align_lines = list()
        self.align_headers = list()
        f = open(alignment_file, "r")
        for l in f.readlines():
            if not l.startswith("#") and len(l) != 0:
                self.align_lines.append(l.split("\t"))
            elif l.startswith("#>0"):
                self.align_headers += l.split("\t")
        f.close()
        self.align_info = list()
        for i in range(0, len(self.align_lines), 5):
            current_info = dict()
            match_info = dict()
            for j in range(len(self.align_headers)):
                match_info[self.align_headers[j]] = self.align_lines[i][j]
            current_info["match_info"] = match_info
            current_info["matches_A"] = [int(x.strip()) for x in self.align_lines[i+1] if x.strip()]
            current_info["matches_B"] = [int(x.strip()) for x in self.align_lines[i+2] if x.strip()]
            current_info["scores_A"] = [float(x.strip()) for x in self.align_lines[i+3] if x.strip()]
            current_info["scores_B"] = [float(x.strip()) for x in self.align_lines[i+4] if x.strip()]
            self.align_info.append(current_info)
    
    def get_matched_pairs(self):
        res = set()
        for entry in self.align_info:
            res.add(tuple(sorted([int(entry["match_info"]["Mol0ID"])-1,
                                  int(entry["match_info"]["Mol1ID"])-1])))
        return res

    def get_all_ids(self):
        res = set()
        for entry in self.align_info:
            res.add(int(entry["match_info"]["Mol0ID"])-1)
            res.add(int(entry["match_info"]["Mol1ID"])-1)
        return res

    def refaligner_to_graph_edges(self, thr, mshift=536798):
        # mshift = max(list(self.get_all_ids()))
        template = "%s\t%s\t%s\n"
        edge_list = list()
        for i in range(len(self.align_info)):
            current_info = self.align_info[i]["match_info"]
            id1, id2, reverse, shift, score = current_info["Mol0ID"],current_info["Mol1ID"],current_info["Orientation"],current_info["Offset\n"],float(current_info["PvalueLog10"])
            id1 = int(id1)-1; id2 = int(id2)-1
            if (id1 > 536798) or (id2 > 536798):
                continue
            if score < thr:
                continue
            else:
                if int(reverse) == 1:
                    reverse = False
                else:
                    reverse = True
                shift = float(shift.strip())
                overhang = shift * 2
                if overhang > 0 and not reverse:
                    l = template % (id1, id2, overhang)
                    l2 = template % (id2+mshift, id1+mshift, overhang)
                elif overhang > 0 and reverse:
                    l = template % (id1+mshift, id2, overhang)
                    l2 = template % (id2+mshift, id1, overhang)
                elif overhang < 0 and not reverse:
                    l = template % (id2, id1, abs(overhang))
                    l2 = template % (id1+mshift, id2+mshift, abs(overhang))
                else:
                    l = template % (id2, id1+mshift, abs(overhang))
                    l2 = template % (id1+mshift, id2, abs(overhang))
                edge_list.append(l)
        return edge_list

    def combine_edgelist_with_current_file(self, fname, thr=15):
        edge_list = self.refaligner_to_graph_edges(thr)
        f = open(fname, "a")
        f.write("".join(edge_list))
        f.close()
