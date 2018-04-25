from OptiScan import *
from scipy import signal
from scipy import ndimage
__author__ = "MAkdel"



"""
Functions for signal pairs. These should include bnx maker.
"""


def get_bnx_info(nick_signal: np.ndarray, snr: float) -> dict:
    bnx_dict = dict()
    median = np.median(nick_signal)
    if median < 1:
        median = 1
    bnx_dict["nick_distances"] = get_peaks(nick_signal, snr, median)
    bnx_dict["nick_snrs"] = [nick_signal[x] / median for x in bnx_dict["nick_distances"]]
    bnx_dict["nick_intensities"] = [nick_signal[x] for x in bnx_dict["nick_distances"]]
    bnx_dict["median"] = median
    return bnx_dict


def get_bnx_info_force_median(nick_signal: np.ndarray, snr: float, median: float) -> dict:
    bnx_dict = dict()
    bnx_dict["nick_distances"] = get_peaks(nick_signal, snr, median)
    bnx_dict["nick_snrs"] = [nick_signal[x] / median for x in bnx_dict["nick_distances"]]
    bnx_dict["nick_intensities"] = [nick_signal[x] for x in bnx_dict["nick_distances"]]
    bnx_dict["median"] = median
    return bnx_dict


def get_peaks(nick_signal: np.ndarray, snr: float, median: float) -> list:
    return [x for x in signal.argrelextrema(nick_signal, np.greater)[0] if nick_signal[x] >= median*snr]


def get_bnx_info_with_upsampling(nick_signal: np.ndarray, snr: float, zoom_ratio: float, final_ratio: float) -> dict:
    final_ratio /= zoom_ratio
    nick_signal = ndimage.zoom(nick_signal, zoom_ratio)
    bnx_info = get_bnx_info(nick_signal, snr)
    bnx_info["nick_distances"] = [x*final_ratio for x in bnx_info["nick_distances"]]
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
        self.bnx_lines = [x for x in self.bnx_lines if not x.startswith('#')]
        self.bnx_headers = [x.strip() for x in header_line.split()][1:]
        for r in range(0, len(self.bnx_lines), 4):
            self.bnx_label_snr.append(self.bnx_lines[r+2].split()[1:])
            self.bnx_labels.append(self.bnx_lines[r+1].split()[1:])
            self.bnx_arrays.append({"info": [x.strip() for x in self.bnx_lines[r].split("\t")],
                                    "labels": [float(x) for x in self.bnx_lines[r+1].split()[1:]],
                                    "label_snr": [float(x) for x in self.bnx_lines[r+2].split()[1:]],
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

    def write_bnx_for_list_of_molecules(self, list_of_molecule_ids, output_file_name, bnx_head_file="./bnx_head.txt"):
        """
        Filters bnx file by given list of molecule ids.
        :param list_of_molecule_ids: List of bnx molecule ids
        :param output_file_name:
        :param bnx_head_file:
        :return:
        """
        bnx_header_file = open(bnx_head_file, "r")
        bnx_head_text = bnx_header_file.read()
        bnx_header_file.close()
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
        self.read_and_load_cmap_file()
        self.get_intervals()
        self.obtain_molecule_match_info()

    def read_and_load_cmap_file(self):
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
