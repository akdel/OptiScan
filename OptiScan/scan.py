from OptiScan import *
from OptiScan.rotation_optimisation import *
from scipy import ndimage
from OptiScan.folder_searcher import FolderSearcher, SaphyrFileTree
from os import listdir
import imagecodecs as imageio
import gc
from PIL import Image

"""
TODO 1. split each col frame into frames and then stich them.
TODO 2. add the third set of frames as the second channel.
TODO 3. focus both channels on the backbone channel.
"""

class SaphyrColumn:
    def __init__(self):
        pass

class Scan:
    """
    This is for reading and recording the scan frames. Add lower level functionality only.
    """
    def __init__(self, tiff_file_location: str, chip_dimension=(12, 95), scan_no=1, load_frames=True,
                 saphyr=False, saphyr_folder_searcher=False):
        """
        Parameters
        ----------
        tiff_file_location : Only compulsory input; the BNG Image frames for a scan.
        chip_dimension : This is to be set if the BNG chip dimensions are not ordinary.
                         Default/ordinary dimension is (12, 95)
        scan_no : Scan id.
        load_frames : If True, loads frames when the class is initiated.
        """
        self.saphyr = saphyr
        self.scan_id = scan_no
        self.tiff_file_location = tiff_file_location
        self.channel_num = 2
        if not saphyr:
            self.im = Image.open(self.tiff_file_location)
            self.nof_frames = self.im.n_frames
            self.frames = list()
            self.frame_dimension = None
            self.chip_dimension = chip_dimension
            print("number of frames: %s" % self.nof_frames)
            if load_frames:
                self._record_frames()
        else:
            self.frames = saphyr_folder_searcher
            self.frames["CH1"] = sorted(self.frames["CH1"])
            self.frames["CH2"] = sorted(self.frames["CH2"])
            if len(self.frames["CH3"]):
                self.frames["CH3"] = sorted(self.frames["CH3"])
                self.channel_num = 3
            self.nof_frames = len(self.frames["CH1"])
            self.chip_dimension = (1, self.nof_frames)

    def _record_frames(self):
        """
        Loads all the frames to memory (stored as numpy arrays in self.frames).
        """
        dim_get = True
        for frame_no in range(self.nof_frames):
            self.im.seek(frame_no)
            scan_array = np.array(self.im)
            self.frames.append(scan_array)
            if dim_get:
                self.frame_dimension = scan_array.shape
                dim_get = False


class AnalyzeScan(Scan):
    """
    This class is for molecule boundary detection and signal extraction from a single BNG scan data-set.
    """
    def __init__(self, tiff_file_location, chip_dimension=(12, 95), scan_no=1,
                 load_frames=True, saphyr=False, saphyr_folder_searcher=False):
        """
        Initiation function. Tiff image frames from BNG is the only compulsory data to be provided.
        Parameters
        ----------
        tiff_file_location : path to the tiff file.
        chip_dimension : This is to be set if the BNG chip dimensions are not ordinary.
                         Default/ordinary dimension is (12, 95)
        scan_no : Scan id.
        """
        Scan.__init__(self, tiff_file_location, chip_dimension=chip_dimension,
                      scan_no=scan_no, load_frames=load_frames, saphyr=saphyr, 
                      saphyr_folder_searcher=saphyr_folder_searcher)
        self.current_column_id = 0
        self.current_lab_column = None
        self.current_lab_column2 = None
        self.current_mol_column = None
        self.mol_slices = list()
        self.lab_slices = list()
        self.lab_slices2 = list()
        self.slice_coordinates = list()

        self.molecules = {i: list() for i in range(self.chip_dimension[1])}
        self.column_info = {i: {"extracted": False,
                                "number_of_molecules": 0,
                                "mean_molecule_length": 0,
                                "abstract_threshold": 0,
                                "minimum_allowed_length": 0,
                                "longest_shape": 0,
                                "mem_start": 0,
                                "mem_end": 0
                                } for i in range(self.chip_dimension[1])}
        self.memmap_path = "%s_scan.mem" % scan_no
        self.memmap_status = False

    def create_molecule(self, molecule_id: str, nick_label: np.ndarray, backbone_label: np.ndarray,
                        check_function=lambda x,y:(x,y)):
        """
        Molecule check point. Adds checked molecule to self.molecules.
        Parameters
        ----------
        molecule_id : Id of the molecule.
        nick_label : Nick label 1d signal
        backbone_label : Backbone label 1d signal
        check_function : Filtering function which can be used optionally to filter out unfit signals.
        """
        if nick_label.shape != backbone_label.shape:
            print("backbone and nicklength did not match!")
            return None
        outcome = check_function(nick_label, backbone_label)
        if not outcome:
            return None
        else:
            nick_label, backbone_label = outcome
            self.molecules[self.current_column_id].append((molecule_id, nick_label, backbone_label))

    def create_molecule_additional(self, molecule_id: str, nick_label: np.ndarray, nick_label2: np.ndarray,
                                   backbone_label: np.ndarray,
                                   check_function=lambda x,y,z:(x,y,z)):
        if nick_label.shape != backbone_label.shape or nick_label.shape != nick_label2.shape:
            print("backbone and nicklength did not match!")
            return None
        outcome = check_function(nick_label, nick_label2, backbone_label)
        if not outcome:
            return None
        else:
            nick_label, nick_label2, backbone_label = outcome
            self.molecules[self.current_column_id].append((molecule_id, nick_label, backbone_label, nick_label2))

    def define_molecules(self, minimum_molecule_length=50*5, abstraction_threshold=100, m=1):
        """
        Defines molecule boundaries in the current column and updates the corresponding metadata.
        Parameters
        ----------
        minimum_molecule_length : Minimum length of a detected molecule in pixels.
        abstraction_threshold : The threshold for abstraction. This shold be defined according to the dataset.
        Returns
        -------
        Molecule slices i.e.// mol_slices, lab_slices, slice_coordinates
        """
        self.column_info[self.current_column_id]["minimum_allowed_length"] = int(minimum_molecule_length)
        self.column_info[self.current_column_id]["abstract_threshold"] = int(abstraction_threshold)
        print(abstraction_threshold, "absthr")
        if not self.saphyr:
            print("irys column")
            molecule_abstract = np.where(self.current_mol_column > abstraction_threshold, self.current_mol_column, 0)
            mask = np.array(list(signal.ricker(16, 1)) * 5).reshape((5, -1))
            molecule_abstract = ndimage.convolve(molecule_abstract.astype(float), mask.astype(float))
            molecule_abstract = np.where(molecule_abstract > (abstraction_threshold * m * mask.shape[0]), 1, 0)
            # molecule_abstract = np.zeros(self.current_mol_column.shape)
            # mask = np.array([[1., -1.], [1., -1.], [1., -1.], [1., -1.]])
            # molecule_abstract = np.maximum(ndimage.convolve(self.current_mol_column, mask), molecule_abstract)
            # molecule_abstract = np.where(molecule_abstract > abstraction_threshold, 1, 0)
        else:
            print("saphyr column")
            molecule_abstract = np.where(self.current_mol_column > abstraction_threshold, self.current_mol_column, 0)
            mask = np.array(list(signal.ricker(16, 1)) * 5).reshape((5, -1))
            molecule_abstract = ndimage.convolve(molecule_abstract.astype(float), mask.astype(float))
            molecule_abstract = np.where(molecule_abstract > (abstraction_threshold*m*mask.shape[0]), 1, 0)
            # self.something = molecule_abstract
        labels = ndimage.label(molecule_abstract)
        slices = ndimage.find_objects(labels[0], max_label=labels[1])
        print("%s potential molecule fragments are found..." % len(slices))
        slices = [x for x in slices if (x[0].stop - x[0].start) >= minimum_molecule_length]
        self.slice_coordinates = slices
        self.mol_slices = list()
        self.mol_slices = list() # TO DO look at this!
        print("%s of these are longer than %s pixels (about %s kbp) and taken as molecule set." % \
              (len(slices), minimum_molecule_length, (minimum_molecule_length * 500)/1000))
        for _slice in slices:
            if self.saphyr:
                new_slice = np.s_[_slice[0].start:_slice[0].stop, _slice[1].start:_slice[1].stop] #+3]
            else:
                new_slice = np.s_[_slice[0].start:_slice[0].stop, _slice[1].start:_slice[1].stop] #+3]
            self.mol_slices.append(self.current_mol_column[new_slice])
            self.lab_slices.append(self.current_lab_column[new_slice])
            if self.channel_num == 3:
                self.lab_slices2.append(self.current_lab_column2[new_slice])
            else:
                continue

    def obtain_column_summary(self, column_id: int):
        """
        Obtains summary values for the current column.
        Parameters
        ----------
        column_id : Column id in the scan
        Returns
        -------
        fills in number of molecules and mean molecule length fields in self.column_info.
        """
        self.column_info[column_id]["number_of_molecules"] = len(self.molecules[column_id])
        if self.channel_num == 3:
            mean_length = np.mean([len(mol) for mol_id, nick, mol, nick2 in self.molecules[column_id]])
        else:
            mean_length = np.mean([len(mol) for mol_id, nick, mol in self.molecules[column_id]])
        try:
            self.column_info[column_id]["mean_molecule_length"] = int(mean_length)
        except ValueError:
            self.column_info[column_id]["mean_molecule_length"] = 0
        self.column_info[column_id]["extracted"] = True

    def _get_1d_signal(self, molecule_id: int) -> (np.ndarray, np.ndarray):
        """
        Obtains a single given signal for a given molecule id inside the column.
        Parameters
        ----------
        molecule_id : Column molecule id (different from scan molecule id)
        Returns
        -------
        (nick signal, backbone signal)
        """
        label_slice = self.lab_slices[molecule_id]
        backbone_slice = self.mol_slices[molecule_id]
        signal_set = label_slice[:,:]
        try:
            signal_averages = np.average(signal_set, axis=0)
        except ZeroDivisionError:
            print("zero division error")
            return np.array([]), np.array([])
        index = np.argmax(signal_averages)
        top_signal = signal_set[:, index]
        backbone_signal = backbone_slice[:, index]
        return top_signal, backbone_signal

    def _get_1d_signal_additional(self, molecule_id: int):
        label_slice = self.lab_slices2[molecule_id]
        backbone_slice = self.mol_slices[molecule_id]
        signal_set = label_slice[:,:]
        try:
            signal_averages = np.average(signal_set, axis=0)
        except ZeroDivisionError:
            print("zero division error")
            return np.array([]), np.array([])
        index = np.argmax(signal_averages)
        top_signal = signal_set[:, index]
        return top_signal

    def get_all_signals(self):
        """
        Obtains all signals from the current column.
        Returns
        -------
        Signals for the current column.
        """
        for i in range(len(self.lab_slices)):
            molecule_id = str(self.current_column_id) + str(i)
            nick_label, backbone_label = self._get_1d_signal(i)
            if self.channel_num == 3:
                nick_label2 = self._get_1d_signal_additional(i)
                self.create_molecule_additional(molecule_id, nick_label, nick_label2, backbone_label)
            else:
                self.create_molecule(molecule_id, nick_label, backbone_label)
        self.obtain_column_summary(self.current_column_id)

    def _clean_and_progress_column(self):
        """
        Progresses column to next.
        """
        if self.current_column_id == self.chip_dimension[1]:
            pass
        else:
            self.current_column_id += 1
        self.current_lab_column = None
        self.current_mol_column = None
        self.current_lab_column2 = None
        self._clean_slices()
        gc.collect()

    def _clean_slices(self):
        """
        Cleans molecule slices.
        """
        self.mol_slices = []
        self.lab_slices = []
        self.lab_slices2 = []
        self.slice_coordinates = []

    def stitch_and_load_column(self):
        """
        Stitches the current columns and loads them for molecule detection.
        Returns
        -------
        Stitched nick and backbone columns for the current column id.
        """
        if not self.saphyr:
            backbone_label_column, nick_label_column = return_column(self.frames, self.current_column_id,
                                                                    self.chip_dimension)
            backbone_label_column, [nick_label_column] = stitch_column(backbone_label_column, [nick_label_column])
            self.current_lab_column = ndimage.zoom(nick_label_column, 0.5)
            self.current_mol_column = ndimage.zoom(backbone_label_column, 0.5)
        else:
            print(self.frames["CH1"][self.current_column_id])
            concat_mol_column = imageio.imread(self.frames["CH1"][self.current_column_id]).astype(float)
            concat_lab_column = imageio.imread(self.frames["CH2"][self.current_column_id]).astype(float)
            if self.channel_num == 3:
                concat_lab_column_2nd_channel = imageio.imread(self.frames["CH3"][self.current_column_id]).astype(float)
            backbone_frames = [ndimage.white_tophat(concat_mol_column[i:i+2048], structure=np.ones((1,3))) for i in range(0, concat_mol_column.shape[0], 2048)]
            if self.channel_num == 3:
                lab_frames = [[ndimage.white_tophat(concat_lab_column[i:i + 2048], structure=disk(6)) for i in
                              range(0, concat_lab_column.shape[0], 2048)],
                              [ndimage.white_tophat(concat_lab_column_2nd_channel[i:i + 2048], structure=disk(6)) for i in
                               range(0, concat_lab_column_2nd_channel.shape[0], 2048)]]
            else:
                lab_frames = [[ndimage.white_tophat(concat_lab_column[i:i + 2048], structure=disk(6)) for i in
                               range(0, concat_lab_column.shape[0], 2048)]]
            # print([f.shape for f in backbone_frames],[ff.shape for ff in lab_frames])
            if self.channel_num == 3:
                self.current_mol_column, (self.current_lab_column, self.current_lab_column2) = stitch_column(backbone_frames, lab_frames, saphyr=self.saphyr)
            else:
                self.current_mol_column, [self.current_lab_column] = stitch_column(backbone_frames, lab_frames, saphyr=self.saphyr)
            
            # self.current_mol_column = ndimage.white_tophat(self.current_mol_column, structure=disk(12))
            # self.current_lab_column = ndimage.white_tophat(self.current_lab_column , structure=disk(12))

    def annotate_column(self, intensity=1000) -> np.ndarray:
        """
        Annotates molecule backbone image by the defined molecule boundaries.
        Parameters
        ----------
        intensity : Intensity of annotation lines.
        Returns
        -------
        Annotated column image.
        """
        assert len(self.mol_slices) != 0
        assert len(self.lab_slices) != 0
        annotated_mol = np.array(self.current_mol_column)
        for _slice in self.slice_coordinates:
            y, x = _slice
            annotated_mol[y.start:y.stop, x.start:x.stop] = intensity
        return annotated_mol

    def stitch_extract_molecules_in_column(self, minimum_molecule_length=50*5, abstraction_threshold=100):
        """
        Stitches the current column and extracts the molecule signals.
        Parameters
        ----------
        minimum_molecule_length : minimum molecule length to be extracted
        abstraction_threshold : this parameter is used to define molecules from the background. Higher is more strict
                                and may split molecules. Lower can lead to false positives and merged molecules.
        Returns
        -------
        Molecule signals in the current column. Can be accessed from self.molecules.
        """
        self.stitch_and_load_column()
        if self.saphyr:
            minimum_molecule_length *= 0.6
        else:
            pass
        self.define_molecules(minimum_molecule_length=int(minimum_molecule_length),
                              abstraction_threshold=abstraction_threshold)
        # frame func call
        self.get_all_signals()

    def stitch_extract_molecules_in_scan(self, minimum_molecule_length=50*5, abstraction_threshold=100):
        """
        This is the pipeline function to extract all molecules from the scan.
        Parameters
        ----------
        minimum_molecule_length : minimum molecule length to be extracted
        abstraction_threshold : this parameter is used to define molecules from the background. Higher is more strict
                                and may split molecules. Lower can lead to false positives and merged molecules.
        Returns
        -------
        Molecule signals per column. Can be accessed from self.molecules.
        """
        for _ in range(self.chip_dimension[1]):
            self.stitch_extract_molecules_in_column(minimum_molecule_length=minimum_molecule_length,
                                                    abstraction_threshold=abstraction_threshold)
            self._clean_and_progress_column()

    def initiate_memmap(self):
        """
        Initiates the memmap file.
        Returns
        -------
        Writes 8 bytes into memmap file and closes it. Turns on self.memmap_status to True.
        """
        assert not self.memmap_status
        self.memmap_status = True
        memmap = np.memmap(self.memmap_path, mode="w+", shape=1, dtype=np.float64)
        memmap[:] = 0.

    def add_column_into_memmap(self):
        """
        Adds molecules found in the self.current_column into memmap file. This should be done sequentially.
        i.e.// start writing from column 0 and proceed until last column. This is stepwise to solve interraptions..
        Returns
        -------
        Column molecules in memmap file.
        """
        assert self.column_info[self.current_column_id]["extracted"]
        dsize = np.array([0], dtype=float).itemsize
        if not self.memmap_status:
            self.initiate_memmap()
        if self.current_column_id == 0:
            start_mem = 8
        else:
            start_mem = self.column_info[self.current_column_id - 1]["mem_end"]

        if self.column_info[self.current_column_id]["number_of_molecules"] == 0:
            self.column_info[self.current_column_id]["mem_start"] = start_mem
            self.column_info[self.current_column_id]["mem_end"] = start_mem
            return None

        longest_shape = max([len(x[1]) for x in self.molecules[self.current_column_id]]) + 1
        self.column_info[self.current_column_id]["longest_shape"] = longest_shape
        number_of_molecules = self.column_info[self.current_column_id]["number_of_molecules"]
        end_mem = (dsize * 2 * number_of_molecules * longest_shape) + start_mem
        self.column_info[self.current_column_id]["mem_start"] = start_mem
        self.column_info[self.current_column_id]["mem_end"] = end_mem
        shift = 0
        for mol_id, nick_label, back_label in self.molecules[self.current_column_id]:
            memmap = np.memmap(self.memmap_path, mode="r+",
                               offset=start_mem + dsize * shift * longest_shape * 2,
                               shape=longest_shape * 2, dtype=np.float64)
            memmap[:nick_label.shape[0]] = nick_label
            memmap[longest_shape:longest_shape + back_label.shape[0]] = back_label
            memmap[-1] = np.float64(back_label.shape[0])
            shift += 1
            del memmap

    def _read_column(self, column_id: int) -> np.ndarray:
        """
        Reads a column with a given column id.
        Parameters
        ----------
        column_id : column id.

        Returns
        -------
        An array containing column data as float64.
        """
        assert self.column_info[column_id]["mem_start"] != 0
        start_mem = self.column_info[column_id]["mem_start"]
        end_mem = self.column_info[column_id]["mem_end"]
        assert (end_mem - start_mem) % np.float64().itemsize == 0
        column_shape = int((end_mem - start_mem) / np.float64().itemsize)
        new_array = np.zeros(column_shape, dtype=np.float64)
        memmap = np.memmap(self.memmap_path, dtype=np.float64, mode="r", offset=start_mem, shape=column_shape)
        new_array[:] = memmap[:]
        del memmap
        return new_array

    def return_signals_in_column(self, column_id: int):
        """
        Extracts signals from column data for a given column id.
        Parameters
        ----------
        column_id : column id

        Returns
        -------
        Stream of signal pairs.
        """
        column_data = self._read_column(column_id)
        longest_mol_shape = self.column_info[column_id]["longest_shape"]
        number_of_molecules = self.column_info[column_id]["number_of_molecules"]
        column_data = column_data.reshape((number_of_molecules * 2, longest_mol_shape))
        for i in range(0, column_data.shape[0], 2):
            molecule_length = int(column_data[i + 1, -1])
            yield column_data[i, :molecule_length], column_data[i + 1, :molecule_length]

    def restart_scan(self):
        """
        Restarts the scan by removing all data in class.
        Returns
        -------
        Fresh class.
        """
        self.__init__(self.tiff_file_location, self.chip_dimension, self.scan_id, load_frames=False)


class SaphyrExtract:
    def __init__(self, run_folder: str):
        self.saphyr_tree = SaphyrFileTree.from_run_folder(run_folder)
        self.bank_info = dict()
        for bank_id in self.saphyr_tree.get_bank_level_ids():
            if bank_id in self.bank_info:
                self.bank_info[bank_id] += self.saphyr_tree.get_image_paths_for_bank_path(bank_id)
            else:
                self.bank_info[bank_id] = self.saphyr_tree.get_image_paths_for_bank_path(bank_id)

    def split_channels_for_bankid(self, bank_id: str):
        images = self.bank_info[bank_id]
        channels = {"CH1": list(),
                    "CH2": list(),
                    "CH3": list()}
        for image in images:
            for channel in channels.keys():
                if channel in image:
                    channels[channel].append(self.saphyr_tree.root + image)
                else:
                    continue
        return channels

    def get_molecules_for_bank_id(self, bank_id: str, abstraction_threshold=100):
        print(self.bank_info[bank_id])
        bank = AnalyzeScan(self.saphyr_tree.root, chip_dimension=(12, 95),
                           scan_no=self.saphyr_tree.root + bank_id,
                           load_frames=True, saphyr=True,
                           saphyr_folder_searcher=self.split_channels_for_bankid(bank_id))
        bank.stitch_extract_molecules_in_scan(abstraction_threshold=abstraction_threshold)
        for col_id in sorted(list(bank.column_info.keys())):
            bank.current_column_id = col_id
            if bank.channel_num == 3:
                np.save(f"{self.saphyr_tree.root + bank_id}/{col_id}_label2.npy",
                        [x[3] for x in bank.molecules[col_id]])
            np.save(f"{self.saphyr_tree.root + bank_id}/{col_id}_label.npy", [x[1] for x in bank.molecules[col_id]])
            np.save(f"{self.saphyr_tree.root + bank_id}/{col_id}_backbone.npy", [x[2] for x in bank.molecules[col_id]])

        del bank

    def save_molecules_for_all_bank_ids(self, abstraction_threshold=100):
        for bank_id in self.bank_info:
            self.get_molecules_for_bank_id(bank_id, abstraction_threshold=abstraction_threshold)

    def save_molecules_for_bank_ids(self, bank_ids: list, abstraction_threshold=100):
        for bank_id in bank_ids:
            print(bank_id)
            print()
            self.get_molecules_for_bank_id(bank_id, abstraction_threshold=abstraction_threshold)


class Run(FolderSearcher):
    """
    This class extends FolderSearcher to locate BNG files and split them into scan. Each scan is then loaded
    with AnalyzeScan.
    """
    def __init__(self, run_directory, chip_dimension=(12, 95), saphyr=False):
        """
        Initiates class.
        Parameters
        ----------
        run_directory : BNG run directory which holds at least a set of scan tiffs.
        chip_dimension : This is user defined according to chip version.
        """
        self.saphyr = saphyr
        if saphyr:
            FolderSearcher.__init__(self, run_directory, saphyr=True)
            print(self.scans.keys())
            self.analyzed_scans = {i:AnalyzeScan(run_directory, chip_dimension=chip_dimension,
                                                scan_no=run_directory + str(i),
                                                load_frames=False, saphyr=saphyr,
                                                saphyr_folder_searcher=self.scans[i]["tiff_location"]) for i in self.scans.keys()}
        else:
            FolderSearcher.__init__(self, run_directory)
            self.analyzed_scans = {i:AnalyzeScan(self.scans[i]["tiff_location"], chip_dimension=chip_dimension,
                                                scan_no=run_directory + str(i),
                                                load_frames=False) for i in self.scans.keys()}

    def read_frames_in_scan(self, scan_id: int):
        """
        Loads tiff images into scan.
        Parameters
        ----------
        scan_id : scan id

        Returns
        -------
        tiff images in given scan class. ie// scan.frames
        """
        if self.saphyr:
            pass
        else:
            self.analyzed_scans[scan_id]._record_frames()

    def read_frames_in_all_scans(self):
        """
        Reads frames in all scans.
        Returns
        -------
        tiff images in all scan classes.
        """
        for scan_id in self.analyzed_scans.keys():
            self.read_frames_in_scan(scan_id)

    def extract_molecules_for_db(self, scan_id: int, abstraction_threshold=100):
        """
        Molecules are extracted in a way that OptiScan database will accept.
        Parameters
        ----------
        scan_id : scan id in run

        Returns
        -------

        """
        scan = self.analyzed_scans[scan_id]
        scan.stitch_extract_molecules_in_scan(abstraction_threshold=abstraction_threshold)
        scan.initiate_memmap()
        for col_id in sorted(list(scan.column_info.keys())):
            scan.current_column_id = col_id
            scan.add_column_into_memmap()


class Runs:
    """
    This collects all runs under one hood.
    """
    def __init__(self, run_directories: [str], chip_dimension=(12, 95), saphyr=False):
        """
        Initiates class.
        Parameters
        ----------
        run_directories : BNG run directories in which each one of them at least have a set of scan tiffs.
        chip_dimension : This is user defined according to chip version.
        """
        self.saphyr = saphyr
        self.analyzed_runs = {x: Run(x, chip_dimension=chip_dimension, saphyr=saphyr) for x in run_directories}

    def extract_molecules_from_scan_in_run(self, run_id: str, scan_id: int, abstraction_threshold=100):
        """
        See Run.extract_molecules_for_db
        Parameters
        ----------
        run_id : run id
        scan_id : scan id in given run id

        Returns
        -------

        """
        run = self.analyzed_runs[run_id]
        run.read_frames_in_scan(scan_id)
        run.extract_molecules_for_db(scan_id, abstraction_threshold=abstraction_threshold)

    def extract_molecules_in_run(self, run_id: str):
        """
        Extracts all molecules in a run by going over all scans in that run.
        Parameters
        ----------
        run_id : run id

        Returns
        -------

        """
        run = self.analyzed_runs[run_id]
        for scan_id in run.scans.keys():
            self.extract_molecules_from_scan_in_run(run_id, scan_id)

    def _flatten_scans_for_parallel(self):
        for run_id in self.analyzed_runs.keys():
            for scan_id in self.analyzed_runs[run_id].analyzed_scans.keys():
                yield (run_id, scan_id)

    def split_scans_to_threads(self, number_of_threads: int):
        scan_jobs = list(self._flatten_scans_for_parallel())
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


def return_column(image_frames: [np.ndarray], column_no: int, dimensions: (int, int)) -> ([np.ndarray], [np.ndarray]):
    """
    Returns a list of nick and a list of backbone frames in a column from a provided image frame set.
    Parameters
    ----------
    image_frames : Full set of image frames
    column_no : Column number
    dimensions : Chip dimension in frames
    Returns
    -------
    Sequential set of backbone and nick frames.
    """
    assert dimensions[1] >= column_no
    assert len(image_frames) / 2 == dimensions[0] * dimensions[1]
    bstart = column_no * (dimensions[0] * 2)
    bend = bstart + (dimensions[0] * 2)
    nstart = bstart + 1
    nend = bend + 1
    if column_no % 2 != 0:
        nick_frames = [ndimage.zoom(image_frames[i].astype(float), 2).astype(int) for i in range(nstart, nend, 2)]
        backbone_frames = [ndimage.zoom(image_frames[i].astype(float), 2).astype(int) for i in range(bstart, bend, 2)]
    else:
        nick_frames = [ndimage.zoom(image_frames[i][::-1].astype(float), 2).astype(int) for i in range(nstart, nend, 2)]
        backbone_frames = [ndimage.zoom(image_frames[i][::-1].astype(float), 2).astype(int) for i in range(bstart, bend, 2)]
    return backbone_frames, nick_frames


def stitch_column(backbone_frames: [np.ndarray], nick_frames: list, saphyr=False) -> (np.ndarray, np.ndarray):
    """
    Wrapper function which
    Parameters
    ----------
    backbone_frames : A sequential set of backbone frames.
    nick_frames :A sequential set of nick frames.
    Returns
    -------
    Stitched column.
    """
    return merging_with_rotation_optimisation_and_xshift(backbone_frames,
                                                         additional_sets=nick_frames,
                                                         y_shift=True,
                                                         tophat=True,
                                                         magnification_optimisation=True, saphyr=saphyr)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    # sc = AnalyzeScan("/Users/akdel/PycharmProjects/projects3/OptiScan/data/", saphyr=True,
    #           saphyr_folder_searcher=["/Users/akdel/PycharmProjects/projects3/OptiScan/data/Bank/B1_CH1_C100.tif",
    #                                   "/Users/akdel/PycharmProjects/projects3/OptiScan/data/Bank/B1_CH2_C100.tif"])
    # sc = AnalyzeScan("/Users/akdel/PycharmProjects/projects3/OptiScan/data/ses/Exp2066_BspQ1_BbvCI_2014-07-01_20_00_Scan01.tiff",
    #                  saphyr=False, chip_dimension=(12, 108))
    sc = AnalyzeScan(
        "/Volumes/Seagate Expansion Drive/2014-09/Exp2066_DNA_Bs+Bb_2014-09-01_10_32/Exp2066_DNA_Bs+Bb_2014-09-01_10_32_Scan10.tiff",
        saphyr=False, chip_dimension=(12, 108))
    sc.current_column_id = 50
    sc.stitch_extract_molecules_in_column(minimum_molecule_length=250, abstraction_threshold=110)
    plt.figure(figsize=(10, 100))
    plt.imshow(sc.annotate_column(1000))
    plt.show()
    # sc.add_column_into_memmap()
    # molecules = [x for x in sc.return_signals_in_column(0)]
    # print(len(molecules))
    # for back, nick in molecules:
    #     plt.plot(back)
    #     plt.plot(nick)
    #     plt.show()