from Photo_map.database_orms import *
from Photo_map.scan import Runs
from Photo_map import *
import types
from Photo_map.utils import get_bnx_info_with_upsampling


class MoleculeDB(Runs):
    """
    Molecule database creator. Molecule extraction can be performed in batch as it inherits Runs class.
    """
    def __init__(self, database_location: str, bionano_directories: [str], organism: str,
                 chip_dimension=(12, 95), date="2016"):
        """
        Initiates class.
        Parameters
        ----------
        database_location :   Location of database.
        bionano_directories : List of bionano directories.
        organism :            Organism name
        chip_dimension :      Chip dimension
        date :                Date
        """
        bionano_directories = [fix_directory_end(x) for x in bionano_directories]
        Runs.__init__(self, bionano_directories, chip_dimension=chip_dimension)
        self.db_location = database_location
        self.db = None
        self.organism = organism
        self.date = date
        self.ids = {x: dict() for x in bionano_directories}

    def create_tables_first_time(self):
        """
        Returns
        -------

        """
        engine = sq.create_engine("sqlite:///%s" % self.db_location, connect_args={"timeout": 1000})
        Base.metadata.create_all(engine)
        self.db = sessionmaker(bind=engine)
        self.db = self.db()
        self.set_tables_in_db()

    def connect_db(self):
        engine = sq.create_engine("sqlite:///%s" % self.db_location, connect_args={"timeout": 1000})
        Base.metadata.bind = engine
        self.db = sessionmaker(bind=engine)
        self.db = self.db()

    def set_tables_in_db(self):
        unique_id = int()
        for run in self.analyzed_runs:
            print(run)
            _run = RunLineage(id=run, organism=self.organism, date=self.date)
            self.db.add(_run)
            self.db.commit()
            for scan_id in self.analyzed_runs[run].scans.keys():
                print(scan_id)
                print(self.analyzed_runs[run].scans[scan_id])
                scan_ana = self.analyzed_runs[run].analyzed_scans[scan_id]
                current_scan = self.analyzed_runs[run].scans[scan_id]
                self.ids[run][scan_id] = unique_id
                scan = ScanLinage(id=unique_id, scan=scan_id, run_id=run,
                                  bnx_path=current_scan["RawMolecules_file_location"],
                                  mol_path=current_scan["Molecules_file_location"],
                                  stitch_path=current_scan["Stitch_file_location"],
                                  lab_path=current_scan["Labels_file_location"],
                                  tiff_path=current_scan["tiff_location"],
                                  mem_path=scan_ana.memmap_path)
                unique_id += 1
                self.db.add(scan)
                self.db.commit()

    def _load_unique_scan_ids_from_db(self):
        assert self.db
        self.ids = {x.run_id: {} for x in self.db.query(ScanLinage).all()}
        for entry in self.db.query(ScanLinage).all():
            self.ids[entry.run_id][entry.scan] = entry.id

    def insert_column(self, run_id: str, scan_id: int, column_id: int):
        column_info = self.analyzed_runs[run_id].analyzed_scans[scan_id].column_info[column_id]
        assert column_info["extracted"]
        assert column_info["mem_start"] != 0
        try:
            unique_scan_id = self.ids[run_id][scan_id]
        except KeyError:
            self._load_unique_scan_ids_from_db()
            unique_scan_id = self.ids[run_id][scan_id]
        column = ColArray(id="%s_%s_%s" % (run_id, scan_id, column_id),
                          scan_id=unique_scan_id, col_id=column_id,
                          mean_molecule_length=column_info["mean_molecule_length"],
                          molecule_count=column_info["number_of_molecules"],
                          minimum_molecule_length=column_info["minimum_allowed_length"],
                          abstract_threshold=column_info["abstract_threshold"],
                          memmap_start=column_info["mem_start"],
                          memmap_end=column_info["mem_end"],
                          longest_shape=column_info["longest_shape"]
                          )
        self.db.add(column)
        self.db.commit()
        col_shape = (column_info["mem_end"] - column_info["mem_start"]) / np.float64().itemsize
        assert col_shape % 2 == 0
        try:
            number_of_molecules = col_shape / (column_info["longest_shape"] * 2)
        except ZeroDivisionError:
            number_of_molecules = 0
        assert number_of_molecules == int(number_of_molecules)

    def insert_columns_for_scan(self, run_id: str, scan_id: int):
        run = self.analyzed_runs[run_id]
        scan = run.analyzed_scans[scan_id]
        for col_id in scan.column_info.keys():
            self.insert_column(run_id, scan_id, col_id)


class MoleculeConnector:
    """
    This is for loading molecules from a OptiScan database.
    """
    def __init__(self, database_location: str):
        """
        Initiates by connecting to database and loading all metadata.
        Parameters
        ----------
        database_location : path to database.
        """
        self.db_location = database_location
        self.db = None
        self.connect_db()
        self.db_runs = dict()
        self.date = None
        self.organism = None
        self.load_db_to_class()

    def connect_db(self):
        """
        Connects to database.
        Returns
        -------
        self.db as an sqlalchemy session.
        """
        engine = sq.create_engine("sqlite:///%s" % self.db_location)
        Base.metadata.bind = engine
        self.db = sessionmaker(bind=engine)
        self.db = self.db()

    def load_run_info(self):
        """
        Loads run metadata.
        Returns
        -------
        Stored in self.db_runs
        """
        run_info = [(x.id, x.organism, x.date) for x in self.db.query(RunLineage).all()]
        for run in run_info:
            self.db_runs[run[0]] = dict()
            self.date = run[2]
            self.organism = run[1]

    def load_scan_and_column_info(self):
        """
        Loads scan and column metadata.
        Returns
        -------
        Stored in self.db_runs
        """
        for scan_element in self.db.query(ScanLinage).all():
            _id = scan_element.id
            scan_id = scan_element.scan
            run_id = scan_element.run_id
            self.db_runs[run_id][scan_id] = dict()
            self.db_runs[run_id][scan_id]["id"] = _id
            self.db_runs[run_id][scan_id]["mem_path"] = scan_element.mem_path
            self.db_runs[run_id][scan_id]["columns"] = {x.col_id: {
                "mean_molecule_length": x.mean_molecule_length,
                "molecule_count": x.molecule_count,
                "minimum_molecule_length": x.minimum_molecule_length,
                "abstract_threshold": x.abstract_threshold,
                "memmap_start": x.memmap_start,
                "memmap_end": x.memmap_end,
                "longest_shape": x.longest_shape,
                "unique_id": x.id} for x in self.db.query(ColArray).all() if x.scan_id == _id}

    def load_db_to_class(self):
        """
        Loads all metadata from db.
        Returns
        -------
        Stored in self.db_runs
        """
        self.load_run_info()
        self.load_scan_and_column_info()

    def yield_molecule_signals_in_column(self, run_id: str, scan_id: int, column_id: int):
        """
        Streams out molecules.
        Parameters
        ----------
        run_id : run id
        scan_id : scan id in given run id
        column_id : column id in given scan id

        Returns
        -------
        Stream of molecule pairs
        """
        column_data = self._read_column(run_id, scan_id, column_id)
        column_info = self.db_runs[run_id][scan_id]["columns"]
        longest_mol_shape = column_info[column_id]["longest_shape"]
        number_of_molecules = column_info[column_id]["molecule_count"]
        column_data = column_data.reshape((number_of_molecules * 2, longest_mol_shape)) # TODO use this in the molecule retrival part
        for i in range(0, column_data.shape[0], 2):
            molecule_length = int(column_data[i + 1, -1])
            yield column_data[i, :molecule_length], column_data[i + 1, :molecule_length]

    def yield_molecule_signals_in_scan(self, run_id: str, scan_id: int):
        """
        Streams out molecules.
        Parameters
        ----------
        run_id : run id
        scan_id : scan id in given run id

        Returns
        -------
        Stream of molecule pairs
        """
        column_info = self.db_runs[run_id][scan_id]["columns"]
        for i in column_info:
            for signal_pair in self.yield_molecule_signals_in_column(run_id, scan_id, i):
                yield signal_pair

    def yield_molecule_signals_in_run(self, run_id: str):
        """
        Streams out molecules.
        Parameters
        ----------
        run_id : run id

        Returns
        -------
        Stream of molecule pairs
        """
        run_info = self.db_runs[run_id]
        for scan_id in run_info:
            for signal_pair in self.yield_molecule_signals_in_scan(run_id, scan_id):
                yield signal_pair

    def yield_molecule_signals_in_all_runs(self):
        """
        Streams out molecules for all runs.

        Returns
        -------
        Stream of molecule pairs
        """
        for run_id in self.db_runs:
            for signal_pair in self.yield_molecule_signals_in_run(run_id):
                yield signal_pair

    def yield_molecule_info_in_database(self):
        n = 0
        for run_id in self.db_runs:
            for scan_id in self.db_runs[run_id]:
                scan_info = self.db_runs[run_id][scan_id]
                for col_id in scan_info["columns"]:
                    index_in_col = -1
                    column_info = scan_info["columns"][col_id]
                    unique_col_id = column_info["unique_id"]
                    for molecule in range(column_info["molecule_count"]):
                        n += 1
                        index_in_col += 1
                        yield (n, unique_col_id, index_in_col)

    def _read_column(self, run_id: str, scan_id: int, column_id: int) -> np.ndarray:
        """
        Reads column data from memmap.
        Parameters
        ----------
        run_id : run id
        scan_id : scan id in given run id
        column_id : column id in given scan id

        Returns
        -------

        """
        column_info = self.db_runs[run_id][scan_id]["columns"]
        assert column_info[column_id]["memmap_start"] != 0
        start_mem = column_info[column_id]["memmap_start"]
        end_mem = column_info[column_id]["memmap_end"]
        assert (end_mem - start_mem) % np.float64().itemsize == 0
        column_shape = int((end_mem - start_mem) / np.float64().itemsize)
        new_array = np.zeros(column_shape, dtype=np.float64)
        memmap_path = self.db_runs[run_id][scan_id]["mem_path"]
        memmap = np.memmap(memmap_path, dtype=np.float64, mode="r", offset=start_mem, shape=column_shape)
        new_array[:] = memmap[:]
        del memmap
        return new_array

    def write_molecule_metadata_to_disk(self):
        """
        To be changed..
        :return:
        :rtype:
        """
        for molecule_id, unique_col_id, index_in_col in self.yield_molecule_info_in_database():
            molecule = Molecule(id=molecule_id, col_id=unique_col_id, index_in_col=index_in_col)
            self.db.add(molecule)
            self.db.commit()

    def _get_number_of_molecules_in_run(self, run):
        n = 0
        for scan_id in self.db_runs[run]:
            scan_info = self.db_runs[run][scan_id]
            for col_id in scan_info["columns"]:
                col_info = scan_info["columns"][col_id]
                n += int(col_info["molecule_count"])
        return n

    def _get_number_of_molecules_in_db(self):
        n = 0
        for run_id in self.db_runs:
            n += self._get_number_of_molecules_in_run(run_id)
        return n


    def get_single_molecule_from_database(self, molecule_id: int):
        """
        To be changed..
        :param molecule_id:
        :type molecule_id:
        :return:
        :rtype:
        """
        molecule_info = self.db.query(Molecule).get(molecule_id)
        column_id = molecule_info.col_id
        index_in_column = molecule_info.index_in_col
        column_info = self.db.query(ColArray).get(column_id)
        scan_info = self.db.query(ScanLinage).get(column_info.scan_id)
        memmap_path = scan_info.mem_path
        assert column_info.memmap_start != 0
        start_mem = column_info.memmap_start
        end_mem = column_info.memmap_end
        assert (end_mem - start_mem) % np.float64().itemsize == 0
        column_shape = int((end_mem - start_mem) / np.float64().itemsize)
        column_data = np.zeros(column_shape, dtype=np.float64)
        memmap = np.memmap(memmap_path, dtype=np.float64, mode="r", offset=start_mem, shape=column_shape)
        column_data[:] = memmap[:]
        del memmap
        column_data = column_data.reshape((column_info.molecule_count * 2, column_info.longest_shape))
        molecule_length = int(column_data[index_in_column*2 + 1, -1])
        return column_data[index_in_column*2, :molecule_length], column_data[index_in_column*2 + 1, :molecule_length]


def fix_directory_end(path: str):
    if not path.endswith("/"):
        return path + "/"
    else:
        return path


def molecules_to_bnx(molecule_stream: types.GeneratorType, zoom_ratio: float, final_ratio: float, bnx_filename: str,
                     bnx_template_path="bnx_head.txt", filter_function=lambda x, y: (x, y),
                     signal_to_noise_ratio=3) -> None:
    """
    Converts molecules to BNX file with a given signal to noise ratio and filter function for removing unwanted
    signals.
    Parameters
    ----------
    molecule_stream :       Stream of molecule signal pairs from the memmap.
                             See MoleculeConnector.yield_molecule_signals..
    zoom_ratio :            This is the zoom factor to be applied before label detection. Improves the method.
                             Suggested to be set to 10.
    final_ratio :           The final zoom factor which indicates bases per pixel. This is usually around 480 - 520.
    bnx_filename :          Output file name.
    bnx_template_path :     Path to Template for BNX headers. Default is set to the standard template provided by
                             the OptiScan tool.
    filter_function :       Any function which filters molecule pairs and returns None or modified molecule pairs.
    signal_to_noise_ratio : Signal to noise ratio filter for peak finding method.

    Returns
    -------
    Bnx file to disk.
    """
    bnx_head_file = open(bnx_template_path, "r")
    bnx_head = bnx_head_file.read()
    bnx_head_file.close()
    out_file = open(bnx_filename, "w")
    out_file.write(bnx_head)
    mol_id = 1
    for nick_signal, backbone_signal in molecule_stream:
        nick_signal, backbone_signal = filter_function(nick_signal, backbone_signal)
        assert nick_signal.shape == backbone_signal.shape
        if nick_signal.shape[0] == 0 or backbone_signal.shape[0] == 0:
            continue
        bnx_info = get_bnx_info_with_upsampling(nick_signal, signal_to_noise_ratio, zoom_ratio, final_ratio)
        if len(bnx_info["nick_snrs"]) > 1:
            backbone_signal_avg = np.mean(backbone_signal)
            first_line = "0\t" + "\t".join([str(mol_id), str(backbone_signal.shape[0]*final_ratio),
                                            str(round(backbone_signal_avg, 1)), str(round(backbone_signal_avg, 1)),
                                            str(len(bnx_info["nick_distances"])), str(mol_id), "1", "-1",
                                            "20249,12205,10/27/2015,850020394", "1", "1", "1"]) + "\n"
            second_line = "1\t" + "\t".join([str(round(x, 1)) for x in bnx_info["nick_distances"]]) + "\n"
            third_line = "QX11\t" + "\t".join([str(round(x, 1)) for x in bnx_info["nick_snrs"]]) + "\n"
            fourth_line = "QX12\t" + "\t".join([str(round(x, 1)) for x in bnx_info["nick_intensities"]]) + "\n"
            out_file.write(first_line + second_line + third_line + fourth_line)
            mol_id += 1
        else:
            continue
    out_file.close()
