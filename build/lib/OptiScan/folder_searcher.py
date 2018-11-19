from os import walk

__author__ = "MAkdel"


class FolderSearcher:
    """
    This class traverses all the folders below the given main dirctory and extracts relevant bionano files.
     The main directory should correspond to a single run.
    """
    def __init__(self, main_directory):
        """
        Initiates the class. Uses all methods in the class. The scans dictionary contains files specific to each scan.
        The mqr_map_file, scafforld_cmap_file, main_bnx_file and the autofocus file directories are for all scans.
        :param main_directory: Directory of the run with possibly multiple scans; str
        """
        self.dir = main_directory
        self.fasta = str()
        self.file_dirs = list()
        [[self.file_dirs.append(x[0] + "/" + y) for y in x[-1]] for x in walk(main_directory) if x[2] != []]
        self.file_names = [x.split("/")[-1] for x in self.file_dirs]
        self.scans = dict()
        self.autofocus_file = str()

        self.mqr_map_file = str()
        self.scaffold_cmap_file = str()
        self.main_bnx_file = str()
        self.get_central_files()
        self.fill_scans_dictionary()
        self._fill_missing_keys()

    def _search_file_exact(self, file_name):
        """
        Can be used to search a specific filename and return their full paths.
        :param file_name: File name as query; case sensitive; str
        :return: An array of full file paths found for the specified filename;list
        """
        return [self.file_dirs[i] for i in range(len(self.file_names)) if file_name == self.file_names[i]]

    def _fill_missing_keys(self):
        all_keys = set()
        essential_db_keys = ["RawMolecules_file_location", "Molecules_file_location", "Stitch_file_location",
                             "Labels_file_location", "tiff_location"]

        [[all_keys.add(y) for y in self.scans[x].keys()] for x in self.scans.keys()]
        [all_keys.add(x) for x in essential_db_keys]

        for _key in list(all_keys):
            for i in self.scans:
                if _key not in self.scans[i].keys():
                    self.scans[i][_key] = ""

    def _search_file(self, query, include_dirs=False):
        """
        Can be used to search filenames which contains a given query string. Reports back the files with their full
        directories
        :param query: String query which will be searched for; str
        :return: returns an  array of file paths; list
        """
        if not include_dirs:
            return [self.file_dirs[i] for i in range(len(self.file_names)) if query in self.file_names[i]]
        else:
            return [self.file_dirs[i] for i in range(len(self.file_names)) if query in self.file_dirs[i]]

    def _get_scans(self):
        """
        Obtains the scan ids by locating all tiff files which contains scan ids in their filenames. Updates the scan
        dictionary with the info.
        :return:
        """
        tiff_scan_files = self._search_file("Scan")
        for i in range(len(tiff_scan_files))[::-1]:
            scan_name_end_index = tiff_scan_files[i].find("Scan")
            scan_id = tiff_scan_files[i][scan_name_end_index+4:]
            scan_id = scan_id[:scan_id.find(".tiff")]
            try:
                scan_id = int(scan_id)
            except ValueError:
                continue
            self.scans[scan_id] = {"tiff_location": tiff_scan_files[i]}

    def _update_corresponding_id_with_a_file(self, file_prefix, file_extension, files):
        """
        Updates scan dictionary with filepaths corresponding to the scan specific files.
        :param file_prefix: The general file prefix; str
        :param file_extension: The file extension; str
        :param files: Array of files to choose files from; str
        :return:
        """
        prefix_length = len(file_prefix)
        for _file in files:
            scan_name_end_index = _file.find(file_prefix) + prefix_length
            file_extension_index = _file.find(file_extension)
            if not file_extension_index == scan_name_end_index:
                try:
                    scan_id = int(_file[scan_name_end_index:file_extension_index])
                    try:
                        self.scans[scan_id].update({file_prefix+"_file_location": _file})
                    except KeyError:
                        pass
                except ValueError:
                    try:
                        scan_id = int(_file[scan_name_end_index + 10 :file_extension_index])
                    except ValueError:
                        continue
                    try:
                        self.scans[scan_id].update({file_prefix + "_file_location": _file})
                    except KeyError:
                        pass

    def fill_scans_dictionary(self):
        """
        Completely fills the scan dictionary by using methods from the class.
        :return:
        """
        self._get_scans()
        stitch_files = self._search_file("Stitch")
        raw_bnx_files = self._search_file("RawMolecules")
        mol_files = self._search_file(".mol")[1:]
        label_files = self._search_file("Labels")
        self._update_corresponding_id_with_a_file("RawMolecules", ".bnx", raw_bnx_files)
        self._update_corresponding_id_with_a_file("Stitch", ".fov", stitch_files)
        self._update_corresponding_id_with_a_file("Molecules", ".mol", mol_files)
        self._update_corresponding_id_with_a_file("Labels", ".lab", label_files)

    def get_central_files(self):
        """
        Obtains central files. See __init__ for description of the central files.
        :return:
        """
        try:
            self.autofocus_file = self._search_file("AutofocusData.csv")[0]
        except IndexError:
            self.autofocus_file = ""
        try:
            self.scaffold_cmap_file = self._search_file(".fasta.cmap")[0]
        except IndexError:
            self.scaffold_cmap_file = ""
        try:
            self.main_bnx_file = self._search_file_exact("Molecules.bnx")[0]
        except IndexError:
            self.main_bnx_file = ""
        try:
            self.mqr_map_file = self._search_file_exact("MoleculeQualityReport.map")[0]
        except IndexError:
            self.mqr_map_file = ""

if __name__ == "__main__":
    folder = "/Users/akdel/Bionano/onescan/apple/2016-01/apple_run1_fc1_2016-05-17_16_35/"
    sr = FolderSearcher(folder)
    print(sr.scans)
    print(sr.scans.keys())

