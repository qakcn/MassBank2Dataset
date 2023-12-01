# This is the 6th step of generating data set from MassBank database.
# This script generates molecular fingerprints for the compounds corresponding 
# to the mass spectra.
#
# Author: qakcn
# Email: qakcn@hotmail.com
# Date: 2023-11-29

# Copyright 2023 qakcn
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# PSL imports
from pathlib import Path
from multiprocessing import Process, Queue, freeze_support

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from classes import *
from classes.Functions import df_add_columns
from classes.Fingerprint import fingerprint_row, start_jvm, shutdown_jvm

class Fingerprint(Process):
    def __init__(self, spectra_slice, slice: int, slice_total,  fp_slices_path, print_queue: Queue):
        super().__init__()
        self.spectra_slice = spectra_slice
        self.slice = slice
        self.slice_total = slice_total
        self.fp_slices_path = fp_slices_path
        self.print_queue = print_queue

    def run(self):
        self.print_queue.put(("start", {'slice': self.slice}))

        self.print_queue.put(("start_fp", {'slice': self.slice, 'slice_total': self.slice_total}))
        start_jvm()
        fp_data = self.spectra_slice.apply(fingerprint_row, axis=1, slice=self.slice, print_queue=self.print_queue)
        shutdown_jvm()
        self.print_queue.put(("end_fp", {'slice': self.slice}))

        self.print_queue.put(("saving", {'slice': self.slice}))
        fp_data = fp_data.to_pickle(self.fp_slices_path / f"{self.slice}.pkl")
        self.print_queue.put(("saved", {'slice': self.slice}))

        self.print_queue.put(("end", {'slice': self.slice}))

def info_printer(print_queue: Queue):
    templates = {
        "slice": "Slice " + TS.blue("{slice}") + ": ",
        "start": "Start slice.",
        "start_fp": "Start fingerprinting spectra.",
        "end_fp": "End fingerprinting spectra.",
        "saving": "Saving spectra file.",
        "saved": "Saved spectra file.",
        "end": "End slice.",
    }
    tqdm_bars = {}
    TS.register_tqdm(tqdm)
    while True:
        msgtype, msgdata = print_queue.get()
        if msgtype == "end_all":
            break
        if msgtype in templates:
            TS.p(templates["slice"].format(**msgdata) + templates[msgtype])
            if msgtype == "start_fp":
                tqdm_bars[msgdata['slice']] = tqdm(desc=f"Slice {TS.blue(msgdata['slice'])}: fingerprinting", total=msgdata['slice_total'], ncols=100, position=msgdata['slice'])
            if msgtype == "end_fp":
                tqdm_bars[msgdata["slice"]].close()
        elif msgtype == "progress":
            if msgdata["slice"] in tqdm_bars:
                tqdm_bars[msgdata["slice"]].update()
            else:
                print_queue.put((msgtype, msgdata))
            

if __name__ == "__main__":
    freeze_support()

    ##################################################
    # Parameters that can be edited by the user      #
    ##################################################
    # Paths
    output_path = Path("outputs")
    intermediate_path = Path("intermediates")

    fp_slices_path = intermediate_path / "fp_slices"

    # Files
    spectra_tested_file = intermediate_path / "spectra.tested.pkl"
    spectra_fp_file = intermediate_path / "spectra.fp.pkl"

    slice_num = 12

    ##################################################
    # End of parameters, do not edit below this line #
    ##################################################

    is_spectra_fped = False

    TS.register_tqdm(tqdm)

    if len(list(fp_slices_path.glob("*.pkl"))) == slice_num:
        is_spectra_fped = True
    if spectra_tested_file.is_file():
        TS.ip(f"Read tested spectra file...")
        spectra = pd.read_pickle(spectra_tested_file)
        TS.p(TS.green(f"Done."))
    

    if not is_spectra_fped:

        fp_slices_path.mkdir(parents=True, exist_ok=True)

        procs = []
        print_queue = Queue()

        total = spectra.shape[0]
        slice_size = total // slice_num

        for slice in range(slice_num):
            start_loc = slice * slice_size
            end_loc = start_loc + slice_size
            if slice == slice_num - 1:
                end_loc = total
            slice_total = end_loc - start_loc

            spectra_slice = spectra.iloc[start_loc:end_loc].copy()
            p = Fingerprint(spectra_slice, slice, slice_total, fp_slices_path, print_queue)
            p.start()
            procs.append(p)
        
        print_p = Process(target=info_printer, args=(print_queue,))
        print_p.start()

        for p in procs:
            p.join()
        
        print_queue.put(("end_all", {}))
        print_p.join()

        is_spectra_fped = True

        TS.p(TS.green("All slices fingerprinted."))

    if is_spectra_fped:
        TS.ip(f"Concatenate fingerprints...")
        columns = [
            "FP2",
            "AtomPair",
            "Avalon",
            "MACCS",
            "Morgan",
            "TopologicalTorsion",
            "RDKitFingerprint",
            "CDKFingerprint",
            "PubChemFingerprint",
            "Klekota-Roth"
        ]

        df_add_columns(spectra, columns)

        TS.p(TS.green("Done."))

        for slice in range(slice_num):
            TS.ip(f"Read slice {slice}...")
            fp_slice = pd.read_pickle(fp_slices_path / f"{slice}.pkl")
            TS.p(TS.green("Done."))

            TS.ip(f"Update spectra for slice {slice}...")
            spectra.update(fp_slice)
            TS.p(TS.green("Done."))

        TS.ip(f"Save fingerprinted spectra...")
        spectra.to_pickle(spectra_fp_file)
        TS.p(TS.green("Done."))
