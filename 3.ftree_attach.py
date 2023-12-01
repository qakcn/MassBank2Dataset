# This is the 3rd step of generating data set from MassBank database.
# This script attaches the fragment tree json to the data.
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
from multiprocessing import freeze_support
from pathlib import Path
# from hashlib import sha256

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from classes import *
from classes.FragmentTree import process_ftree_file

##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
input_path = Path("inputs")
output_path = Path("outputs")
intermediate_path = Path("intermediates")

# Ftree input path
ftrees_input_path = input_path / "ftrees"

# ftree slices path
spectra_slices_path = intermediate_path / "spectra_slices"

# Files
spectra_peaked_named_file = intermediate_path / "spectra.peaked.named.pkl"
spectra_ftrees_file = intermediate_path / "spectra.ftrees.pkl"

slice_num = 12

##################################################
# End of parameters, do not edit below this line #
##################################################

# Initialization
TS.register_tqdm(tqdm)
is_spectra_peaked_named = is_spectra_ftrees = False


if spectra_ftrees_file.is_file():
    TS.ip("Loading ftrees spectra from pickle file...")
    spectra = pd.read_pickle(spectra_ftrees_file)
    is_spectra_peaked_named = is_spectra_ftrees = True
    TS.p(TS.green("Done."))
elif spectra_peaked_named_file.is_file():
    TS.ip("Loading peaked & named spectra from pickle file...")
    spectra = pd.read_pickle(spectra_peaked_named_file)
    is_spectra_peaked_named = True
    TS.p(TS.green("Done."))
    

if not is_spectra_peaked_named:
    TS.p(TS.red(TS.bold("Error: peaked & named spectra pickle file not found.")))
elif not is_spectra_ftrees:
    TS.ip("Getting ftree numbers...")
    ftree_num = dict(
        orbitrap = len(list((ftrees_input_path/"orbitrap").glob("*.json"))),
        qtof = len(list((ftrees_input_path/"qtof").glob("*.json"))),
    )
    TS.p(f"orbitrap: {TS.cyan(ftree_num['orbitrap'])}, qtof: {TS.cyan(ftree_num['qtof'])}")

    tqdm_bars = dict(
        orbitrap = tqdm(desc="orbit ftrees", total=ftree_num["orbitrap"], ncols=100, position=1),
        qtof = tqdm(desc="qtof ftrees", total=ftree_num["qtof"], ncols=100, position=2),
    )

    tqdm.pandas(desc="Reading ftree files", total=spectra.shape[0], ncols=100, position=0)
    spectra["ftree"] = spectra.progress_apply(process_ftree_file, axis=1, ftrees_path=ftrees_input_path, tqdm_bars=tqdm_bars)
    tqdm_bars["orbitrap"].close()
    tqdm_bars["qtof"].close()

    TS.ip("Dropping spectra without ftrees...")
    spectra.dropna(subset=["ftree"], inplace=True)
    TS.p(TS.green("Done")+f", {TS.cyan(spectra.shape[0])} spectra left.")

    TS.ip("Saving spectra ftrees to pickle file...")
    spectra.to_pickle(spectra_ftrees_file)
    TS.p(TS.green("Done."))
    
    is_spectra_ftrees = True

# if is_spectra_ftrees:
#     total = spectra.shape[0]
#     slice_size = total // slice_num
#     spectra_slices_path.mkdir(parents=True, exist_ok=True)
#     TS.p(f"Splitting spectra into {slice_num} slices...")
#     for slice in range(slice_num):
#         start = slice * slice_size
#         end = (slice + 1) * slice_size
#         if slice == slice_num - 1:
#             end = total
#         TS.ip(f"    {TS.cyan(slice)}: {TS.yellow(start)}-{TS.yellow(end)} total: {TS.magenta(end-start)}...")
#         spectra.iloc[start:end].to_pickle(spectra_slices_path / f"spectra.ftrees.{slice}.pkl")
#         TS.p(TS.green("Done."))
#     TS.p(TS.green("All done."))
