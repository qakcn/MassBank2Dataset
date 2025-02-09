# This is the 8th step of generating data set from MassBank database.
# This script calculates the feature values and generates the data set file.
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
import pickle

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from classes import *
from classes.Fingerprint import calc_dataset_row

##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
output_path = Path("outputs")
intermediate_path = Path("intermediates")

# Files
spectra_fp_without_casmi_file = intermediate_path / "spectra.fp.no_casmi.pkl"
spectra_fp_within_casmi_file = intermediate_path / "spectra.fp.casmi.pkl"
counter_file = intermediate_path / "counter.nopd.pkl"

datasets_without_casmi_file = output_path / "dataset.no_casmi.pkl"
datasets_within_casmi_file = output_path / "dataset.casmi.pkl"

##################################################
# End of parameters, do not edit below this line #
##################################################

TS.register_printer(tqdm.write)

elemset = set()

TS.ip("Reading files...")
spectra_no_casmi = pd.read_pickle(spectra_fp_without_casmi_file)
spectra_casmi = pd.read_pickle(spectra_fp_within_casmi_file)
with open(counter_file, "rb") as f:
    counter = pickle.load(f)
TS.p(TS.green("Done."))

TS.p("Calculating for spectra within CASMI...")
tqdm.pandas(desc=f"    Calculating", total=spectra_casmi.shape[0], ncols=100)
dataset_within_casmi = spectra_casmi.progress_apply(calc_dataset_row, axis=1, counter=counter, elemset=elemset)
TS.p(TS.green("Done."))

TS.p("Calculating for spectra without CASMI...")
tqdm.pandas(desc="    Calculating", total=spectra_no_casmi.shape[0], ncols=100)
dataset_without_casmi = spectra_no_casmi.progress_apply(calc_dataset_row, axis=1, counter=counter, elemset=elemset)
TS.p(TS.green("Done."))

TS.ip("Saving datasets...")
with open(output_path/"element_casmi.nopd.pkl", "wb") as f:
    pickle.dump(elemset, f)
dataset_within_casmi.to_pickle(datasets_within_casmi_file)
dataset_without_casmi.to_pickle(datasets_without_casmi_file)
TS.p(TS.green("Done."))
