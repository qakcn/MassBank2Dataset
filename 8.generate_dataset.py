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

##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
output_path = Path("outputs")
intermediate_path = Path("intermediates")

# Files
spectra_fp_without_csi_fingerid_file = intermediate_path / "spectra.fp.no_csi_fid.pkl"
spectra_fp_within_csi_fingerid_file = intermediate_path / "spectra.fp.csi_fid.pkl"
counter_file = intermediate_path / "counter.nopd.pkl"

datasets_without_csi_file = output_path / "datasets.no_csi_fid.pkl"
datasets_within_csi_file = output_path / "datasets.csi_fid.pkl"

##################################################
# End of parameters, do not edit below this line #
##################################################

TS.register_tqdm(tqdm)

elemset = set()

TS.ip("Reading files...")
spectra_no_csi = pd.read_pickle(spectra_fp_without_csi_fingerid_file)
spectra_csi = pd.read_pickle(spectra_fp_within_csi_fingerid_file)
with open(counter_file, "rb") as f:
    counter = pickle.load(f)
TS.p(TS.green("Done."))

TS.p("Calculating for spectra within CSI:FingerID...")
tqdm.pandas(desc=f"    Calculating", total=spectra_csi.shape[0], ncols=100)
datasets_within_csi = spectra_csi.progress_apply(calc_dataset_row, axis=1, counter=counter, elemset=elemset)
TS.p(TS.green("Done."))

TS.p("Calculating for spectra without CSI:FingerID...")
tqdm.pandas(desc="    Calculating", total=spectra_no_csi.shape[0], ncols=100)
datasets_without_csi = spectra_no_csi.progress_apply(calc_dataset_row, axis=1, counter=counter, elemset=elemset)
TS.p(TS.green("Done."))

TS.ip("Saving datasets...")
with open(output_path/"element.nopd.pkl", "wb") as f:
    pickle.dump(elemset, f)
datasets_within_csi.to_pickle(datasets_within_csi_file)
datasets_without_csi.to_pickle(datasets_without_csi_file)
TS.p(TS.green("Done."))
