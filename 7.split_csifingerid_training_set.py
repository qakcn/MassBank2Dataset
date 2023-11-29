# This is the 7th step of generating data set from MassBank database.
# This script splits the data according to the training set of CSI:FingerID.
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

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from classes import *



##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
input_path = Path("inputs")
output_path = Path("outputs")
intermediate_path = Path("intermediates")

fp_slices_path = intermediate_path / "fp_slices"

# Files
csi_fingerid_file = intermediate_path / "csi_fid.pkl"
spectra_fp_file = intermediate_path / "spectra.fp.pkl"
spectra_fp_standardized_file = intermediate_path / "spectra.fp.standardized.pkl"
spectra_fp_without_csi_fingerid_file = intermediate_path / "spectra.fp.no_csi_fid.pkl"
spectra_fp_within_csi_fingerid_file = intermediate_path / "spectra.fp.csi_fid.pkl"

##################################################
# End of parameters, do not edit below this line #
##################################################

is_spectra_standardized = is_csi_fingerid = False
TS.register_tqdm(tqdm)

if csi_fingerid_file.is_file():
    TS.ip(f"Read CSI:FingerID file...")
    csi_fingerid = pd.read_pickle(csi_fingerid_file)
    is_csi_fingerid = True
    TS.p(TS.green(f"Done."))
if spectra_fp_standardized_file.is_file():
    TS.ip(f"Read standardized spectra file...")
    spectra = pd.read_pickle(spectra_fp_standardized_file)
    is_spectra_standardized = True
    TS.p(TS.green(f"Done."))

if not is_csi_fingerid:
    TS.ip(f"Read CSI:FingerID traning set file...")
    pos_csi = pd.read_table(input_path / "csi_fingerid-trainingstructures-positive", sep="\t", names=["inchikey", "inchi"])
    neg_csi = pd.read_table(input_path / "csi_fingerid-trainingstructures-positive", sep="\t",  names=["inchikey", "inchi"])
    TS.p(TS.green(f"Done."))

    pos_csi["sign"] = "positive"
    neg_csi["sign"] = "negative"
    csi_fingerid = pd.concat([pos_csi, neg_csi], ignore_index=True)
    del pos_csi, neg_csi

    TS.ip(f"Fixing InChI strings...")
    fix_list = [
        (r"InChI=1S/C22H35NO2/c1-14-7-6-9-17(23(14)3)11-12-19-18-10-5-4-8-16(18)13-20-21(19)15(2)25-22(20)24/h11-12,14", r"InChI=1S/C22H35NO2/c1-14-7-6-9-17(23(14)3)11-12-19-18-10-5-4-8-16(18)13-20-21(19)15(2)25-22(20)24/h11-12,14-21H,4-10,13H2,1-3H3")
    ]
    for old, new in fix_list:
        csi_fingerid["inchi"] = csi_fingerid["inchi"].str.replace(old, new)
    TS.p(TS.green(f"Done."))

    tqdm.pandas(desc="Standardizing SMILES", total=csi_fingerid.shape[0], ncols=100)
    csi_fingerid["can_smiles"] = csi_fingerid.progress_apply(standardize_row, axis=1, molstr_type="inchi")
    TS.p(TS.green(f"Done standardizing SMILES."))

    TS.ip(f"Saving CSI:FingerID file...")
    csi_fingerid.to_pickle(csi_fingerid_file)
    TS.p(TS.green(f"Done."))

    is_csi_fingerid = True

if not is_spectra_standardized:
    TS.ip(f"Read spectra file...")
    spectra = pd.read_pickle(spectra_fp_file)
    TS.p(TS.green(f"Done."))

    tqdm.pandas(desc="Standardizing SMILES", total=spectra.shape[0], ncols=100)
    spectra["can_smiles"] = spectra.progress_apply(standardize_row, axis=1, molstr_type="smiles")
    TS.p(TS.green(f"Done standardizing SMILES."))

    TS.ip(f"Saving standardized spectra file...")
    spectra.to_pickle(spectra_fp_standardized_file)
    TS.p(TS.green(f"Done."))

    is_spectra_standardized = True

if is_csi_fingerid and is_spectra_standardized:
    TS.ip(f"Splitting CSI:FingerID training set...")
    mask = spectra["can_smiles"].isin(csi_fingerid["can_smiles"])
    spectra_within_csi = spectra[mask]
    spectra_without_csi = spectra[~mask]
    TS.p(TS.green(f"Done."))

    TS.ip(f"Saving spectra split files...")
    spectra_within_csi.to_pickle(spectra_fp_within_csi_fingerid_file)
    spectra_without_csi.to_pickle(spectra_fp_without_csi_fingerid_file)
    TS.p(TS.green(f"Done."))
