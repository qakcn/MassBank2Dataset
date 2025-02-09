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
from classes.Fingerprint import standardize_row


##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
input_path = Path("inputs")
output_path = Path("outputs")
intermediate_path = Path("intermediates")

fp_slices_path = intermediate_path / "fp_slices"
casmi_path = input_path / "casmi2016"

# Files
casmi_file = intermediate_path / "casmi.pkl"
spectra_fp_file = intermediate_path / "spectra.fp.pkl"
spectra_fp_standardized_file = intermediate_path / "spectra.fp.standardized.pkl"
spectra_fp_without_casmi_file = intermediate_path / "spectra.fp.no_casmi.pkl"
spectra_fp_within_casmi_file = intermediate_path / "spectra.fp.casmi.pkl"

##################################################
# End of parameters, do not edit below this line #
##################################################

is_spectra_standardized = is_casmi = False
TS.register_printer(tqdm.write)

if casmi_file.is_file():
    TS.ip(f"Read CASMI file...")
    casmi = pd.read_pickle(casmi_file)
    is_casmi = True
    TS.p(TS.green(f"Done."))
if spectra_fp_standardized_file.is_file():
    TS.ip(f"Read standardized spectra file...")
    spectra = pd.read_pickle(spectra_fp_standardized_file)
    is_spectra_standardized = True
    TS.p(TS.green(f"Done."))

if not is_casmi:
    TS.ip(f"Read CASMI traning set file...")
    casmi = pd.DataFrame(columns=["compound", "formula", "parentmass", "inchikey", "inchi", "smiles", "instrument", "sign"])

    for sign in ("positive", "negative"):
        for file in (casmi_path / sign).glob("*.ms"):
            msdata = {"sign": sign}
            with open(file, "r") as f:
                for line in f:
                    if line.startswith(">compound"):
                        msdata["compound"] = line.split(" ", 1)[1].strip()
                    elif line.startswith("#formula"):
                        msdata["formula"] = line.split(" ", 1)[1].strip()
                    elif line.startswith(">parentmass"):
                        msdata["parentmass"] = line.split(" ", 1)[1].strip()
                    elif line.startswith("#inchikey"):
                        msdata["inchikey"] = line.split(" ", 1)[1].strip()
                    elif line.startswith("#inchi"):
                        msdata["inchi"] = line.split(" ", 1)[1].strip()
                    elif line.startswith("#smiles"):
                        msdata["smiles"] = line.split(" ", 1)[1].strip()
                    elif line.startswith(">instrument"):
                        msdata["instrument"] = line.split(" ", 1)[1].strip()
            casmi.loc[casmi.shape[0]] = msdata

    TS.p(TS.green(f"Done."))

    tqdm.pandas(desc="Standardizing SMILES", total=casmi.shape[0], ncols=100)
    casmi["can_smiles"] = casmi.progress_apply(standardize_row, axis=1, molstr_type="inchi")
    TS.p(TS.green(f"Done standardizing SMILES."))

    TS.ip(f"Saving CSI:FingerID file...")
    casmi.to_pickle(casmi_file)
    TS.p(TS.green(f"Done."))

    is_casmi = True

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

if is_casmi and is_spectra_standardized:
    TS.ip(f"Splitting CSI:FingerID training set...")
    mask = spectra["can_smiles"].isin(casmi["can_smiles"])
    spectra_within_csi = spectra[mask]
    spectra_without_csi = spectra[~mask]
    TS.p(TS.green(f"Done."))

    TS.ip(f"Saving spectra split files...")
    spectra_within_csi.to_pickle(spectra_fp_within_casmi_file)
    spectra_without_csi.to_pickle(spectra_fp_without_casmi_file)
    TS.p(TS.green(f"Done."))
