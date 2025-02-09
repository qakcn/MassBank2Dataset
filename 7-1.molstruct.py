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
from classes.Fingerprint import standardize_row, gen_molstruct_row


##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
input_path = Path("inputs")
output_path = Path("outputs")
intermediate_path = Path("intermediates")

# Files
spectra_fp_file = intermediate_path / "spectra.fp.pkl"
spectra_fp_standardized_file = intermediate_path / "spectra.fp.standardized.pkl"

spectra_fp_molstruct_file = output_path / "dataset.molstruct.pkl"

##################################################
# End of parameters, do not edit below this line #
##################################################

is_spectra_standardized  = False
TS.register_printer(tqdm.write)

if spectra_fp_standardized_file.is_file():
    TS.ip(f"Read standardized spectra file...")
    spectra = pd.read_pickle(spectra_fp_standardized_file)
    is_spectra_standardized = True
    TS.p(TS.green(f"Done."))

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

if is_spectra_standardized:
    tqdm.pandas(desc="Generating molecular structures", total=spectra.shape[0], ncols=100)
    spectra["mol_struct"] = spectra.progress_apply(gen_molstruct_row, axis=1)
    TS.p(TS.green(f"Done generating molecular structures."))

    spectra = spectra.drop_duplicates(subset=["can_smiles"], keep="first")

    spectra = spectra[["spectrum_id", "pubchem", "can_smiles", "FP2", "AtomPair", "Avalon", "MACCS", "Morgan", "TopologicalTorsion", "RDKitFingerprint", "CDKFingerprint", "PubChemFingerprint", "Klekota-Roth", "mol_struct"]]

    TS.ip(f"Saving molecular structures file...")
    spectra.to_pickle(spectra_fp_molstruct_file)
    TS.p(TS.green(f"Done."))
    