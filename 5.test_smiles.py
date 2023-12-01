# This is the 5th step of generating data set from MassBank database.
# This script tests whether the SMILES strings are valid or not.
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
from classes.Functions import df_add_columns
from classes.Fingerprint import test_molstr_row, start_jvm, shutdown_jvm

##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
output_path = Path("outputs")
intermediate_path = Path("intermediates")

# Files
spectra_ftrees_parsed_file = intermediate_path / "spectra.ftrees.parsed.pkl"

spectra_tested_file = intermediate_path / "spectra.tested.pkl"
spectra_fp_file = intermediate_path / "spectra.fp.pkl"

##################################################
# End of parameters, do not edit below this line #
##################################################

TS.register_printer(tqdm.write)

if spectra_ftrees_parsed_file.is_file():
    TS.ip(f"Read ftree parsed spectra file...")
    spectra = pd.read_pickle(spectra_ftrees_parsed_file)
    TS.p(TS.green(f"Done."))


TS.p(f"Test spectra...")

total = spectra.shape[0]

start_jvm()
tqdm.pandas(desc=f"Test", total=total, ncols=100)
new_data = spectra.progress_apply(test_molstr_row, axis=1)
TS.p(TS.green(f"Tested."))
shutdown_jvm()

TS.ip(f"Concatenate test results...")
df_add_columns(spectra, ["test_smile", "test_inchi", "test_results"])
spectra.update(new_data)
TS.p(TS.green("Done."))

spectra_smiles = spectra[spectra["test_smile"].astype(bool)]
spectra_nosmiles = spectra[~spectra["test_smile"].astype(bool)]

TS.ip(f"Save tested spectra...")
spectra_smiles.to_pickle(spectra_tested_file)
spectra_nosmiles.to_pickle(intermediate_path / "spectra.tested.nosmiles.pkl")
TS.p(TS.green("Done."))
