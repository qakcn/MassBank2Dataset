# This is the 1st step of generating data set from MassBank database.
# This script query data from MassBank database and generate .ms file.
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
import sqlalchemy as sa
from tqdm import tqdm

# Local imports
from classes import *
from classes.Functions import peak_spectra_row, name_spectra_row, generate_ms_file_row

##################################################
# Parameters that can be edited by the user      #
##################################################
## Paths
input_path = Path("inputs")
output_path = Path("outputs")
intermediate_path = Path("intermediates")

ms_output_path = output_path / "ms"

## Files
spectra_file = intermediate_path / "spectra.pkl"
spectra_peaked_file = intermediate_path / "spectra.peaked.pkl"
spectra_peaked_named_file = intermediate_path / "spectra.peaked.named.pkl"

## Database
sa_mysql_str = "mysql+mysqldb://username:password@localhost:3306/massbank"
##################################################
# End of parameters, do not edit below this line #
##################################################

# Initialization
is_spectra = is_spectra_peaked = is_spectra_peaked_named = False
TS.register_tqdm(tqdm)

if spectra_peaked_named_file.is_file():
    TS.ip("Loading peaked & named spectra from pickle file...")
    spectra = pd.read_pickle(spectra_peaked_named_file)
    is_spectra = is_spectra_peaked = is_spectra_peaked_named = True
    TS.p(TS.green("Done"))
elif spectra_peaked_file.is_file():
    TS.ip("Loading peaked spectra from pickle file...")
    spectra = pd.read_pickle(spectra_peaked_file)
    is_spectra = is_spectra_peaked = True
    TS.p(TS.green("Done"))
elif spectra_file.is_file():
    TS.ip("Loading spectra from pickle file...")
    spectra = pd.read_pickle(spectra_file)
    is_spectra = True
    TS.p(TS.green("Done"))


if not is_spectra_peaked_named:
    engine = sa.create_engine(sa_mysql_str)

    with engine.connect() as conn:
        if not is_spectra:
            TS.ip("Loading spectra from database...")
            qs =  r"SELECT a.*,b.* FROM `msms_spectrum` AS a JOIN `ms_compound` AS b ON a.`compound_id`=b.`compound_id` WHERE a.`ms_level`='MS2'"
            spectra = pd.read_sql_query(qs, conn)
            TS.p(TS.green("Done"))

            TS.ip("Saving spectra to pickle file...")
            spectra.to_pickle(spectra_file)
            TS.p(TS.green("Done"))

            is_spectra = True

        if not is_spectra_peaked:
            spectra_peaks_file = intermediate_path / "spectra_peaks.pkl"

            if spectra_peaks_file.is_file():
                TS.ip("Loading spectra peaks from pickle file...")
                spectra_peaks = pd.read_pickle(spectra_peaks_file)
                TS.p(TS.green("Done"))
            else:
                TS.ip("Loading spectra peaks from database...")
                qs = r"SELECT `spectrum_id`,`mz`,`intensity` FROM `msms_spectrum_peak`"
                spectra_peaks = pd.read_sql_query(qs, conn)
                TS.p(TS.green("Done"))

                TS.ip("Saving spectrum peaks to pickle file...")
                spectra_peaks.to_pickle(spectra_peaks_file)
                TS.p(TS.green("Done"))

            tqdm.pandas(desc="Peaking spectra", total=spectra.shape[0], ncols=100)
            spectra["peaks"] = spectra.progress_apply(peak_spectra_row, axis=1, spectra_peaks=spectra_peaks)

            TS.ip("Saving peaked spectra to pickle file...")
            spectra.to_pickle(spectra_peaked_file)
            TS.p(TS.green("Done"))

            is_spectra_peaked = True

        
        # if not is_spectra_peaked_named:
        compound_names_file = intermediate_path / "compound_names.pkl"

        if compound_names_file.is_file():
            TS.ip("Loading compound names from pickle file...")
            compound_names = pd.read_pickle(compound_names_file)
            TS.p(TS.green("Done"))
        else:
            TS.ip("Loading compound names from database...")
            qs = r"SELECT * FROM `compound_name` JOIN `name` ON `compound_name`.`NAME`=`name`.`ID`"
            compound_names = pd.read_sql_query(qs, conn)
            TS.p(TS.green("Done"))

            TS.ip("Saving compound names to pickle file...")
            compound_names.to_pickle(compound_names_file)
            TS.p(TS.green("Done"))
        
        tqdm.pandas(desc="Naming compound", total=spectra.shape[0], ncols=100)
        spectra["compound_name"] = spectra.progress_apply(name_spectra_row, axis=1, compound_names=compound_names)

        TS.ip("Saving named spectra to pickle file...")
        spectra.to_pickle(spectra_peaked_named_file)
        TS.p(TS.green("Done"))

        is_spectra_peaked_named = True

if is_spectra_peaked_named:
    TS.p(TS.green(TS.bold("Spectra are ready for generation of MS files.")))
    ms_output_path.mkdir(parents=True, exist_ok=True)
    (ms_output_path / "qtof").mkdir(parents=True, exist_ok=True)
    (ms_output_path / "orbitrap").mkdir(parents=True, exist_ok=True)
    tqdm.pandas(desc="Generating MS files", total=spectra.shape[0], ncols=100)
    spectra.progress_apply(generate_ms_file_row, axis=1, output_path=ms_output_path)