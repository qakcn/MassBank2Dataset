# This script contains functions functioning in other scripts.
#
# Author: qakcn
# Email: qakcn@hotmail.com
# Date: 2023-12-01

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


if __name__ == "__main__":
    raise SystemExit("This script is not meant to be run directly")

# PSL imports
from pathlib import Path
from hashlib import sha256
from typing import List, Union

# Third party imports
import pandas as pd

def peak_spectra_row(row, spectra_peaks: pd.DataFrame):
    spectrum_id = row["spectrum_id"]
    peaks = spectra_peaks.loc[spectra_peaks["spectrum_id"] == spectrum_id, ["mz", "intensity"]]
    peaks_dict = peaks.to_dict(orient="split")
    return peaks_dict["data"]

def name_spectra_row(row, compound_names: pd.DataFrame):
    compound_id = row["compound_id"]
    names = compound_names.loc[compound_names["COMPOUND"] == compound_id, ["CH_NAME"]]
    names_dict = names.to_dict(orient="split")
    return [i for j in names_dict["data"] for i in j]

def generate_ms_file_row(row, output_path: Path):
    spectrum_id = row["spectrum_id"]
    compound_id = row["compound_id"]
    
    metadata = dict(
        compound = row["compound_name"][0],
        formula = row["formula"],
        instrumentation = f"{row['instrument']} ({row['instrument_type']})",
        parentmass = row["precursor_mz_text"],
    )

    if row["adduct"] is not None:
        metadata["ionization"] = row["adduct"].replace("[", r"\[").replace("]", r"\]")

    comments = dict(
        spectrum_id = spectrum_id,
        compound_id = compound_id,
        compound_name = ", ".join(row["compound_name"]),
        profile = inst_profile(row["instrument_type"], row["instrument"])
    )

    peaks = row["peaks"]

    ms_lines = []

    for key, value in metadata.items():
        ms_lines.append(f">{key} {value}")
    
    ms_lines.append("")

    for key, value in comments.items():
        ms_lines.append(f"#{key} {value}")
    
    ms_lines.append("")

    ms_lines.append(f">ms2")

    for mz, intensity in peaks:
        ms_lines.append(f"{mz} {intensity}")
    
    ms_content = "\n".join(ms_lines)
    ms_file = output_path / comments["profile"] / f"{spectrum_id}_{compound_id}.ms"

    with open(ms_file, "w") as f:
        f.write(ms_content)
    
def inst_profile(inst_type: str, inst: str) -> str:
    p2i = {
        "qtof": ["ESI-QTOF", "LC-ESI-QTOF", "GC-APCI-QTOF", "LC-ESI-QIT", "LC-ESI-QQ", "ESI-TOF", "MALDI-TOFTOF", "LC-ESI-IT", "FAB-EBEB", "LC-ESI-Q", "LC-ESI-TOF", "LC-ESI-ITTOF", "ESI-ITTOF", "MALDI-QIT", "APCI-ITTOF", "MALDI-QITTOF", "ESI-QQ", "GC-EI-QQ", "LC-APPI-QQ", "LC-ESI-QQQ", "ESI-QIT", "ESI-ITFT"],
        "orbitrap": ["LC-ESI-QFT", "LC-ESI-ITFT", "LC-APCI-QFT", "LC-ESI-FT", "LC-APCI-ITFT", "APCI-ITFT"],
        }

    i2p = {i: j for j in p2i for i in p2i[j]}

    if inst_type == "ESI-ITFT" and inst.find("FT-ICP") != -1:
        return "fticr"
    else:
        return i2p[inst_type]

def df_add_columns(df: pd.DataFrame, columns: List[str]):
    for column in columns:
        if column not in df.columns:
            df[column] = None

def get_hashkey(data: Union[str, List[str]]) -> str:
    str = "|".join(data) if isinstance(data, list) else data
    hasher = sha256()
    hasher.update(str.encode("utf-8"))
    return hasher.hexdigest()


    