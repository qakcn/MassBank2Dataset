# This script contains functions functioning in other scripts.
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


if __name__ == "__main__":
    raise SystemExit("This script is not meant to be run directly")

# PSL imports
import json
from pathlib import Path
from hashlib import sha256
from typing import List, Union

# Third party imports

# Local imports
from .FragmentTree import *
from .TerminalStyle import *
from .FingerprintUtils import *

def peak_spectra_row(row, spectra_peaks):
    spectrum_id = row["spectrum_id"]
    peaks = spectra_peaks.loc[spectra_peaks["spectrum_id"] == spectrum_id, ["mz", "intensity"]]
    peaks_dict = peaks.to_dict(orient="split")
    return peaks_dict["data"]

def name_spectra_row(row, compound_names):
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
    
def inst_profile(inst_type, inst):
    p2i = {
        "qtof": ["ESI-QTOF", "LC-ESI-QTOF", "GC-APCI-QTOF", "LC-ESI-QIT", "LC-ESI-QQ", "ESI-TOF", "MALDI-TOFTOF", "LC-ESI-IT", "FAB-EBEB", "LC-ESI-Q", "LC-ESI-TOF", "LC-ESI-ITTOF", "ESI-ITTOF", "MALDI-QIT", "APCI-ITTOF", "MALDI-QITTOF", "ESI-QQ", "GC-EI-QQ", "LC-APPI-QQ", "LC-ESI-QQQ", "ESI-QIT", "ESI-ITFT"],
        "orbitrap": ["LC-ESI-QFT", "LC-ESI-ITFT", "LC-APCI-QFT", "LC-ESI-FT", "LC-APCI-ITFT", "APCI-ITFT"],
        }

    i2p = {i: j for j in p2i for i in p2i[j]}

    if inst_type == "ESI-ITFT" and inst.find("FT-ICP") != -1:
        return "fticr"
    else:
        return i2p[inst_type]

def process_ftree_file(row, ftrees_path: Path, tqdm_bars):
    spectrum_id = row["spectrum_id"]
    
    profile = inst_profile(row["instrument_type"], row["instrument"])
    profile_path = ftrees_path / profile

    file = list(profile_path.glob(f"*_{spectrum_id}_*.json"))

    if len(file) == 0:
        return None
    
    tqdm_bars[profile].update(1)
    
    with open(file[0], "r") as f:
        ftree_json = json.load(f)
    
    return ftree_json

def parse_ftree_row(row, vocab, counter, orphan_list):
    spectrum_id = row["spectrum_id"]
    compound_id = row["compound_id"]
    compound_name = row["compound_name"][0]

    ftree_json = row["ftree"]

    if len(ftree_json["fragments"]) == 1:
        orphan_list.append(row)
        TS.p(f"Orphan {TS.blue(len(orphan_list))}: {TS.cyan(spectrum_id)}::{TS.yellow(compound_id)}::{TS.magenta(compound_name)}")

        return None
    
    sample_data=row.to_dict()
    del sample_data["ftree"]
    hashkey = get_hashkey(spectrum_id)
    vocab["sample"][hashkey] = sample_data
    del sample_data

    fragments = ftree_json["fragments"]
    losses = ftree_json["losses"]
    del ftree_json["fragments"]
    del ftree_json["losses"]

    compound_columns = ["compound_id", "compound_name", "formula", "exactmass", "smiles", "inchi", "inchikey", "cas", "pubchem"]

    ftree = FragmentTree(spectrum_id, ftree_json, row[compound_columns].to_dict())

    fragment_hashkey = {}

    for fragment in fragments:
        if fragment["id"] == 0:
            ftree.getRootNode().setAttrs(fragment)
        else:
            ftree.addNode(fragment["id"], fragment)
            hashkey = get_hashkey(fragment["molecularFormula"])
            fragment_hashkey[fragment["id"]] = hashkey
            vocab["fragment"][hashkey] = fragment["molecularFormula"]
            ftree.getNode(fragment["id"]).setAttr("hashkey", hashkey)

    for hashkey in set(fragment_hashkey.values()):
        key = (hashkey, "fragment")
        counter[key] = 1 if key not in counter else counter[key]+1

    for loss in losses:
        from_id = loss["source"]
        to_id = loss["target"]
        ftree.addEdge(from_id, to_id, loss)
        if from_id != 0:
            loss_data = [loss["molecularFormula"], fragment_hashkey[from_id], fragment_hashkey[to_id]]
            hashkey = get_hashkey(loss_data)
            vocab["loss"][hashkey] = loss_data
            ftree.getEdge(from_id, to_id).setAttr("hashkey", hashkey)
            key = (hashkey, "loss")
            counter[key] = 1 if key not in counter else counter[key]+1
            key = (fragment_hashkey[from_id], "fragment_in_loss")
            counter[key] = 1 if key not in counter else counter[key]+1
            key = (fragment_hashkey[to_id], "fragment_in_loss")
            counter[key] = 1 if key not in counter else counter[key]+1

            hashkey = get_hashkey("loss_total")
            key = (hashkey, "loss_total")
            counter[key] = 1 if key not in counter else counter[key]+1
    return ftree

def df_add_columns(df: pd.DataFrame, columns: List[str]):
    for column in columns:
        if column not in df.columns:
            df[column] = None

def test_molstr_row(row):
    smiles = row["smiles"]
    inchi = row["inchi"]

    succ_smiles, res_smiles = FingerprintUtils.isSmiles(smiles)
    succ_inchi, res_inchi = FingerprintUtils.isINCHI(inchi)

    return pd.Series({
        "test_smile": succ_smiles,
        "test_inchi": succ_inchi,
        "test_results": {"smiles": res_smiles, "inchi": res_inchi}
    })

def fingerprint_row(row, slice = None, print_queue = None):
    molstr = row["smiles"]
    fptype = [
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
    fplist = {}
    for fp in fptype:
        fplist[fp] = FingerprintUtils.getFingerprint(molstr, fp)
    if print_queue is not None:
        print_queue.put(("progress", {"slice": slice}))
    return pd.Series(fplist)

def standardize_row(row, molstr_type, slice = None, print_queue = None):
    molstr = row[molstr_type]
    std_molstr = FingerprintUtils.standardizeSmiles(molstr, molstr_type)
    if print_queue is not None:
        print_queue.put(("progress", {"slice": slice}))
    return pd.Series({"std_smiles": std_molstr})

def get_hashkey(data: Union[str, List[str]]):
    str = "|".join(data) if isinstance(data, list) else data
    hasher = sha256()
    hasher.update(str.encode("utf-8"))
    return hasher.hexdigest()

def calc_link_weight(loss: FragmentTreeEdge, counter):
    from_id = loss.getAttr("source")
    to_id = loss.getAttr("target")

    if from_id == to_id:
        return 1
    elif from_id == 0:
        return calc_tf_idf(loss, counter)
    else:
        return calc_PMI(loss, counter)

def calc_PMI(loss: FragmentTreeEdge, counter):
    from_frag = loss.getFromNode()
    to_frag = loss.getToNode()
    hashkey_from = from_frag.getAttr("hashkey")
    hashkey_to = to_frag.getAttr("hashkey")
    hashkey_loss = loss.getAttr("hashkey")

    from_frag_in_loss = counter[(hashkey_from, "fragment_in_loss")]
    to_frag_in_loss = counter[(hashkey_to, "fragment_in_loss")]
    loss_total = counter[(get_hashkey("loss_total"), "loss_total")]
    loss_in_sample = counter[(hashkey_loss, "loss")]
    return math.log((loss_in_sample * loss_total) / (from_frag_in_loss * to_frag_in_loss))

def calc_tf_idf(loss: FragmentTreeEdge, counter):
    frag = loss.getToNode()
    hashkey_frag = frag.getAttr("hashkey")

    tf = frag.getAttr("relativeIntensity")
    sample_total = counter[(get_hashkey("sample_total"), "sample_total")]
    frag_in_sample = counter[(hashkey_frag, "fragment")]
    idf = math.log(sample_total / frag_in_sample)
    return tf * idf

def formula2dict(formula):
    formula=formula + '='
    elems = {}
    elemstr=''
    numstr=''
    for x in range(len(formula)):
        l=formula[x]
        if l == '[':
            continue
        if l in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ]=':
            if elemstr != '':
                if numstr != '':
                    num = int(numstr)
                    numstr = ''
                else:
                    num = 1
                elems[elemstr] = num
            elemstr = l
            if l == ']':
                break
        if l in '0123456789':
            numstr += l
        if l in 'abcdefghijklmnopqrstuvwxyz':
            elemstr += l
    return elems

def calc_dataset_row(row, counter, elemset: set):
    ftree = row["ftree"]
    nodes = ftree.getAllNodes()
    edges = ftree.getAllEdges()

    nodes_features = {}

    for idx, node in nodes.items():
        formula = node.getAttr("molecularFormula")
        mz = node.getAttr("mz")
        rel_int = node.getAttr("relativeIntensity")
        elems = formula2dict(formula)
        for elem in elems.keys():
            elemset.add(elem)
        nodes_features[idx]=dict(
            formula=formula,
            mz=mz,
            rel_int=rel_int,
            elems=elems
        )

    edge_features = {}
    
    for idx, edge in edges.items():
        feature = calc_link_weight(edge, counter)
        edge_features[idx] = feature

    fingerprints = {}

    fptype = [
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

    for fp in fptype:
        fingerprints[fp] = row[fp]
    
    return pd.Series(dict(
        nodes = nodes_features,
        edges = edge_features,
        fingerprints = fingerprints
    ))
    