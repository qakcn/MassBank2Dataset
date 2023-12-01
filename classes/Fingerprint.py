# This script contains functions to validate SMILES and INCHI strings and 
# generate fingerprint form those strings.
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
from typing import Optional, Tuple, Dict
import math

# Third-party imports
import pandas as pd
from bitarray import bitarray
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Avalon import pyAvalonTools

# Local imports
from .CDKImporter import *
from .Functions import get_hashkey
from .FragmentTree import FragmentTreeEdge

class FingerprintUtils:
    instance_container = {}
    instance_counter = {}

    @classmethod
    def getCDKInstance(cls, name: str):
        cdk = import_cdk()
        cdklist = {
            "CDKFingerprint": lambda: cdk.fingerprint.Fingerprinter(),
            "PubChemFingerprint": lambda: cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance()),
            "Klekota-Roth": lambda: cdk.fingerprint.KlekotaRothFingerprinter(),
            "SmilesParser": lambda: cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
        }

        if name in cls.instance_container and cls.instance_counter[name] < 500:
            instance = cls.instance_container[name]
            cls.instance_counter[name] += 1
        else:
            instance = cdklist[name]()
            cls.instance_container[name] = instance
            cls.instance_counter[name] = 1
        return instance

    @classmethod
    def getFingerprint(cls, molstr: str, fptype: str) -> bitarray:
        rdkitlist = {
            "smiles": Chem.MolFromSmiles,
            "inchi": Chem.MolFromInchi,
            "AtomPair": AllChem.GetHashedAtomPairFingerprintAsBitVect,
            "Avalon": pyAvalonTools.GetAvalonFP,
            "MACCS": AllChem.GetMACCSKeysFingerprint,
            "Morgan": lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2),
            "TopologicalTorsion": AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect,
            "RDKitFingerprint": Chem.RDKFingerprint
        }
        cdklist = {
            "CDKFingerprint": 1024,
            "PubChemFingerprint": 881,
            "Klekota-Roth": 4860
        }

        if molstr.startswith("InChI="):
            molstr_type = "inchi"
        else:
            molstr_type = "smiles"
        
        hash2bit = lambda x, n: [1 if i in x else 0 for i in range(n)]
        if fptype == "FP2":
            mol = pybel.readstring(molstr_type, molstr)
            fp = mol.calcfp(fptype)
            fplist = hash2bit(fp.bits, 1024)
        elif fptype in rdkitlist:
            mol = rdkitlist[molstr_type](molstr)
            fp = rdkitlist[fptype](mol)
            fplist = fp.ToList()
            if fptype == "MACCS":
                fplist = fplist[1:]
        elif fptype in cdklist:
            if molstr_type == "inchi":
                cdk = import_cdk()
                mol = cdk.inchi.InChIToStructure(molstr, cdk.DefaultChemObjectBuilder.getInstance()).getAtomContainer()
            elif molstr_type == "smiles":
                mol = FingerprintUtils.getCDKInstance("SmilesParser").parseSmiles(molstr)
            fp = FingerprintUtils.getCDKInstance(fptype).getBitFingerprint(mol)
            fplist = hash2bit(list(fp.getSetbits()), cdklist[fptype])
        return bitarray(fplist)

    @classmethod
    def isSmiles(cls, molstr: str) -> Tuple[bool, Dict[str, bool]]:
        smiles_pybel_fail = False
        smiles_rdkit_fail = False
        smiles_cdk_fail = False
        smiles_rdkit2cdk_fail = False
        smiles_rdkit2pybel_fail = False

        try:
            mol = pybel.readstring("smiles", molstr)
        except:
            smiles_pybel_fail = True
        try:
            mol = Chem.MolFromSmiles(molstr)
            if mol is None:
                raise Exception
        except:
            smiles_rdkit_fail = True
            smiles_rdkit2cdk_fail = True
            smiles_rdkit2pybel_fail = True
        try:
            mol = cls.getCDKInstance("SmilesParser").parseSmiles(molstr)
        except:
            smiles_cdk_fail = True
        if not smiles_rdkit_fail and smiles_cdk_fail:
            try:
                mol = Chem.MolFromSmiles(molstr)
                molstr2 = Chem.MolToSmiles(mol)
                mol = cls.getCDKInstance("SmilesParser").parseSmiles(molstr2)
            except:
                smiles_rdkit2cdk_fail = True
        if not smiles_rdkit_fail and smiles_pybel_fail:
            try:
                mol = Chem.MolFromSmiles(molstr)
                molstr2 = Chem.MolToSmiles(mol)
                mol = pybel.readstring("smiles", molstr2)
            except:
                smiles_rdkit2pybel_fail = True
        
        status = dict(pybel_fail=smiles_pybel_fail, rdkit_fail=smiles_rdkit_fail, cdk_fail=smiles_cdk_fail, rdkit2cdk_fail=smiles_rdkit2cdk_fail, rdkit2pybel_fail=smiles_rdkit2pybel_fail)
        return not any(status.values()), status

    @classmethod
    def isINCHI(cls, molstr: str) -> Tuple[bool, Dict[str, bool]]:
        inchi_pybel_fail = False
        inchi_rdkit_fail = False
        inchi_cdk_fail = False
        try:
            mol = pybel.readstring("inchi", molstr)
        except:
            inchi_pybel_fail = True
        try:
            mol = Chem.MolFromInchi(molstr)
            if mol is None:
                raise Exception
        except:
            inchi_rdkit_fail = True
        try:
            cdk=import_cdk()
            mol = cdk.inchi.InChIToStructure(molstr, cdk.DefaultChemObjectBuilder.getInstance()).getAtomContainer()
        except:
            inchi_cdk_fail = True
        status = dict(pybel_fail=inchi_pybel_fail, rdkit_fail=inchi_rdkit_fail, cdk_fail=inchi_cdk_fail)
        
        return not any(status.values()), status

    @classmethod
    def standardizeSmiles(cls, molstr: str, molstr_type: str = "smiles") -> Optional[str]:
        mol = pybel.readstring(molstr_type, molstr)
        new_str = mol.write("can").strip()
        return None if new_str == "" else new_str


def test_molstr_row(row) -> pd.Series:
    smiles = row["smiles"]
    inchi = row["inchi"]

    succ_smiles, res_smiles = FingerprintUtils.isSmiles(smiles)
    succ_inchi, res_inchi = FingerprintUtils.isINCHI(inchi)

    return pd.Series({
        "test_smile": succ_smiles,
        "test_inchi": succ_inchi,
        "test_results": {"smiles": res_smiles, "inchi": res_inchi}
    })

def fingerprint_row(row, slice: int = None, print_queue = None) -> pd.Series:
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

def standardize_row(row, molstr_type: str) -> pd.Series:
    molstr = row[molstr_type]
    std_molstr = FingerprintUtils.standardizeSmiles(molstr, molstr_type)
    return pd.Series({"std_smiles": std_molstr})

def calc_link_weight(loss: FragmentTreeEdge, counter: Dict[Tuple[str, str], int]) -> float:
    from_id = loss.getAttr("source")
    to_id = loss.getAttr("target")

    if from_id == to_id:
        return 1
    elif from_id == 0:
        return calc_tf_idf(loss, counter)
    else:
        return calc_PMI(loss, counter)

def calc_PMI(loss: FragmentTreeEdge, counter: Dict[Tuple[str, str], int]) -> float:
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

def calc_tf_idf(loss: FragmentTreeEdge, counter: Dict[Tuple[str, str], int]) -> float:
    frag = loss.getToNode()
    hashkey_frag = frag.getAttr("hashkey")

    tf = frag.getAttr("relativeIntensity")
    sample_total = counter[(get_hashkey("sample_total"), "sample_total")]
    frag_in_sample = counter[(hashkey_frag, "fragment")]
    idf = math.log(sample_total / frag_in_sample)
    return tf * idf

def formula2dict(formula: str) -> Dict[str, int]:
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

def calc_dataset_row(row, counter: Dict[Tuple[str, str], int], elemset: Dict[str, int]) -> pd.Series:
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