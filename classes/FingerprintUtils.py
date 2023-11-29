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

from typing import Optional, Union, Tuple, Dict

import pandas as pd
from bitarray import bitarray
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Avalon import pyAvalonTools
from .CDKImporter import *

class FingerprintUtils:
    instance_container = {}
    instance_counter = {}

    @classmethod
    def getCDKInstance(cls, name):
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
    def getFingerprint(cls, molstr, fptype) -> bitarray:
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
    def isSmiles(cls, molstr) -> Tuple[bool, Dict[str, bool]]:
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
    def isINCHI(cls, molstr) -> Tuple[bool, Dict[str, bool]]:
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
    def standardizeSmiles(cls, molstr, molstr_type = "smiles") -> Optional[str]:
        mol = pybel.readstring(molstr_type, molstr)
        new_str = mol.write("can").strip()
        return None if new_str == "" else new_str

