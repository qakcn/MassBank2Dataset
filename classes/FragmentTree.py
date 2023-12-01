# This script contains classes that representing fragment tree and its node and 
# edge.
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
from typing import Optional, Dict
from pathlib import Path
import json

# Local imports
from .Functions import get_hashkey, inst_profile

# Classes
class FragmentTree:
    def __init__(self, spectrum_id, attrs, compound: dict):
        self.spectrum_id = spectrum_id
        self.nodes = {}
        self.edges = {}
        self.compound = compound
        self.attrs = attrs

        self.addNode(0)

    def addNode(self, idx: int, attrs: Optional[Dict] = None):
        self.nodes[idx] = FragmentTreeNode(self, attrs)
        return self
    
    def getRootNode(self):
        return self.getNode(0)
    
    def getNode(self, idx: int):
        return self.nodes[idx]

    def addEdge(self, from_idx: int, to_idx: int, attrs: dict = {}):
        fromNode = self.getNode(from_idx)
        toNode = self.getNode(to_idx)
        self.edges[(from_idx, to_idx)] = FragmentTreeEdge(self, fromNode, toNode, attrs)
        return self
    
    def getEdge(self, from_idx: int, to_idx: int):
        if (from_idx, to_idx) in self.edges:
            return self.edges[(from_idx, to_idx)]
    
    def getAllNodes(self):
        return self.nodes

    def getAllEdges(self):
        return self.edges

class FragmentTreeNode:
    def __init__(self, rootFT : "FragmentTree", attrs: dict = {}):
        self.rootFT = rootFT
        self.edgefromthis = []
        self.edgetothis = []
        self.attrs = attrs if attrs is not None else {}

    def setAttrs(self, attrs: dict):
        self.attrs.update(attrs)
        return self
    
    def setAttr(self, key, value):
        self.attrs[key] = value
        return self

    def getAttrs(self):
        return self.attrs
    
    def getAttr(self, key):
        if key in self.attrs:
            return self.attrs[key]
        else:
            raise ValueError(f"key {key} not exists in attrs")
    
    def getRootFT(self):
        return self.rootFT
    
    def addEdgeFromThisNode(self, edge: "FragmentTreeEdge"):
        self.edgefromthis.append(edge)
        return self
    
    def addEdgeToThisNode(self, edge: "FragmentTreeEdge"):
        self.edgetothis.append(edge)
        return self
    
    def getEdgeFromThisNode(self):
        return self.edgefromthis
    
    def detEdgesToThisNode(self):
        return self.edgetothis


class FragmentTreeEdge:
    def __init__(self, rootFT : "FragmentTree", fromNode: "FragmentTreeNode", toNode: "FragmentTreeNode", attrs: dict = {}):
        self.rootFT = rootFT
        self.fromNode = fromNode
        self.toNode = toNode
        self.attrs = attrs

        fromNode.addEdgeFromThisNode(self)
        toNode.addEdgeToThisNode(self)

    def setAttrs(self, attrs: dict):
        self.attrs.update(attrs)
        return self

    def setAttr(self, key, value):
        self.attrs[key] = value
        return self

    def getAttrs(self):
        return self.attrs
    
    def getAttr(self, key):
        if key in self.attrs:
            return self.attrs[key]
        else:
            raise ValueError(f"key {key} not exists in attrs")
    
    def getFromNode(self):
        return self.fromNode
    
    def getToNode(self):
        return self.toNode
    
    def getRootFT(self):
        return self.rootFT

# Functions
def process_ftree_file(row, ftrees_path: Path, tqdm_bars: dict):
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

def parse_ftree_row(row, vocab: dict, counter: dict, orphan_list: list, print_queue = None) -> FragmentTree:
    spectrum_id = row["spectrum_id"]
    compound_id = row["compound_id"]
    compound_name = row["compound_name"][0]

    ftree_json = row["ftree"]

    if len(ftree_json["fragments"]) == 1:
        orphan_list.append(row)
        print_queue.put((spectrum_id, compound_id, compound_name))

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
