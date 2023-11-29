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
import math
from typing import Optional, Dict

# Local imports
from .VocabUtils import *

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
    
    def getNode(self, idx):
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
