# This is the 4th step of generating data set from MassBank database.
# This script parses the fragment tree from the json data.
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
import pickle
from threading import Thread
from queue import Queue

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from classes import *
from classes.Functions import get_hashkey
from classes.FragmentTree import parse_ftree_row

##################################################
# Parameters that can be edited by the user      #
##################################################
# Paths
output_path = Path("outputs")
intermediate_path = Path("intermediates")

# Files
spectra_ftrees_file = intermediate_path / "spectra.ftrees.pkl"

vocab_file = output_path / "vocab.nopd.pkl"
orphan_list_file = output_path / "orphan.nopd.pkl"

counter_file = intermediate_path / "counter.nopd.pkl"

spectra_ftrees_parsed_file = intermediate_path / "spectra.ftrees.parsed.pkl"
spectra_ftrees_unparsed_file = intermediate_path / "spectra.ftrees.unparsed.pkl"

##################################################
# End of parameters, do not edit below this line #
##################################################

def print_orphan(print_queue):
    cnt = 0
    while True:
        spectrum_id, compound_id, compound_name = print_queue.get()
        if spectrum_id is None:
            break
        cnt += 1
        TS.p(f"Orphan {TS.blue(cnt)}: {TS.cyan(spectrum_id)}::{TS.yellow(compound_id)}::{TS.magenta(compound_name)}")

# Initialization
TS.register_printer(tqdm.write)

if not spectra_ftrees_file.is_file():
    TS.p(TS.red("Error loading ftrees spectra from pickle file...file not exists."))
else:
    TS.ip("Loading ftrees spectra from pickle file...")
    spectra = pd.read_pickle(spectra_ftrees_file)
    TS.p(TS.green("Done."))

    total = spectra.shape[0]

    vocab = {
        "orphan": {},
        "sample": {},
        "fragment": {},
        "loss": {},
    }
    counter={}
    orphan_list=[]

    counter[(get_hashkey("sample_total"), "sample_total")] = total

    print_queue = Queue()

    pot = Thread(target=print_orphan, args=(print_queue,))
    pot.start()

    tqdm.pandas(desc=f"Parsing", total=total, ncols=100)
    spectra["ftree"] = spectra.progress_apply(parse_ftree_row, axis=1, vocab=vocab, counter=counter, orphan_list=orphan_list, print_queue=print_queue)
    TS.p(TS.green(f"Parsed."))

    print_queue.put((None,None,None))
    pot.join()

    TS.ip("Saving spectra...")
    mask = spectra["ftree"].apply(lambda x: x.__class__.__name__ == "FragmentTree")
    spectra_parsed = spectra[mask].copy()
    spectra_unparsed = spectra[~mask].copy()

    spectra_parsed.to_pickle(spectra_ftrees_parsed_file)
    spectra_unparsed.to_pickle(spectra_ftrees_unparsed_file)
    TS.p(TS.green("Done."))

    # After multi-threads parsing
    TS.ip("Saving vocabulary...")
    with open(vocab_file, "wb") as f:
        pickle.dump(vocab, f)
    TS.p(TS.green("Done."))

    TS.ip("Saving counter...")
    with open(counter_file, "wb") as f:
        pickle.dump(counter, f)
    TS.p(TS.green("Done."))

    TS.ip("Saving orphan list...")
    with open(orphan_list_file, "wb") as f:
        pickle.dump(orphan_list, f)
    TS.p(TS.green("Done."))
