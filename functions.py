import MySQLdb
import pickle
from pathlib import Path
import json
from bitarray import bitarray
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import rdmolops, AllChem
from rdkit.Avalon import pyAvalonTools
from classes import *


def get_compounds(conn: MySQLdb.Connection, path: Path, level='MS2', msfile=False):
    qs = r"SELECT COUNT(*) FROM `msms_spectrum` WHERE `ms_level`='MS2'"
    c = conn.cursor()
    c.execute(qs)
    r = c.fetchone()
    c.close()
    num = r[0]

    cycle = num // 1000

    compoundlist = []
    cnt = 0
    for i in range(cycle + 1):
        qst = r"SELECT a.`spectrum_id`,a.`spectrum_name`,a.`compound_id`,a.`precursor_mz_text`,a.`adduct`," \
              r"a.`polarity`,a.`collision_energy_text`,b.`formula`,b.`exactmass`,b.`smiles`,b.`inchi` FROM " \
              r"`msms_spectrum` AS a JOIN `ms_compound` AS b ON a.`compound_id`=b.`compound_id` WHERE " \
              r"a.`ms_level`='MS2' LIMIT {page:d},1000"
        qs = qst.format(page=i * 1000)
        column_name = ('spectrum_id', 'spectrum_name', 'compound_id', 'precursor_mz_text', 'adduct', 'polarity',
                       'collision_energy_text', 'formula', 'exactmass', 'smiles', 'inchi')
        c = conn.cursor()
        c.execute(qs)
        rs = c.fetchall()
        c.close()
        for r in rs:
            compound = {}
            cnt += 1
            print("\rProcessing Compound {cnt:>5d}/{num:>5d}.".format(cnt=cnt, num=num), flush=True, end='')
            for j in range(len(r)):
                compound[column_name[j]] = r[j] if r[j] is not None else ''

            compound['spectrum_peak'] = get_spectrum_peak(compound['spectrum_id'], conn)
            compound['compound_name'] = get_compound_name(compound['compound_id'], conn)
            compoundlist.append(compound)
            if msfile:
                gen_ms_file(compound, path / level)
    print('')
    (path / ('all.' + level + '.pkl')).write_bytes(pickle.dumps(compoundlist, 4))
    print("All Finished!")
    return compoundlist


def get_spectrum_peak(sid, conn: MySQLdb.Connection):
    qst = r"SELECT `mz`,`intensity` FROM `msms_spectrum_peak` WHERE `spectrum_id`='{id:s}'"
    qs = qst.format(id=sid)
    c = conn.cursor()
    c.execute(qs)
    rs = c.fetchall()
    c.close()
    return rs


def get_compound_name(cid, conn: MySQLdb.Connection):
    qst = r"SELECT `NAME`.`CH_NAME` FROM `COMPOUND_NAME` JOIN `NAME` ON `COMPOUND_NAME`.`NAME`=`NAME`.`ID` WHERE " \
          r"`COMPOUND_NAME`.`COMPOUND`={id:d}"
    qs = qst.format(id=cid)
    c = conn.cursor()
    c.execute(qs)
    rs = c.fetchall()
    c.close()
    return [n for n, *_ in rs]


def gen_ms_file(compound, path: Path):
    path.mkdir(parents=True, exist_ok=True)
    file = path / (str(compound['compound_id']) + '.ms')
    tpl = """>compound {name}
>formula {form}
>ionization {ion}
>charge {charge}
>parentmass {premass}
    
>ms2
{peaks}"""
    if compound['polarity'] == 'POSITIVE':
        charge = '1'
    elif compound['polarity'] == 'NEGATIVE':
        charge = '-1'
    else:
        charge = ''
    peaks = ''
    for mz, its in compound['spectrum_peak']:
        peaks += str(mz) + ' ' + str(its) + '\n'

    file.write_text(tpl.format(
        name=compound['compound_name'][0],
        form=compound['formula'],
        ion=compound['adduct'].replace('[', '\\[').replace(']', '\\]'),
        charge=charge,
        premass=compound['precursor_mz_text'],
        peaks=peaks
    ))


def increase(i, ilist: dict):
    if i in ilist:
        ilist[i] += 1
    else:
        ilist[i] = 1


def parse_ftrees(ftree_path: Path, compounds_file: Path, output_path: Path):
    compounds = pickle.loads(compounds_file.read_bytes())
    cids = {}
    for i in range(len(compounds)):
        cids[compounds[i]['compound_id']] = i
    output_path.mkdir(parents=True, exist_ok=True)
    cnt = 0
    orplist = []
    errlist = []
    allnum = sum([f.name[-5:] == '.json' for f in ftree_path.iterdir()])
    print()

    if (output_path / 'vocabs.pkl').exists():
        Vocab.load(pickle.loads((output_path / 'vocabs.pkl').read_bytes()))

    for f in ftree_path.iterdir():
        if f.name[-5:] == '.json':
            cnt += 1
            sn, cid, *_ = f.name.split('_')

            compound = compounds[cids[int(cid)]]

            Vocab.get('root').add(cid, compound)
            rootid = Vocab.get('root').get_index(cid)

            print('\rProcessing [{cnt:>5d}/{all:>5d}] {sn:>5s}_{id:>5s}_{name:30s}'.format(
                cnt=cnt, all=allnum, sn=sn, id=cid, name=compound['compound_name'][0]
            ), end='', flush=True)

            try:
                graph = json.loads(f.read_bytes())
                root = graph['root']

                if len(graph['fragments']) == 1:
                    orplist.append((sn, cid, compound['compound_name'][0], root))
                    print('Orphan', len(orplist))
                    continue

                tree = Ftree.gen_ftree(rootid)

                for frag in graph['fragments']:
                    formula = frag['molecularFormula']
                    if formula == root:
                        tree.add_frag(frag['id'], rootid, True, **frag)
                    else:
                        Vocab.get('frag').add(formula)
                        fragid = Vocab.get('frag').get_index(formula)
                        tree.add_frag(frag['id'], fragid, **frag)

                for loss in graph['losses']:
                    fromidx = loss['source']
                    toidx = loss['target']
                    fromfragid = tree.get_fragid(fromidx)
                    tofragid = tree.get_fragid(toidx)
                    Vocab.get('link').add((fromfragid, tofragid))
                    linkid = Vocab.get('link').getIndex((fromfragid, tofragid))
                    tree.add_link(fromidx, toidx, fromfragid, tofragid, linkid, **loss)

            except Exception:
                errlist.append((sn, cid, compound['compound_name'][0]))
                print('Error', len(errlist))
                Ftree.del_ftree(rootid)
                continue

        

    (output_path / 'ftrees.pkl').write_bytes(pickle.dumps(Ftree.save(), 4))
    (output_path / 'vocabs.pkl').write_bytes(pickle.dumps(Vocab.save(), 4))
    (output_path / 'orphan.ftrees.pkl').write_bytes(pickle.dumps(orplist, 4))
    (output_path / 'error.ftrees.pkl').write_bytes(pickle.dumps(errlist, 4))

    print('\nFinished.')


def gen_dataset(outputpath: Path, doload=False):
    if doload:
        Ftree.load(pickle.loads((outputpath / 'ftrees.pkl').read_bytes()))
        Vocab.load(pickle.loads((outputpath / 'vocabs.pkl').read_bytes()))

    ftree2fp_dataset = []
    mol2fp_dataset = []
    trees = Ftree.get_all_ftrees()
    allnum = len(trees)
    cnt = 0

    nosmileslist = []

    print()

    for tree in trees:
        cnt += 1
        rootid = tree.rootid
        compound = Vocab.get('root').getDataByIndex(rootid)

        node, node_attr, edge_index, edge_attr = tree.get_data()

        smiles = compound['smiles']

        if smiles == 'N/A':
            nosmileslist.append(dict(
                rootid=rootid,
                compound_id=compound['compound_id'],
                compound_name=compound['compound_name'][0]
            ))
            continue

        FP_LIST ={
            "FP2" : 1024,
            "AtomPair" : 2048,
            "Avalon" : 512,
            "MACCS" : 166,
            "Morgan" : 2048,
            "RDKitFingerprint" : 2048,
            "TopologicalTorsion" : 2048,
            "CDKFingerprint": 1024,
            "PubChemFingerprint": 881,
            "Klekota-Roth": 4860,
        }

        fingerprint={}

        for fptype in FP_LIST.keys():
            fingerprint[fptype]=get_fingerprint(smiles, fptype)

        ftree2fp_dataset.append(dict(
            node=node,
            node_attr=node_attr,
            edge_index=edge_index,
            edge_attr=edge_attr,
            fingerprint=fingerprint
        ))

        element, bond_index, bond_attr = get_mol_data(smiles)

        mol2fp_dataset.append(dict(
            element=element,
            bond_index=bond_index,
            bond_attr=bond_attr,
            fingerprint=fingerprint
        ))

        print('\rProcessing [{cnt:>5d}/{all:>5d}] {name:55s}'.format(
            cnt=cnt, all=allnum, name=compound['compound_name'][0]
        ), end='', flush=True)

    (outputpath / 'ftree2fp.dataset.pkl').write_bytes(pickle.dumps(ftree2fp_dataset, 4))
    (outputpath / 'mol2fp.dataset.pkl').write_bytes(pickle.dumps(mol2fp_dataset, 4))
    (outputpath / 'vocabs_2.pkl').write_bytes(pickle.dumps(Vocab.save(), 4))
    (outputpath / 'nosmiles.pkl').write_bytes(pickle.dumps(nosmileslist, 4))

    print('\nFinished.')


def get_mol_data(smiles):
    mol = Chem.MolFromSmiles(smiles)
    adj = rdmolops.GetAdjacencyMatrix(mol, useBO=True)
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    from_idx=[]
    to_idx=[]
    weight=[]
    elem_idx=[]
    for i in range(len(adj)):
        for j in range(len(adj)):
            if adj[i][j] != 0:
                from_idx.append(i)
                to_idx.append(j)
                weight.append(adj[i][j])
    for s in symbols:
        Vocab.get('elem').add(s)
        eid = Vocab.get('elem').getIndex(s)
        elem_idx.append(eid)
    return elem_idx, [from_idx, to_idx], weight
    

def get_fingerprint(molstr, fptype):
    rdkitlist = {
        'smiles': Chem.MolFromSmiles,
        'inchi': Chem.MolFromInchi,
        'AtomPair': AllChem.GetHashedAtomPairFingerprintAsBitVect,
        'Avalon': pyAvalonTools.GetAvalonFP,
        'MACCS': AllChem.GetMACCSKeysFingerprint,
        'Morgan': lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2),
        'TopologicalTorsion': AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect,
        'RDKitFingerprint': Chem.RDKFingerprint
    }
    cdklist = {
        'CDKFingerprint': lambda: cdk.fingerprint.Fingerprinter(),
        'PubChemFingerprint': lambda: cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance()),
        'Klekota-Roth': lambda: cdk.fingerprint.KlekotaRothFingerprinter(),
    }

    if molstr.startswith('InChI='):
        molstr_type = 'inchi'
    else:
        molstr_type = 'smiles'
    
    hash2bit = lambda x, n: [1 if i in x else 0 for i in range(n)]
    if fptype == 'FP2':
        mol = pybel.readstring(molstr_type, molstr)
        fp = mol.calcfp(fptype)
        fplist = hash2bit(fp.bits, 1024)
    elif fptype in rdkitlist:
        mol = rdkitlist[molstr_type](molstr)
        fp = rdkitlist[fptype](mol)
        fplist = fp.ToList()
        if fptype == 'MACCS':
            fplist = fplist[1:]
    elif fptype in cdklist:
        cdklen = {"CDKFingerprint": 1024,
              "PubChemFingerprint": 881,
              "Klekota-Roth": 4860}
        if molstr_type == 'inchi':
            mol = cdk.inchi.InChIToStructure(molstr, cdk.DefaultChemObjectBuilder.getInstance()).getAtomContainer()
        elif molstr_type == 'smiles':
            parser = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
            mol = parser.parseSmiles(molstr)
        fp = cdklist[fptype]().getBitFingerprint(mol)
        fplist = hash2bit(list(fp.getSetbits()), cdklen[fptype])
    return bitarray(fplist)
