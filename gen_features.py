import numpy as np
import requests
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.SASA import ShrakeRupley
from crc64iso.crc64iso import crc64
import os
from param import *


def read_p_matrix_dict1(in_file):
    ret = {}
    aa_list = []
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) <= 1:
                continue
            lst = line.split()
            if lst[0] == '.':
                aa_list = lst[1:]
                for aa in aa_list:
                    ret[aa] = []
            else:
                aa = lst[0]
                # ret[aa] = np.array(list(map(float, lst[1:])), dtype='float')
                ret[aa] = np.array(lst[1:], dtype='float')
    return ret


def get_p_features(seq, _p, w1, w2, skip=0):
    w = w1 + w2 + skip
    aa_map = {'-': 20}
    features = {'w1': np.zeros(len(seq), dtype='float'), 'w2': np.zeros(len(seq), dtype='float')}
    for ii, aa in enumerate(_p):
        aa_map[aa] = ii
    _pad = '-' * w
    # print(seq, flush=True)
    p2_seq = list(_pad + seq + _pad)
    for ii, aa in enumerate(p2_seq):
        if aa not in aa_map:
            p2_seq[ii] = '-'
    aa_count = {'w1': np.zeros(21, dtype='int32'), 'w2': np.zeros(21, dtype='int32')}
    aa_count['w1'][20] = w1
    aa_count['w2'][20] = w2
    for ii in range(w1):
        aai = aa_map[p2_seq[w + ii + 1 + skip]]
        aa_count['w1'][aai] += 1
    for ii in range(w1, w1+w2):
        aai = aa_map[p2_seq[w + ii + 1 + skip]]
        aa_count['w2'][aai] += 1
    # print(aa_count['w1'].sum(), aa_count['w1'])
    # print(aa_count['w2'].sum(), aa_count['w2'])
    # print(_p[seq[0]])
    features['w1'][0] = np.dot(aa_count['w1'][:20], _p[seq[0]])
    features['w2'][0] = np.dot(aa_count['w2'][:20], _p[seq[0]])
    # print(features['w1'][0], features['w2'][0])
    for i in range(1, len(seq)):
        ii = i + w
        aai = aa_map[p2_seq[ii - w1 - 1 - skip]]
        aa_count['w1'][aai] -= 1
        aai = aa_map[p2_seq[ii + skip]]
        aa_count['w1'][aai] -= 1
        aai = aa_map[p2_seq[ii - 1 - skip]]
        aa_count['w1'][aai] += 1
        aai = aa_map[p2_seq[ii + w1 + skip]]
        aa_count['w1'][aai] += 1

        aai = aa_map[p2_seq[ii - w - skip]]
        aa_count['w2'][aai] -= 1
        aai = aa_map[p2_seq[ii + w1 + skip]]
        aa_count['w2'][aai] -= 1
        # ==============================================
        aai = aa_map[p2_seq[ii - w1 - 1 - skip]]
        aa_count['w2'][aai] += 1
        # print(len(p2_seq), ii + w, flush=True)
        aai = aa_map[p2_seq[ii + w]]
        aa_count['w2'][aai] += 1
        # print(f"{i}\t{p2_seq[i+w]}\t{seq[i]}", flush=True)
        if p2_seq[i+w] in _p:
            features['w1'][i] = np.dot(aa_count['w1'][:20], _p[p2_seq[i+w]])
            features['w2'][i] = np.dot(aa_count['w2'][:20], _p[p2_seq[i+w]])
        else:
            features['w1'][i] = 0
            features['w2'][i] = 0
    return features


def gen_features_p(ac, seq, _w1=15, _w2=20, _skip=0, dbg=False):
    _p = read_p_matrix_dict1(f"stuff/P.tsv")
    features = get_p_features(seq=seq, _p=_p, w1=_w1, w2=_w2, skip=_skip)
    if dbg:
        features_file = f"features/P.f"
        with open(features_file, 'w') as fout:
            print(f"# AC ...:\t{ac}\n# length:\t{len(seq):,}\n# Seq ..:\t{seq}\n#", file=fout)
            print(f"# Feature\tw1\t{_w1}", file=fout)
            print(f"# Feature\tw2\t{_w2}", file=fout)
            print("#", file=fout)
            for ii in range(len(features['w1'])):
                print(f"{ii + 1}\t{features['w1'][ii]:.4}\t{features['w2'][ii]:.4}", file=fout)
    return features


def read_counts(in_file):
    ret = {}
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            lst = line.split()
            if lst[0] not in ret:
                ret[lst[0]] = {}
            ret[lst[0]][lst[1]] = float(lst[2]) / 2.0
    return ret


def gen_features_counts(ac, seq, dbg=False):
    features_counts = {'PDB': [], 'IDR': [], 'Linker': [], 'P_bind': [], 'N_bind': []}
    _pp_features = read_counts(in_file="stuff/LNPPI.tsv")
    if dbg:
        features_file = "features/aac.f"
        with open(features_file, 'w') as fout:
            print(f"# AC ...:\t{ac}\n# length:\t{len(seq):,}\n# Seq ..:\t{seq}\n#", file=fout)
            for tg in ['PDB', 'IDR', 'Linker', 'P_bind', 'N_bind']:
                print(f"# Feature\tf{tg}", file=fout)
            print("#", file=fout)
            for ii, aa in enumerate(seq):
                print(f"{ii + 1}", end='', file=fout)
                for tg in ['PDB', 'IDR', 'Linker', 'P_bind', 'N_bind']:
                    if aa in _pp_features[tg]:
                        print(f"\t{_pp_features[tg][aa]:.5}", end='', file=fout)
                        features_counts[tg].append(_pp_features[tg][aa])
                    else:
                        print(f"\t50.0", end='', file=fout)
                        features_counts[tg].append(50.0)
                print(file=fout)
    else:
        for ii, aa in enumerate(seq):
            for tg in ['PDB', 'IDR', 'Linker', 'P_bind', 'N_bind']:
                if aa in _pp_features[tg]:
                    features_counts[tg].append(_pp_features[tg][aa])
                else:
                    features_counts[tg].append(50.0)
    return features_counts


def get_url(url, **kwargs):
    response = requests.get(url, **kwargs);
    if not response.ok:
        print(response.text)
        return None
        # response.raise_for_status()
        # sys.exit()
    return response


def get_ac_list(seq):
    ac_list = []
    checksum = crc64(seq)
    r = get_url(f"https://rest.uniprot.org/uniparc/search?query=checksum: {checksum}")
    if r is not None:
        data = r.json()
        # for xx in data["results"][0]['uniParcCrossReferences']:
            # if 'UniProtKB' in xx['database']:
            #     ac_list.append(xx['id'])
        for xx in data["results"][0]['uniProtKBAccessions']:
            ac_list.append(xx)
    return ac_list


def cif_file_name(ac):
    return f"{ac}.cif"


def get_af2(ac, part=1):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{ac}-F{part}-model_v4.cif"
    try:
        print(f"\t\tExtracting AlphaFold-2 structure for {ac}:\t", end='')
        response = requests.get(url)
        if response.status_code == 200:
            print('Found', flush=True)
            with open(f"{af2_cif_path}{cif_file_name(ac)}", "wb") as file:
                file.write(response.content)
                return True
        else:
            print('Not found', flush=True)
            return False  # File does not exist or other error
    except requests.ConnectionError:
        print(f"Connection error: {ac}")
        return False  # Connection error, file may not exist or server is unreachable


def gen_features_alpha_fold(in_ac, in_seq, dbg=False):
    max_sasa = read_max_sasa("stuff/AminoAcids.tsv")
    cif_file = f"{af2_cif_path}{in_ac}.cif"  # cif_file_name(ac)
    if os.path.exists(cif_file):
        seq, dta = get_cif_data(cif_file, max_sasa)
        if seq != in_seq:
            print(f"(Inconsistent {in_ac} cif file)", end='\t')
            return None
        features_cif = {'pLDDt': [], 'RSA': []}
        if dbg:
            features_file = f"features/AlphaFold.f"
            with open(features_file, 'w') as fout:
                print(f"# AC ...:\t{in_ac}\n# length:\t{len(seq):,}\n# Seq ..:\t{seq}\n#", file=fout)
                print(f"# Feature\tpLDDT", file=fout)
                print(f"# Feature\tRSA", file=fout)
                print("#", file=fout)
                for ii, ln_dict in enumerate(dta):
                    print(f"{ii+1}\t{ln_dict['pLDDT']}\t{ln_dict['RSA']:.5}", file=fout)
                    features_cif['pLDDt'].append(ln_dict['pLDDT'])
                    features_cif['RSA'].append(ln_dict['RSA'] * 20)
        else:
            for ii, ln_dict in enumerate(dta):
                features_cif['pLDDt'].append(ln_dict['pLDDT'])
                features_cif['RSA'].append(ln_dict['RSA'] * 20)
        return features_cif
    return None


def assemble_af2_features(in_ac, in_seq, dbg=True):
    features = {}
    ftr_lst = []
    _af = False
    ac_list = get_ac_list(seq=in_seq)
    # print(ac_list)
    tmp_list = list(ac_list)
    for ac in tmp_list:
        if '.' in ac:
            ac_list.remove(ac)
            continue
        if '-' in ac:
            ac0 = ac.split('-')[0]
            if ac0 not in ac_list:
                ac_list.append(ac0)
    for ac in ac_list:
        if '.' in ac:
            continue
        if not os.path.isfile(f"{af2_cif_path}{ac}.cif"):
            if download_cif:
                get_af2(ac, part=1)
        features_cif = gen_features_alpha_fold(in_ac=ac, in_seq=in_seq)
        if features_cif is not None:
            ftr_lst.append(features_cif)
            _af = True
            break
    if not _af:
        return False, None
    features_p = gen_features_p(in_ac, in_seq, _w1=20, _w2=25, _skip=1)
    ftr_lst.append(features_p)
    features_counts = gen_features_counts(in_ac, in_seq)
    ftr_lst.append(features_counts)
    for f_dict in ftr_lst:
        for ff in f_dict:
            features[ff] = f_dict[ff]
    for ff in features:
        if ff not in ['w1', 'w2']:
            features[ff] = np.array(features[ff], dtype='float32')
    return _af, features


def assemble_features(ac, seq, dbg=True):
    features = {}
    ftr_lst = []

    features_p = gen_features_p(ac, seq, _w1=20, _w2=25, _skip=1)
    ftr_lst.append(features_p)
    features_counts = gen_features_counts(ac, seq)
    ftr_lst.append(features_counts)
    for f_dict in ftr_lst:
        for ff in f_dict:
            features[ff] = f_dict[ff]
    for ff in features:
        if ff not in ['w1', 'w2']:
            features[ff] = np.array(features[ff], dtype='float32')
    return seq, features


def read_fasta(fasta_file):
    _seq = ''
    _ac = ''
    with open(fasta_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                _ac = line[1:].split('|')[0]
                continue
            if _ac != '':
                _seq = _seq + line
    return _ac, _seq


def read_max_sasa(in_file):
    # print('===================', flush=True)
    ret = {}
    with open(in_file, 'r') as fin:
        for line in fin:
            lst = line.strip().split()
            ret[lst[0]] = {'AA': lst[1], 'max': int(lst[2]), 'name': lst[3]}
    return ret


def get_cif_data(cif_file, max_sasa):
    dta = []
    lst = []
    parser = MMCIFParser()
    structure = parser.get_structure("xx", cif_file)
    mmcif_dict = MMCIF2Dict(cif_file)
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    for pos in range(1, len(structure[0]['A'])+1):
        aaa = structure[0]["A"][pos].get_resname()
        aa = max_sasa[aaa]['AA']
        sasa = structure[0]['A'][pos].sasa
        rsa = sasa/max_sasa[aaa]['max']
        plddt = float(mmcif_dict['_ma_qa_metric_local.metric_value'][pos-1])
        dta.append({'pLDDT': plddt, 'RSA': rsa, 'AA': aa})
        lst.append(aa)
    return ''.join(lst), dta
