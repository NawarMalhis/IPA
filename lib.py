import argparse
import torch
import os
import math

from param import *
# from annotated_fasta import *
from Bayes import *
from gen_features import *
from Models import C3F2

sig = torch.nn.Sigmoid()


def get_arguments():
    st_ubc = 'Michael Smith Laboratories, The University of British Columbia (2024)'
    xx = 'IDR Probabilistic Annotation (IPA)\n'
    parser = argparse.ArgumentParser(description=f"{xx}By: Nawar Malhis, {st_ubc}.")
    parser.add_argument('-p', "--path", type=str, default='./',
                        help="Path for input & output files, default='./'")
    parser.add_argument('-in', '--in_file', type=str, help='Input file in multi sequences fasta format',
                        required=True)
    parser.add_argument('-af2', "--AlphaFold2", action='store_true',
                        help="Utilizing AlphaFold-2 RSA and pLDDT features")
    parser.add_argument('-cpu', "--cpu", action='store_true',
                        help="Force executing on a CPU")
    # parser.add_argument('-o', "--org_seq", action='store_true',
    #                     help="Do not substitute non-standard amino acids with Alanines (A)")
    arguments = parser.parse_args()
    return arguments


def load_models_info_dict():
    m_info_dict = {}
    in_file = 'stuff/models.tsv'
    with open(in_file, 'r') as fin:
        for line in fin:
            if len(line) < 5:
                continue
            lst = line.strip().split()
            if lst[0] not in m_info_dict:
                m_info_dict[lst[0]] = []
            dd = {'in_chl': int(lst[2]), 'chl_inc': int(lst[3]), 'ck_sz': int(lst[4]), 'pk_sz': int(lst[5]),
                  'fc1_out': int(lst[6])}
            m_info_dict[lst[0]].append(dd)
    return m_info_dict


def load_models(device='cpu', af2=False, verbose=False):
    m_info = load_models_info_dict()
    models = {}
    m_list = mdl_list
    if af2:
        m_list = af2_mdl_list
    for mm in m_list:
        models[mm] = {'f_used': mdl_data[mm]['FU'], 'models': []}
        for cv in [0, 1, 2, 3]:
            m_name = f'{models_path}mm_{mm}_{cv}.cnn.x'
            # print(m_name, flush=True)
            models[mm]['models'].append(C3F2(in_chl=m_info[mm][cv]['in_chl'], chl_inc=m_info[mm][cv]['chl_inc'],
                                             c1k_sz=m_info[mm][cv]['ck_sz'], c2k_sz=m_info[mm][cv]['ck_sz'],
                                             c3k_sz=m_info[mm][cv]['ck_sz'], pk_sz=m_info[mm][cv]['pk_sz'],
                                             fc1_out=m_info[mm][cv]['fc1_out']).to(device))
            models[mm]['models'][-1].load_state_dict(torch.load(m_name, map_location=device, weights_only=True))
            models[mm]['models'][-1].eval()
    if verbose:
        for mdl in models:
            print(f'================ {mdl}:')
            for mm in models[mdl]['models']:
                print(mm)
                print()
    return models


def load_priors():
    pr_dict_3 = {}
    priors_file = f"{models_path}priors.tsv"
    with open(priors_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            lst = line.split()
            tag = lst[0]
            if tag not in pr_dict_3:
                pr_dict_3[tag] = {}
            af = lst[1]
            if af not in pr_dict_3[tag]:
                pr_dict_3[tag][af] = {}
            cv = lst[2]
            pr_dict_3[tag][af][cv] = float(lst[3])
    return pr_dict_3


def assemble_data_1(ac, seq, features, features_used=None, w_out=10, pad=250):
    if features_used is None:
        features_used = ['RSA', 'pLDDt', 'w1', 'w2', 'PDB', 'IDR', 'Linker', 'P_bind', 'N_bind']
    w_in = w_out + (2 * pad)
    r_ftrs = []
    r_seq = []
    r_ac_st = []
    sz = len(seq)
    features_list = list(features.keys())
    wind_count = math.ceil(sz / w_out)

    tmp_dict = {}
    zeros_left = [0] * pad
    zeros_right = [0] * (pad + (wind_count * w_out) - sz)
    tmp_dict['seq'] = list(seq) + ['x'] * (pad + (wind_count * w_out) - sz)
    for ftr in features_list:
        tmp_dict[ftr] = list(zeros_left) + list(features[ftr]) + list(zeros_right)

    # features
    for xst in range(wind_count):
        st = xst * w_out
        r_ac_st.append((ac, st))
        r_seq.append(tmp_dict['seq'][st:st + w_out])
        tmp_lst = []
        for ftr in features_used:
            if ftr in features_list:
                tmp_lst.append(tmp_dict[ftr][st:st + w_in])
        r_ftrs.append(tmp_lst)
    return {'features': torch.tensor(np.array(r_ftrs, dtype='float32'), dtype=torch.float32),
            'seq': r_seq,
            'ac_st': r_ac_st}


def score_af2_ensemble(output_path, in_ac, in_seq, models, priors_dict, device):
    _af, features = assemble_af2_features(in_ac=in_ac, in_seq=in_seq, dbg=True)
    mdl_used = af2_mdl_list
    if not _af:
        mdl_used = []
        print(f'\tNo cif file for {in_ac} was found', flush=True)
    # else:
    #     print(flush=True)
    sz = len(in_seq)
    ym_dict = {}
    for mdl in mdl_used:
        tag = mdl[:-2]
        __af = mdl_data[mdl]['AF']
        dta = assemble_data_1(ac=in_ac, seq=in_seq, features=features, features_used=mdl_data[mdl]['FU'], w_out=10,
                              pad=250)
        scores = []
        for cv in [0, 1, 2, 3]:
            yy = sig(models[mdl]['models'][cv](dta['features'].to(device)).to('cpu').detach()).numpy().reshape(-1)[0:sz]
            yy = bayes_evidence(yy, prior=priors_dict[tag][__af][f'{cv}'])  # priors_dict[mdl]
            scores.append(np.array(list(yy), dtype='float32'))
        ym_dict[mdl] = smooth(merge_mean(scores))
    if len(ym_dict) > 0:
        save_results(output_path, in_ac, ym_dict, mdl_used, in_seq)


def smooth(s2):
    s1 = np.concatenate((np.array([s2[0]], dtype='float32'), s2))
    s0 = np.concatenate((np.array([s1[0]], dtype='float32'), s1))
    s3 = np.concatenate((s2, np.array([s2[-1]], dtype='float32')))
    s4 = np.concatenate((s3, np.array([s3[-1]], dtype='float32')))
    return (s0[:-2] + s1[:-1] + s2 + s3[1:] + s4[2:]) / 5.0


def score_ensemble(output_path, in_ac, in_seq, models, priors_dict, device):
    seq, features = assemble_features(ac=in_ac, seq=in_seq, dbg=True)
    mdl_used = mdl_list
    sz = len(seq)
    y_dict = {}
    for mdl in mdl_used:  # ['P_bind_3', 'Linker_3', 'N_bind_1']
        tag = mdl[:-2]
        __af = mdl_data[mdl]['AF']
        dta = assemble_data_1(ac=in_ac, seq=seq, features=features, features_used=mdl_data[mdl]['FU'], w_out=10,
                              pad=250)
        scores = []
        for cv in [0, 1, 2, 3]:
            yy = sig(models[mdl]['models'][cv](dta['features'].to(device)).to('cpu').detach()).numpy().reshape(-1)[0:sz]
            yy = bayes_evidence(yy, prior=priors_dict[tag][__af][f'{cv}'])  # priors_dict[mdl]
            scores.append(np.array(list(yy), dtype='float32'))
        y_dict[mdl] = smooth(merge_mean(scores))
        # y_dict[mdl] = merge_bayes(scores)
    if len(y_dict) > 0:
        save_results(output_path, in_ac, y_dict, mdl_used, seq)


def save_results(output_path, ac, y_dict, mdl_used, seq):
    for mdl in y_dict:
        m_name = mdl_data[mdl]['Name']
        out_file = f'{output_path}IPA_{m_name}/{ac}.caid'
        with open(out_file, 'w') as fout:
            print(f'>{ac}', file=fout)
            for ii, aa in enumerate(seq):
                print(f"{ii+1}\t{aa}\t{y_dict[mdl][ii]:.5}", file=fout)


def create_out_dirs(_path, af2=False):
    output_path = f"{_path}/output/"
    if not os.path.exists(output_path):
        os.system(f"mkdir {output_path}")
    for mdl in mdl_data:
        if af2:
            if mdl not in af2_mdl_list:
                continue
        else:
            if mdl not in mdl_list:
                continue
        m_name = mdl_data[mdl]['Name']
        dr = f"{output_path}IPA_{m_name}"
        if not os.path.exists(dr):
            os.system(f"mkdir {dr}")
