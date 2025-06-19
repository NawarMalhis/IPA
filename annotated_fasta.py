import copy
import os


def annotated_fasta():
    return {'data': {}, 'metadata': {'tags': [], 'name_tags': [], 'statistics': None}}


def annotated_fasta_load(in_file: str, _mark=None):
    af_sequences = {}
    tags = []
    name_tags = []
    _more_tags = True
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                if _more_tags:
                    _lst = line.split()
                    if len(_lst) > 1:
                        if _lst[1] == 'TAG':
                            tags.append(_lst[2])
                continue
            _more_tags = False
            if len(tags) == 0:
                tags.append('mask')
            if line[0] == '>':
                ac_lst = line[1:].split('|')
                ac = ac_lst[0]
                af_sequences[ac] = {'seq': ''}
                for extra in ac_lst[1:]:
                    ex_lst = extra.split('=')
                    af_sequences[ac][ex_lst[0]] = ex_lst[1]
                    if ex_lst[0] not in name_tags:
                        name_tags.append(ex_lst[0])
                #     if len(ac_lst) == 2:
                #         ox_lst = ac_lst[1].split('=')
                #         if len(ox_lst) > 1:
                #             if ox_lst[0] == 'OX':
                #                 ox = ox_lst[1]
                # af_sequences[ac] = {'mark': _mark, 'OX': ox}
                continue
            af_sz = len(af_sequences[ac])
            if af_sequences[ac]['seq'] == '':
                af_sequences[ac]['seq'] = line
                tg_i0 = af_sz  # len(af_sequences[ac])
            else:
                af_sequences[ac][tags[af_sz - tg_i0]] = line.replace('x', '-')
            continue
    af = {'data': af_sequences, 'metadata': {'tags': tags, 'name_tags': name_tags, 'statistics': None}}
    # print(af['metadata']['name_tags'], flush=True)
    return af


def annotated_fasta_gen_statistics(af):
    # counts = {'seq': 0, 'seg': 0}
    if len(af['data']) == 0:
        af['metadata']['statistics'] = None
        return
    af['metadata']['statistics'] = {}
    for tg in af['metadata']['tags']:
        af['metadata']['statistics'][tg] = {}
        af['metadata']['statistics'][tg]['seq'] = 0
        af['metadata']['statistics'][tg]['seg'] = 0
        for ac in af['data']:
            mask = str(af['data'][ac][tg])
            mask = mask.replace('-', '0')
            cnt = len([xx for xx in mask.split('0') if xx])
            if cnt > 0:
                af['metadata']['statistics'][tg]['seq'] += 1
                af['metadata']['statistics'][tg]['seg'] += cnt
        for cc in ['0', '1', '-']:
            af['metadata']['statistics'][tg][cc] = 0
            for ac in af['data']:
                for i in range(len(af['data'][ac][tg])):
                    if af['data'][ac][tg][i] == cc:
                        af['metadata']['statistics'][tg][cc] += 1


def get_string_stat(af):
    if not af['metadata']['statistics']:
        annotated_fasta_gen_statistics(af)
    _msg = "# Statistics:\n#\t---\ttag\tSeq#\tSeg#\t'0'\t'1'\t'-'"
    for tg in af['metadata']['tags']:
        _msg = _msg + f"\n#\tTAG\t{tg}"
        # print(af['metadata']['statistics'].keys(), flush=True)
        for cc in af['metadata']['statistics'][tg]:
            _msg = _msg + f"\t{af['metadata']['statistics'][tg][cc]:,}"
    return _msg


def annotated_fasta_save_fasta(af, f_name):
    with open(f_name, 'w') as fout:
        for ac in af['data']:
            ac_o = ac
            # for tg in af['data'][ac]:
            #     if tg not in af['metadata']['tags']:
            #         ac_o = f"{ac_o}|{tg}={af['data'][ac][tg]}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)


def annotated_fasta_save(af, f_name, data_name=None,  header_top=None, header_bottom=None):
    with open(f_name, 'w') as fout:
        if data_name:
            print(f"# dataset: {data_name}\n#", file=fout)
        if header_top:
            print(header_top, file=fout)
        print(f"# Sequences:\t{len(af['data']):,}", file=fout)
        print("#", file=fout)
        print("# Format:", file=fout)
        print("#\t>accession", file=fout)
        print("#\tAmino acid sequence", file=fout)
        for tg in af['metadata']['tags']:
            print(f"#\t{tg} annotation", file=fout)
        print("#", file=fout)
        print(get_string_stat(af), file=fout)
        if header_bottom:
            print('#', file=fout)
            print(header_bottom, file=fout)
        print('#', file=fout)
        for ac in af['data']:
            ac_o = ac
            for tg in af['data'][ac]:
                if tg in af['metadata']['name_tags']:
                    ac_o = f"{ac_o}|{tg}={af['data'][ac][tg]}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)
            for tg in af['metadata']['tags']:
                print(f"{af['data'][ac][tg]}", file=fout)
    return


def annotated_fasta_remove_ac_set(af, ac_set):
    for ac in ac_set:
        if ac in af['data']:
            del af['data'][ac]


def annotated_fasta_remove_tag(af, tags_out):
    tag_list = []
    for tg in af['metadata']['tags']:
        if tg not in tags_out:
            tag_list.append(tg)
    for ac in af['data']:
        for tg in af['metadata']['tags']:
            if tg in tags_out:
                del af['data'][ac][tg]
    af['metadata']['tags'] = tag_list


def annotated_fasta_rename_tag(af, old_tag, new_tag):
    for ii in range(len(af['metadata']['tags'])):
        if af['metadata']['tags'][ii] == old_tag:
            af['metadata']['tags'][ii] = new_tag
            break
    for ac in af['data']:
        af['data'][ac][new_tag] = af['data'][ac][old_tag]  # copy.deepcopy
        del af['data'][ac][old_tag]


def annotated_fasta_merge2(af1, af2):
    # print(af2['metadata']['tags'])
    merged2 = {'data': copy.deepcopy(af1['data']), 'metadata': {'tags': af1['metadata']['tags'], 'statistics': None}}
    for ac in af2['data']:
        if ac not in merged2['data']:
            print(f"AC {ac} not in merged2", flush=True)
            merged2['data'][ac] = copy.deepcopy(af2['data'][ac])
        else:
            if len(af2['data'][ac]['seq']) != len(merged2['data'][ac]['seq']):
                print(f"{ac} DELETED: seq size is different among the two input af data", flush=True)
                del merged2['data'][ac]
                continue
            if af2['data'][ac]['seq'] != merged2['data'][ac]['seq']:
                w_msg = f"{ac} sequences are not identical among the two input af data, same size."
                print(f"WARNING: {w_msg}", flush=True)
                if 'warnings' not in merged2:
                    merged2['warnings'] = {}
                merged2['warnings'][ac] = {'message': w_msg, 'alternative': af2['data'][ac]['seq']}
            for tg in af2['metadata']['tags']:
                if tg == 'seq':
                    continue
                if tg not in merged2['data'][ac].keys():
                    merged2['data'][ac][tg] = af2['data'][ac][tg]
                else:
                    m_tag = list(merged2['data'][ac][tg])
                    for i in range(len(af2['data'][ac][tg])):
                        if af2['data'][ac][tg][i] == '1':
                            m_tag[i] = '1'
                        elif af2['data'][ac][tg][i] == '0':
                            if m_tag[i] == '-':
                                m_tag[i] = '0'
                    merged2['data'][ac][tg] = ''.join(m_tag)
    return merged2


def annotated_fasta_merge_list(af_lst):
    if len(af_lst) == 0:
        return None
    elif len(af_lst) == 1:
        return af_lst[0]
    merged = copy.deepcopy(af_lst[0])
    for af in af_lst[1:]:
        merged = annotated_fasta_merge2(merged, af)
    return merged


def annotated_fasta_remove_empty(af, tag):
    if tag not in af['metadata']['tags']:
        return
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if '1' not in af['data'][ac][tag]:
            del af['data'][ac]


def annotated_fasta_remove_empty_all(af):
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        rmv = True
        for tg in af['metadata']['tags']:
            if '1' in af['data'][ac][tg]:
                rmv = False
                break
        if rmv:
            del af['data'][ac]


def load_cdhit_clusters(af, cdhit_clstr_file):
    cluster_dict = {}
    with open(cdhit_clstr_file, 'r') as fin:
        cluster = ''
        p_identity = ''
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue
            lst = line.split()
            if lst[0][0] == '>':
                cluster = lst[1]
                continue
            ac = lst[2].split('|')[0][1:]
            if lst[3] == '*':
                cluster_dict[cluster] = ac
                p_identity = '100.00%'
            elif lst[3] == 'at':
                p_identity = lst[4]
            else:
                print(f"ERROR ===============================================\t{line}", flush=True)
            if ac not in af['data']:
                print(ac, flush=True)
                continue
            af['data'][ac]['cluster'] = int(cluster)
            af['data'][ac]['p_identity'] = p_identity
    ac_list = list(af['data'].keys())
    # print("----", len(ac_list), flush=True)
    # exit(0)
    for ac in ac_list:
        if 'cluster' not in af['data'][ac]:
            del af['data'][ac]
        # af['data'][ac]['c_center'] = cluster_dict[str(af['data'][ac]['cluster'])]
    af['metadata']['name_tags'].append('cluster')
    af['metadata']['name_tags'].append('p_identity')
    af['metadata']['name_tags'].append('c_center')
    pass


def annotated_fasta_load_fasta(in_file, org_seq=True):
    ret_dict = {}
    if not os.path.exists(in_file):
        print(f"Input fasta file not found: {in_file}", flush=True)
        exit(0)
    with open(in_file, 'r') as fin:
        ac = ''
        seq = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                if len(seq) > 5 and len(ac) > 4:
                    if not org_seq:
                        seq = seq.replace('B', 'A').replace('Z', 'A')
                        seq = seq.replace('X', 'A').replace('U', 'A')
                    ret_dict[ac] = seq
                ac = line[1:]
                seq = ''
            else:
                seq = seq + line
        if len(seq) > 5 and len(ac) > 4:
            if not org_seq:
                seq = seq.replace('B', 'A').replace('Z', 'A')
                seq = seq.replace('X', 'A').replace('U', 'A')
            ret_dict[ac] = seq
    return ret_dict
