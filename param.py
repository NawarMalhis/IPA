# input_fasta = 'input.fasta'
# output_path = "output/"
aff_path = '/home/nmalhis/Papers/AFF/'
af2_cif_path = "cif/"
download_cif = True

models_path = 'Models/'

f_used = {'C': ['PDB', 'IDR', 'Linker', 'P_bind', 'N_bind'],
          'CP': ['w1', 'w2', 'PDB', 'IDR', 'Linker', 'P_bind', 'N_bind'],
          'ACP': ['RSA', 'pLDDt', 'w1', 'w2', 'PDB', 'IDR', 'Linker', 'P_bind', 'N_bind']
          }

mdl_data = {'LIP_2': {'FU': f_used['CP'], 'AF': '--', 'Name': 'All_bind', 'Tag': 'LIP'},
            'P_bind_2': {'FU': f_used['CP'], 'AF': '--', 'Name': 'Protein_bind', 'Tag': 'P_bind'},
            'N_bind_1': {'FU': f_used['C'], 'AF': '--', 'Name': 'Nucleotide_bind', 'Tag': 'N_bind'},
            'P_bind_3': {'FU': f_used['ACP'], 'AF': 'AF', 'Name': 'AF2_Protein_bind', 'Tag': 'P_bind'},
            'Linker_3': {'FU': f_used['ACP'], 'AF': 'AF', 'Name': 'AF2_Linker', 'Tag': 'Linker'},
            }

mdl_list = ['P_bind_2', 'N_bind_1', 'LIP_2']
af2_mdl_list = ['P_bind_3', 'Linker_3']
