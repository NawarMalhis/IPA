# IPA: IDR Probabilistic Annotation for amino acid sequences
#
# Developed by Nawar Malhis, 2024
# Michael Smith Laboratories,
# The University of British Columbia

from lib import *
from datetime import datetime
from param import *
from annotated_fasta import annotated_fasta_load_fasta


if __name__ == '__main__':
    arg = get_arguments()
    input_fasta = f"{arg.path}/{arg.in_file}"
    output_path = f"{arg.path}/output/"
    create_out_dirs(_path=arg.path, af2=arg.AlphaFold2)
    if arg.cpu:
        device = torch.device('cpu')
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print("\n\tIDR Probabilistic Annotation (IPA)", flush=True)
    print(f"\tRunning on {device}\n", flush=True)
    models = load_models(device=device, af2=arg.AlphaFold2, verbose=False)
    priors_dict = load_priors()
    fasta = annotated_fasta_load_fasta(input_fasta)

    print(f"Sequence count:\t{len(fasta):,}")
    if arg.AlphaFold2:
        t_file = f"{output_path}ipa-af2-timings.csv"
    else:
        t_file = f"{output_path}ipa-timings.csv"

    # microsecond
    with open(t_file, 'w') as fout:
        print(f"# Protein_Accession , Execution_Time_In_Microseconds", file=fout, flush=True)
        for ac in fasta:
            start_time = datetime.now()
            print(f"Processing: {ac}", flush=True)
            if arg.AlphaFold2:
                score_af2_ensemble(output_path, in_ac=ac, in_seq=fasta[ac], models=models, priors_dict=priors_dict, device=device)
            else:
                score_ensemble(output_path, in_ac=ac, in_seq=fasta[ac], models=models, priors_dict=priors_dict, device=device)
            end_time = datetime.now()
            td = end_time - start_time
            print(f"{ac} , {round(td.total_seconds() * 1000)}", file=fout, flush=True)
