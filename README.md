# IPA: Probabilistic Annotations of Protein Sequences for Intrinsically Disordered Features

This IDR Probabilistic Annotation (IPA) platform predicts 'Linker' regions and 'nucleic', 'protein', and 'all' (protein or
nucleic) IDR binding sites within protein amino acid sequences.

### Please reference the following publications:
Malhis N, "Probabilistic Annotations of Protein Sequences for Intrinsically Disordered Features" bioRxiv (2024). [doi: https://doi.org/10.1101/2024.12.18.629275].

## Minimum Hardware Requirements
OS: Linux (tested on Ubuntu and CentOS7).  
RAM: 8 GB minimum, 16 GB recommended.  
CPU: Multicore with 4+ cores recommended.  

## To install:
```bash
# Clone IPA
git clone https://github.com/NawarMalhis/IPA.git
# Clone the "annotated fasta format" library:	
git clone https://github.com/NawarMalhis/AFF.git
# Change directory:	
cd IPA
# Update the path to the AFF (annotated fasta format) folder in param.py
aff_path = '/xxx/xxx/AFF/'
# Create an ipa_cpu or ipa_cuda environment:
conda env create -f ipa_cuda.yml
```

## Input data:
All input sequences can be in a single fasta file, so each sequence has a fasta header line with a unique accession and a
sequence line. Example of a sequence in a fasta format:  
```bash
>P07766
MQSGTHWRVLGLCLLSVGVWGQDGNEEMGGITQTPYKVSISGTTVILTCPQYPGSEILWQHNDKNIGG…
>P05067
MLPGLALLLLAAWTARALEVPTDGNAGLLAEPQIAMFCGRLNMHMNVQNGKWDSDPSGTKTCIDTKEG…
>P38398
MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDI…
```


## To run:
```bash
# Activate the ipa environment:
conda activate ipa_cpu

python3 ipa.py –p PATH –in INPUT_FILE [-af2] [-cpu]

-af2 is used to run IPA-AF2
-cpu can be used to force running on a CPU in a Cuda environment.
INPUT_FILE: is the input fasta file. This input file can be located in any accessible directory.
PATH: is the path (directory) of the INPUT_FILE
IPA will create a directory "results" inside PATH to store all IPA predictions.
```

## Example:
```bash
python3 ipa.py -p test_data/ -in input.fasta
python3 ipa.py -p test_data/ -in input.fasta –af2
```
