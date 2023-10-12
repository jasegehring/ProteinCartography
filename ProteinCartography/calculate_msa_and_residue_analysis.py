## Generate a fasta file prior to calculating the multiple sequence alignment
import csv
import sys
import argparse
import subprocess
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", nargs='+', required=True, help="Path to the aggregated features file."
    )
    parser.add_argument(
        "-o", "--output", nargs='+', required=True, help="Path of output fasta file."
    )
    args = parser.parse_args()

    return args

# write a fasta file from the aggregated tsv file
def tsv_to_fasta(input_filename, output_filename):
    with open(input_filename, 'r') as tsv_file, open(output_filename, 'w') as fasta_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        next(reader)  # skip the header
        for row in reader:
            protid = row[0]
            sequence = row[11]
            fasta_file.write(f">{protid}\n{sequence}\n")

# use famsa to calculate the multiple sequence alignment
def calculate_MSA(input_filename, output_filename):
    # run famsa
    subprocess.run(["famsa", "-t", "0", input_filename, output_filename])

def get_aa_frequencies(msa_aln, aa_frequencies):
    """
    Input: multiple sequence alignment
    Function: 
        1. Produce a Pandas DataFrame with amino acid frequences for each position in the MSA
    """
    seqs = [(record.id, list(record.seq)) for record in SeqIO.parse(msa_aln, "fasta")]
    MSA_df = pd.DataFrame.from_dict(dict(seqs))
    #delete seqs dictionary from memory
    seqs=None
    
    frequencies_df_temp = pd.DataFrame()
    for index, row in MSA_df.iterrows():
        frequencies_df_temp = frequencies_df_temp.append(ProteinAnalysis(row.str.cat()).count_amino_acids(), ignore_index=True)
    frequencies_df_temp.to_csv(aa_frequencies, sep='\t', index=True, header=True)

# run this if called from the interpreter
def main():
    args = parse_args()
    input_files = args.input
    output_files = args.output

    tsv_to_fasta(input_files[0], output_files[0])
    calculate_MSA(output_files[0], output_files[1])
    get_aa_frequencies(output_files[1], output_files[2])


# check if called from interpreter
if __name__ == "__main__":
    main()