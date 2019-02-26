import re
import csv
from Bio import SeqIO
from Bio.SeqIO import FastaIO

denovo_file = "data.training/aa.hla.bassani.nature_2016.mel_15/feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"
denovo_peptide_file = "data.training/aa.hla.bassani.nature_2016.mel_15/step4.peptide"
fasta_file = "data/uniprot_sprot.human.plus_contaminants.fasta"


def drop_mod(peptide_str):
    peptide_str = re.sub("M\(Oxidation\)", "M", peptide_str)
    peptide_str = re.sub("N\(Deamidation\)", "N", peptide_str)
    peptide_str = re.sub("Q\(Deamidation\)", "Q", peptide_str)
    return peptide_str


def change_I_to_L(string):
    return string.replace('I', "L")



if __name__ == '__main__':
    denovo_peptide_set = set()
    with open(denovo_file, 'r') as fr:
        reader = csv.reader(fr, delimiter='\t')
        names = next(reader)
        seq_index = names.index('predicted_sequence')
        for line in reader:
            if not line[seq_index]:
                continue
            raw_seq = line[seq_index]
            raw_seq = drop_mod(raw_seq)
            if raw_seq in denovo_peptide_set:
                continue
            else:
                denovo_peptide_set.add(raw_seq)
    print("top scoring denovo peptides: {}".format(len(denovo_peptide_set)))
    with open(fasta_file, 'r') as input_fasta_handle:
        record_list = list(SeqIO.parse(input_fasta_handle, "fasta"))
        print("Number of protein sequences: ", len(record_list))
    human_protein_list = [str(record.seq) for record in record_list]

    # with open("./I_to_L/human_protein_origin.cc.txt", 'w') as f:
    #     for protein in human_protein_list:
    #         protein = ','.join(list(protein))
    #         f.write(protein + '\n')
    to_L_protein_list = [change_I_to_L(protein) for protein in human_protein_list]
    pure_denovo_seq_set = set()
    for i, raw_seq in enumerate(denovo_peptide_set):
        peptide_string = change_I_to_L(''.join(raw_seq.split(',')))
        indb = False
        for protein in to_L_protein_list:
            if peptide_string in protein:
                indb = True
                break
        if not indb:
            pure_denovo_seq_set.add(raw_seq)
        if i % 1000 == 0:
            print("processing {}".format(i))
    print("num pure denovo peptide: {}".format(len(pure_denovo_seq_set)))
    with open(denovo_peptide_file, 'w') as fw:
        for peptide in pure_denovo_seq_set:
            fw.write(''.join(peptide.split(',')) + '\n')
