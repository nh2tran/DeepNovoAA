import csv
from Bio import SeqIO
from Bio.SeqIO import FastaIO


# remove denovo peptides that exist in the database fasta file
# combine db and denovo into a peptide list file for PEAKS X DB search round 2
denovo_file = "data.training/aa.hla.bassani.nature_2016.mel_15/feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"
db_fasta_file = "data/uniprot_sprot.human.plus_contaminants.fasta"
labeled_feature_file = "data.training/aa.hla.bassani.nature_2016.mel_15/feature.csv.labeled"
peptide_list_fasta = "data.training/aa.hla.bassani.nature_2016.mel_15/step_4.peptide_list.fasta"


def drop_mod(peptide):
  peptide = peptide.replace("M(Oxidation)", "M")
  peptide = peptide.replace("N(Deamidation)", "N")
  peptide = peptide.replace("Q(Deamidation)", "Q")
  return peptide


def drop_mod_peaks(peptide):
  peptide = peptide.replace("M(+15.99)", "M")
  peptide = peptide.replace("N(+.98)", "N")
  peptide = peptide.replace("Q(+.98)", "Q")
  return peptide


def change_I_to_L(string):
    return string.replace('I', 'L')


if __name__ == '__main__':

  denovo_peptide_set = set()
  with open(denovo_file, 'r') as fr:
      reader = csv.reader(fr, delimiter='\t')
      names = next(reader)
      seq_index = names.index('predicted_sequence')
      for line in reader:
          if not line[seq_index]:
              continue
          peptide = line[seq_index]
          peptide = drop_mod(peptide)
          peptide = ''.join(peptide.split(','))
          if peptide in denovo_peptide_set:
              continue
          else:
              denovo_peptide_set.add(peptide)
  print("Number of top-scoring denovo peptides: {}".format(len(denovo_peptide_set)))

  with open(db_fasta_file, 'r') as input_fasta_handle:
      record_list = list(SeqIO.parse(input_fasta_handle, "fasta"))
      print("Number of protein sequences: ", len(record_list))
  human_protein_list = [str(record.seq) for record in record_list]

  # remove denovo peptides that exist in the database fasta file
  to_L_protein_list = [change_I_to_L(protein) for protein in human_protein_list]
  pure_denovo_seq_set = set()
  for i, peptide in enumerate(denovo_peptide_set):
      peptide_string = change_I_to_L(peptide)
      indb = False
      for protein in to_L_protein_list:
          if peptide_string in protein:
              indb = True
              break
      if not indb:
          pure_denovo_seq_set.add(peptide)
      if i % 1000 == 0:
          print("processing {}".format(i))
  print("Number of denovo peptides not in database: {}".format(len(pure_denovo_seq_set)))

  db_peptide_set = set()
  with open(labeled_feature_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      peptide = drop_mod_peaks(row['seq'])
      db_peptide_set.add(peptide)

  with open(peptide_list_fasta, 'w') as output_handle:
    counter = 0
    for peptide in db_peptide_set:
      counter += 1
      output_handle.write(">DB|db_{}\n".format(counter))
      output_handle.write(peptide + '\n')
    counter = 0
    for peptide in pure_denovo_seq_set:
      counter += 1
      output_handle.write(">DENOVO|denovo_{}\n".format(counter))
      output_handle.write(''.join(peptide) + '\n')

  num_db_peptides = len(db_peptide_set)
  num_denovo_peptides = len(pure_denovo_seq_set)
  print("num_db_peptides =", num_db_peptides)
  print("num_denovo_peptides =", num_denovo_peptides)
