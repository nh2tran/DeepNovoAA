import csv


labeled_feature_file = "data.training/aa.hla.bassani.nature_2016.mel_15/feature.csv.labeled"
denovo_peptide_file = "data.training/aa.hla.bassani.nature_2016.mel_15/step4.output_peptide_list"
output_fasta_file = "data.training/aa.hla.bassani.nature_2016.mel_15/step5.peptide_list.fasta"


def drop_mod_peaks(peptide):
  peptide = peptide.replace("M(+15.99)", "M")
  peptide = peptide.replace("N(+.98)", "N")
  peptide = peptide.replace("Q(+.98)", "Q")
  return peptide


if __name__ == '__main__':

  db_peptide_set = set()
  with open(labeled_feature_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      peptide = drop_mod_peaks(row['seq'])
      db_peptide_set.add(peptide)

  denovo_peptide_set = set()
  with open(denovo_peptide_file, 'r') as input_handle:
    for line in input_handle.readlines():
      peptide = line.strip()
      denovo_peptide_set.add(peptide)

  with open(output_fasta_file, 'w') as output_handle:
    counter = 0
    for peptide in db_peptide_set:
      counter += 1
      output_handle.write(">DB|db_{}\n".format(counter))
      output_handle.write(peptide + '\n')
    counter = 0
    for peptide in denovo_peptide_set:
      counter += 1
      output_handle.write(">DENOVO|denovo_{}\n".format(counter))
      output_handle.write(peptide + '\n')

  num_db_peptides = len(db_peptide_set)
  num_denovo_peptides = len(denovo_peptide_set)
  print("num_db_peptides =", num_db_peptides)
  print("num_denovo_peptides =", num_denovo_peptides)
