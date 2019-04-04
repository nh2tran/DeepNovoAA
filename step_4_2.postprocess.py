import csv


# extract denovo peptides from the PSMs of PEAKS X DB search round_2
psm_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/step_4.DB search psm round_2_FDR_1%.csv"
output_denovo_peptide_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/step_4.output_peptide_list"


def drop_mod_peaks(peptide):
  peptide = peptide.replace("M(+15.99)", "M")
  peptide = peptide.replace("N(+.98)", "N")
  peptide = peptide.replace("Q(+.98)", "Q")
  return peptide


if __name__ == '__main__':

  denovo_peptide_set = set()
  with open(psm_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      peptide = drop_mod_peaks(row['Peptide'])
      accession = drop_mod_peaks(row['Accession'])
      if accession == 'DENOVO':
        denovo_peptide_set.add(peptide)

  with open(output_denovo_peptide_file, 'w') as output_handle:
    for peptide in denovo_peptide_set:
        output_handle.write(peptide + '\n')

  num_denovo_peptides = len(denovo_peptide_set)
  print("num_denovo_peptides =", num_denovo_peptides)
