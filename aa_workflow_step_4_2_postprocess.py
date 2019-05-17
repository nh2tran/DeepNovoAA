from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import csv


def drop_mod_peaks(peptide):
  peptide = peptide.replace("M(+15.99)", "M")
  peptide = peptide.replace("N(+.98)", "N")
  peptide = peptide.replace("Q(+.98)", "Q")
  return peptide


def step_4_2_postprocess(psm_file, output_denovo_peptide_file):
  """Extract denovo peptides from the PSMs of PEAKS X DB search round 2.

     Usage:
       psm_file = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/aa_workflow.step_4.psm.csv"
       output_denovo_peptide_file = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/aa_workflow.step_4.output_peptide_list"
  """

  print("".join(["="] * 80)) # section-separating line
  print("step_4_2_postprocess()")

  print("psm_file =", psm_file)
  print("output_denovo_peptide_file =", output_denovo_peptide_file)

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
