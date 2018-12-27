from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math
import os
import random
import sys
import time
import re

import csv
import numpy as np
random.seed(0)
np.random.seed(0)

from Bio import SeqIO
from Bio.SeqIO import FastaIO

import deepnovo_config






# write multi-line fasta file into single-line format
def write_fasta_1line(input_fasta_file, output_fasta_file):
  with open(input_fasta_file, "r")  as handle:
    record_list = list(SeqIO.parse(handle, "fasta"))
    print(input_fasta_file)
    print("Number of protein sequences: ", len(record_list))
  with open(output_fasta_file, "w") as handle:
    fasta_writer = FastaIO.FastaWriter(handle, wrap=None)
    fasta_writer.write_file(record_list)

# ~ input_fasta_file = "data/uniprot.human_all_isoforms.fasta"
# ~ output_fasta_file = input_fasta_file + ".1line"
# ~ write_fasta_1line(input_fasta_file, output_fasta_file)


# calculate peptide mass = N-terminus + amino acids + C-terminus
def compute_peptide_mass(peptide):
  """TODO(nh2tran): docstring.
  """

  #~ print("".join(["="] * 80)) # section-separating line ===
  #~ print("WorkerDB: _compute_peptide_mass()")

  peptide_mass = (deepnovo_config.mass_N_terminus
                  + sum(deepnovo_config.mass_AA[aa] for aa in peptide)
                  + deepnovo_config.mass_C_terminus)

  return peptide_mass

# ~ peptide = 'AAAAAAALQAK'
# ~ print(peptide)
# ~ print(compute_peptide_mass(peptide))


# calculate ppm of precursor_mz against peptide_mz
def check_mass_shift(input_feature_file):
  print("check_mass_shift()")

  print("input_feature_file = ", os.path.join(input_feature_file))

  with open(input_feature_file, mode="r") as input_handle:
    # header line
    line = input_handle.readline()
    # first feature
    line = input_handle.readline()
    precursor_ppm_list = []
    # ~ test_handle = open("test.txt", mode="w")
    while line:
      line_split = re.split(',|\r|\n', line)

      # parse raw peptide sequence
      raw_sequence = line_split[deepnovo_config.col_raw_sequence]
      raw_sequence_len = len(raw_sequence)
      peptide = []
      index = 0
      while index < raw_sequence_len:
        if raw_sequence[index] == "(":
          if peptide[-1] == "C" and raw_sequence[index:index + 8] == "(+57.02)":
            peptide[-1] = "C(Carbamidomethylation)"
            index += 8
          elif peptide[-1] == 'M' and raw_sequence[index:index + 8] == "(+15.99)":
            peptide[-1] = 'M(Oxidation)'
            index += 8
          elif peptide[-1] == 'N' and raw_sequence[index:index + 6] == "(+.98)":
            peptide[-1] = 'N(Deamidation)'
            index += 6
          elif peptide[-1] == 'Q' and raw_sequence[index:index + 6] == "(+.98)":
            peptide[-1] = 'Q(Deamidation)'
            index += 6
          else:  # unknown modification
            print("ERROR: unknown modification!")
            print("raw_sequence = ", raw_sequence)
            sys.exit()
        else:
          peptide.append(raw_sequence[index])
          index += 1

      precursor_mz = float(line_split[deepnovo_config.col_precursor_mz])
      precursor_charge = float(line_split[deepnovo_config.col_precursor_charge])
      peptide_mass = compute_peptide_mass(peptide)
      peptide_mz = (peptide_mass + precursor_charge * deepnovo_config.mass_H) / precursor_charge
      precursor_ppm = (precursor_mz - peptide_mz) / peptide_mz * 1e6
      # ~ print('precursor_mz = ', precursor_mz)
      # ~ print('peptide_mz = ', peptide_mz)
      # ~ print('precursor_ppm = ', precursor_ppm)
      # ~ print(abc)
      # ~ test_row = '\t'.join(['{0:.6f}'.format(x) for x in [peptide_mz, precursor_mz, precursor_ppm]])
      # ~ print(test_row, file=test_handle, end="\n")
      precursor_ppm_list.append(precursor_ppm)
      line = input_handle.readline()
    # ~ test_handle.close()
    print(np.mean(precursor_ppm_list))

# ~ input_feature_file = "identified_features_corrected.csv"
# ~ check_mass_shift(input_feature_file)


# randomly partition a feature file into train/valid/test files
def partition_feature_file(input_feature_file, prob):
  print("partition_feature_file()")

  print("input_feature_file = ", os.path.join(input_feature_file))
  print("prob = ", prob)

  output_file_train = input_feature_file + ".train"
  output_file_valid = input_feature_file + ".valid"
  output_file_test = input_feature_file + ".test"

  with open(input_feature_file, mode="r") as input_handle:
    with open(output_file_train, mode="w") as output_handle_train:
      with open(output_file_valid, mode="w") as output_handle_valid:
        with open(output_file_test, mode="w") as output_handle_test:
          counter = 0
          counter_train = 0
          counter_valid = 0
          counter_test = 0
          # header line
          line = input_handle.readline()
          output_handle_train.write(line)
          output_handle_valid.write(line)
          output_handle_test.write(line)
          # first feature
          line = input_handle.readline()
          while line:
            counter += 1
            set_num = np.random.choice(a=3, size=1, p=prob)
            if set_num == 0:
              output_handle = output_handle_train
              counter_train += 1
            elif set_num == 1:
              output_handle = output_handle_valid
              counter_valid += 1
            else:
              output_handle = output_handle_test
              counter_test += 1
            output_handle.write(line)
            line = input_handle.readline()

  input_handle.close()
  output_handle_train.close()
  output_handle_valid.close()
  output_handle_test.close()

  print("counter = {0:d}".format(counter))
  print("counter_train = {0:d}".format(counter_train))
  print("counter_valid = {0:d}".format(counter_valid))
  print("counter_test = {0:d}".format(counter_test))


# partition a feature file into train/valid/test files with NO common peptides
def partition_feature_file_nodup(input_feature_file, prob):
  print("partition_feature_file_nodup()")

  print("input_feature_file = ", os.path.join(input_feature_file))
  print("prob = ", prob)

  output_file_train = input_feature_file + ".train" + ".nodup"
  output_file_valid = input_feature_file + ".valid" + ".nodup"
  output_file_test = input_feature_file + ".test" + ".nodup"

  peptide_train_list = []
  peptide_valid_list = []
  peptide_test_list = []
  
  with open(input_feature_file, mode="r") as input_handle:
    with open(output_file_train, mode="w") as output_handle_train:
      with open(output_file_valid, mode="w") as output_handle_valid:
        with open(output_file_test, mode="w") as output_handle_test:
          counter = 0
          counter_train = 0
          counter_valid = 0
          counter_test = 0
          counter_unique = 0
          # header line
          line = input_handle.readline()
          output_handle_train.write(line)
          output_handle_valid.write(line)
          output_handle_test.write(line)
          # first feature
          line = input_handle.readline()
          while line:
            counter += 1
            # check if the peptide already exists in any of the three lists
            # if yes, this new feature will be assigned to that list
            peptide = re.split(',|\r|\n', line)[deepnovo_config.col_raw_sequence]
            if (peptide in peptide_train_list):
              output_handle = output_handle_train
              counter_train += 1
            elif (peptide in peptide_valid_list):
              output_handle = output_handle_valid
              counter_valid += 1
            elif (peptide in peptide_test_list):
              output_handle = output_handle_test
              counter_test += 1
            # if not, this new peptide and its spectrum will be randomly assigned
            else:
              counter_unique += 1
              set_num = np.random.choice(a=3, size=1, p=prob)
              if set_num == 0:
                peptide_train_list.append(peptide)
                output_handle = output_handle_train
                counter_train += 1
              elif set_num == 1:
                peptide_valid_list.append(peptide)
                output_handle = output_handle_valid
                counter_valid += 1
              else:
                peptide_test_list.append(peptide)
                output_handle = output_handle_test
                counter_test += 1
            output_handle.write(line)
            line = input_handle.readline()

  input_handle.close()
  output_handle_train.close()
  output_handle_valid.close()
  output_handle_test.close()

  print("counter = {0:d}".format(counter))
  print("counter_train = {0:d}".format(counter_train))
  print("counter_valid = {0:d}".format(counter_valid))
  print("counter_test = {0:d}".format(counter_test))
  print("counter_unique = {0:d}".format(counter_unique))

#~ folder_path = "data.training/dia.pecan.hela.2018_03_29/"
#~ partition_feature_file_nodup(input_feature_file=folder_path + "training_5mz_4to7.feature.csv",
                             #~ prob=[0.90, 0.05, 0.05])


# merge multiple mgf files into one, adding fraction ID to scan ID
def cat_file_mgf(input_file_list, fraction_list, output_file):
  print("cat_file_mgf()")
  
  # iterate over mgf files and their lines
  counter = 0
  with open(output_file, mode="w") as output_handle:
    for input_file, fraction in zip(input_file_list, fraction_list):
      print("input_file = ", os.path.join(input_file))
      with open(input_file, mode="r") as input_handle:
        for line in input_handle.readlines():
          if "SCANS=" in line: # a spectrum found
            counter += 1
            scan = re.split('=|\n|\r', line)[1]
            # re-number scan id
            output_handle.write("SCANS=F{0}:{1}\n".format(fraction, scan))
          else:
            output_handle.write(line)
  print("output_file = {0:s}".format(output_file))
  print("counter = {0:d}".format(counter))


# merge multiple feature files into one, adding fraction ID to feature & scan ID
def cat_file_feature(input_file_list, fraction_list, output_file):
  print("cat_file_feature()")
  
  # read and write header line
  csv_reader = csv.DictReader(open(input_file_list[0]))
  csv_writer = csv.DictWriter(open(output_file, mode='w'), csv_reader.fieldnames)
  csv_writer.writeheader()
  # iterate over feature files and their rows
  counter = 0
  for input_file, fraction in zip(input_file_list, fraction_list):
    print("input_file = ", os.path.join(input_file))
    csv_reader = csv.DictReader(open(input_file))
    for row in csv_reader:
      counter += 1
      # add fraction to feature id
      feature_id = row['spec_group_id']
      feature_id = "F" + str(fraction) + ":" + feature_id
      row['spec_group_id'] = feature_id
      # add fraction to scan id
      scan_list = re.split(';', row['scans'])
      scan_list = ["F" + str(fraction) + ":" + x for x in scan_list]
      row['scans'] = ";".join(scan_list)
      # join the line back together and write to output
      csv_writer.writerow(row)
  print("output_file = {0:s}".format(output_file))
  print("counter = {0:d}".format(counter))

folder_path = "data.training/aa.hla.bassani.nature_2016.mel_15/"
fraction_list = range(0, 15+1)
cat_file_mgf(
    input_file_list=[folder_path + "export_" + str(i) + ".mgf" for i in fraction_list],
    fraction_list=fraction_list,
    output_file=folder_path + "spectrum.mgf")
cat_file_feature(
    input_file_list=[folder_path + "export_" + str(i) + ".csv" for i in fraction_list],
    fraction_list=fraction_list,
    output_file=folder_path + "feature.csv")

