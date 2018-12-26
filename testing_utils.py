from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math
import os
import random
import sys
import time
import re

import numpy as np

from Bio import SeqIO
from Bio.SeqIO import FastaIO

#~ from tensorflow.models.rnn.translate import data_utils
import deepnovo_config

random.seed(0)
np.random.seed(0)







# Write fasta file in single-line format
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
    print(np.mean(precursor_ppm_list))
    # ~ test_handle.close()

# ~ input_feature_file = "identified_features_corrected.csv"
# ~ check_mass_shift(input_feature_file)


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


#~ folder_path = "data.training/dia.pecan.hela.2018_03_29/"
#~ partition_feature_file_nodup(input_feature_file=folder_path + "training_5mz_4to7.feature.csv",
                             #~ prob=[0.90, 0.05, 0.05])
#~ partition_feature_file(input_feature_file=folder_path + "testing.feature.csv",
                       #~ prob=[0.10, 0.90, 0.00])






def cat_file_mgf(input_file_list, fraction_list, output_file):
  print("cat_file_mgf()")
  
  counter = 0
  with open(output_file, mode="w") as output_handle:
    for index, input_file in enumerate(input_file_list):
      print("input_file = ", os.path.join(input_file))
      with open(input_file, mode="r") as input_handle:
          line = input_handle.readline()
          while line:
            if "SCANS=" in line: # a spectrum found
              counter += 1
              scan = re.split('=|\n', line)[1]
              # re-number scan id
              output_handle.write("SCANS=F{0}:{1}\n".format(fraction_list[index], scan))
            else:
                output_handle.write(line)
            line = input_handle.readline()

  print("output_file = {0:s}".format(output_file))
  print("counter = {0:d}".format(counter))


def cat_file_feature(input_file_list, fraction_list, output_file):
  print("cat_file_feature()")
  
  counter = 0
  with open(output_file, mode="w") as output_handle:
    for index, input_file in enumerate(input_file_list):
      print("input_file = ", os.path.join(input_file))
      with open(input_file, mode="r") as input_handle:
          header_line = input_handle.readline()
          if index == 0 :
            output_handle.write(header_line)
          line = input_handle.readline()
          while line:
            counter += 1
            line = re.split(',|\r|\n', line)[:deepnovo_config.col_num]
            # add fraction to feature id
            feature_id = line[deepnovo_config.col_feature_id]
            feature_id = "F" + str(fraction_list[index]) + ":" + feature_id
            line[deepnovo_config.col_feature_id] = feature_id
            # add fraction to scan id
            scan_list = re.split(';', line[deepnovo_config.col_scan_list])
            scan_list = ["F" + str(fraction_list[index]) + ":" + x for x in scan_list]
            line[deepnovo_config.col_scan_list] = ";".join(scan_list)
            # join the line back together and write to output
            output_handle.write(",".join(line) + "\n")
            line = input_handle.readline()

  print("output_file = {0:s}".format(output_file))
  print("counter = {0:d}".format(counter))


#~ folder_path = "data.training/dia.pecan.hela.2018_03_29/"
#~ fraction_list = range(4, 7+1)
#~ cat_file_mgf(input_file_list=[folder_path + str(i) + "_frac.spectrum.mgf"
                              #~ for i in fraction_list],
             #~ fraction_list=fraction_list,
             #~ output_file=folder_path + "testing_jurkat_oxford.spectrum.mgf")
#~ cat_file_feature(input_file_list=[folder_path + str(i) + "_frac.feature.csv"
                                  #~ for i in fraction_list],
                 #~ fraction_list=fraction_list,
                 #~ output_file=folder_path + "training_5mz_4to7.feature.csv")
#~ cat_file_feature(input_file_list=[folder_path + str(i) + "_frac.unlabeled.csv"
                                  #~ for i in fraction_list],
                 #~ fraction_list=fraction_list,
                 #~ output_file=folder_path + "testing_jurkat_oxford.unlabeled.csv")


# mixed training
#~ fraction_list = range(1, 3+1)
#~ cat_file_mgf(
    #~ input_file_list=["data.training/dia.urine.2018_04_23/training_pain.spectrum.mgf",
                     #~ "data.training/dia.abrf.2018_03_27/training.spectrum.mgf",
                     #~ "data.training/dia.pecan.hela.2018_03_29/training.spectrum.mgf"],
    #~ fraction_list=fraction_list,
    #~ output_file="data.training/dia.urine_abrf_hela.2018_05_31/training.spectrum.mgf")
#~ cat_file_feature(
    #~ input_file_list=["data.training/dia.urine.2018_04_23/training_pain.feature.csv",
                     #~ "data.training/dia.abrf.2018_03_27/training.feature.csv",
                     #~ "data.training/dia.pecan.hela.2018_03_29/training.feature.csv"],
    #~ fraction_list=fraction_list,
    #~ output_file="data.training/dia.urine_abrf_hela.2018_05_31/training.feature.csv")


# mixed training
pain_prefix = ['DIA_A1',
                'DIA_A10',
                'DIA_A11',
                'DIA_A2',
                'DIA_A3',
                'DIA_A4',
                'DIA_A7',
                'DIA_A8',
                'DIA_A9',
                'DIA_B1',
                'DIA_B10',
                'DIA_B11',
                'DIA_B12',
                'DIA_B2',
                'DIA_B3',
                'DIA_B4',
                'DIA_B5',
                'DIA_B6',
                'DIA_B7',
                'DIA_B8',
                'DIA_B9',
                'DIA_C1',
                'DIA_C11',
                'DIA_C3',
                'DIA_C4',
                'DIA_C6',
                'DIA_C7',
                'DIA_C9',
                'DIA_D10',
                'DIA_D12',
                'DIA_D3',
                'DIA_D4',
                'DIA_D5',
                'DIA_D6',
                'DIA_D7',
                'DIA_D8',
                'DIA_E1',
                'DIA_E10',
                'DIA_E12',
                'DIA_E2',
                'DIA_E3',
                'DIA_E4',
                'DIA_E7',
                'DIA_E8',
                'DIA_E9',
                'DIA_F12',
                'DIA_F2',
                'DIA_F3',
                'DIA_F4',
                'DIA_F5',
                'DIA_F7',
                'DIA_F8',
                'DIA_G11',
                'DIA_G12',
                'DIA_G3',
                'DIA_G4',
                'DIA_G5',
                'DIA_G6',
                'DIA_G7',
                'DIA_G8',
                'DIA_G9',
                'DIA_H1',
                'DIA_H2',
                'DIA_H3',]
oc_prefix = ['DIA_A12',
                #~ 'DIA_A6',
                'DIA_C12',
                #~ 'DIA_C5',
                #~ 'DIA_C8',
                'DIA_D11',
                #~ 'DIA_D9',
                'DIA_F9',
                #~ 'DIA_G1',
                'DIA_G2',
                #~ 'DIA_H4',
                'DIA_H5',]
uti_prefix = ['DIA_A5',
              #~ 'DIA_C10',
              'DIA_C2',
              #~ 'DIA_D1',
              'DIA_D2',
              #~ 'DIA_E11',
              'DIA_E5',
              #~ 'DIA_E6',
              'DIA_F10',
              #~ 'DIA_F6',
              'DIA_G10',]
folder_path = "data.training/dia.urine.2018_04_23/"
fraction_list = range(1, len(pain_prefix)+1)
#~ cat_file_mgf(
    #~ input_file_list=[folder_path + x + '.mgf' for x in uti_prefix],
    #~ fraction_list=fraction_list,
    #~ output_file="data.training/dia.urine.2018_04_23/testing_uti.spectrum.mgf")
#~ cat_file_feature(
    #~ input_file_list=[folder_path + x + '_labelled.csv' for x in uti_prefix],
    #~ fraction_list=fraction_list,
    #~ output_file="data.training/dia.urine.2018_04_23/testing_uti.feature.csv")
#~ cat_file_feature(
    #~ input_file_list=[folder_path + x + '.csv' for x in pain_prefix],
    #~ fraction_list=fraction_list,
    #~ output_file="data.training/dia.urine.2018_04_23/training_pain.unlabeled.csv")





def generate_labeled_mgf(fraction_list, input_mgf_list, input_psm_file):
  print("generate_labeled_mgf()")
  
  ### store peptide-spectrum matches (psm) in dictionary
  #   {key = fraction:scan; value = (peptide, score)}
  #   if more than 1 peptide, select the one with highest score
  psm_dict = {}
  # format of peaks.db.psm.csv
  col_peptide = 0
  col_score = 1
  col_fraction = 8
  col_scan = 9
  print("input_psm_file = {0:s}".format(input_psm_file))
  with open(input_psm_file, mode="r") as psm_handle:
    psm_line = psm_handle.readline() # header line
    for psm_line in psm_handle.readlines():
      psm_line_list = re.split(',|\n', psm_line)
      peptide = psm_line_list[col_peptide]
      score = float(psm_line_list[col_score])
      fraction = psm_line_list[col_fraction]
      scan = psm_line_list[col_scan]
      key = "F" + fraction + ":" + scan
      # PEAKS DB report duplicates I/L with same score >> score is useless here?
      if (not key in psm_dict) or (score > psm_dict[key][1]):
        psm_dict[key] = (peptide, score)
  print("psm_dict length = {0:d}".format(len(psm_dict)))

  ### read raw mgf and assign peptide/"" according to fraction:scan
  spectrum_count = 0
  spectrum_count_labeled = 0
  for fraction, input_file in zip(fraction_list, input_mgf_list):
    print("input_file = {0:s}".format(input_file))
    output_file = input_file + ".labeled"
    with open(output_file, mode="w") as output_handle:
      with open(input_file, mode="r") as input_handle:
          line = input_handle.readline()
          while line:
            if "SCANS=" in line:
              spectrum_count += 1
              scan = re.split('=|\n', line)[1]
              key = "F" + str(fraction) + ":" + scan
              peptide = ""
              if key in psm_dict:
                peptide = psm_dict[key][0]
                spectrum_count_labeled += 1
              # write scan line
              output_handle.write(line)
              # write rt line
              line = input_handle.readline()
              assert "RTINSECONDS=" in line, "Error: RTINSECONDS="
              output_handle.write(line)
              # write seq line
              seq_line = "SEQ=" + peptide + "\n"
              output_handle.write(seq_line)
            else:
              output_handle.write(line)
            line = input_handle.readline()
            # sanity check: SEQ= should not be in the raw mgf file
            assert not"SEQ=" in line, "Error: SEQ="
  print("spectrum_count = {0:d}".format(spectrum_count))
  print("spectrum_count_labeled = {0:d}".format(spectrum_count_labeled))


#~ folder_path = "data.training/yeast.high.seidel_2016/"
#~ fraction_list = range(1, 6+1)
#~ generate_labeled_mgf(fraction_list,
                     #~ input_mgf_list=[folder_path + str(x) + "_frac.raw.mgf"
                                     #~ for x in fraction_list],
                     #~ input_psm_file=folder_path + "peaks.db.psm.csv")






def split_dia(input_file):
  """TODO(nh2tran): docstring.
     Split a dia entry (spectrum, peptide list) into multiple entries of
     (spectrum, peptide).
     The spectrum is not changed, the new entries share the same spectrum.
  """

  print("".join(["="] * 80)) # section-separating line
  print("split_dia()")

  output_file = input_file + ".split"
  print("input_file = {0:s}".format(input_file))
  print("output_file = {0:s}".format(output_file))

  input_handle = open(input_file, mode="r")
  output_handle = open(output_file, mode="w")

  # find and count all "BEGIN IONS"
  location_list = []
  keyword = "BEGIN IONS"
  line = True
  while line:
    location = input_handle.tell()
    line = input_handle.readline()
    if keyword in line:
      location_list.append(location)
  spectrum_count_dia = len(location_list)

  # process each dia spectrum
  spectrum_count_split = 0
  spectrum_count_total = 0
  for location in location_list:

    input_handle.seek(location)

    # header BEGIN IONS
    line_begin = input_handle.readline()
    assert "BEGIN IONS" in line_begin, "Error: wrong input BEGIN IONS"

    # header TITLE
    line_title = input_handle.readline()
    assert "TITLE=" in line_title, "Error: wrong input TITLE"

    # header PEPMASS
    line_mass = input_handle.readline()
    assert "PEPMASS=" in line_mass, "Error: wrong input PEPMASS"
    precursor_mz_list = re.split('=|\n|\r', line_mass)[1]

    # header CHARGE
    line_charge = input_handle.readline()
    assert "CHARGE=" in line_charge, "Error: wrong input CHARGE"
    charge_list = re.split('=|\n|\r', line_charge)[1]

    # header SCANS
    line_scan = input_handle.readline()
    assert "SCANS=" in line_scan, "Error: wrong input SCANS"
    scan = re.split('=|\n|\r', line_scan)[1]

    # header RTINSECONDS
    line_rt = input_handle.readline()
    assert "RTINSECONDS=" in line_rt, "Error: wrong input RTINSECONDS"

    # header SEQ
    if deepnovo_config.FLAGS.header_seq:
      line_seq = input_handle.readline()
      assert "SEQ=" in line_seq, "Error: wrong input SEQ"
      raw_sequence_list = re.split('=|\n|\r', line_seq)[1]

    # ion list
    ion_list = []
    line_ion = input_handle.readline()
    while not "END IONS" in line_ion:
      ion_list.append(line_ion)
      line_ion = input_handle.readline()

    # split entries if not empty
    if precursor_mz_list:
      precursor_mz_list = re.split(';', precursor_mz_list)
      charge_list = re.split(';', charge_list)
      raw_sequence_list = re.split(';', raw_sequence_list)

      # sanity check
      assert len(precursor_mz_list) == len(charge_list)
      if deepnovo_config.FLAGS.header_seq:
        assert len(precursor_mz_list) == len(raw_sequence_list)

      entry_count = len(precursor_mz_list)
      spectrum_count_split += 1
      spectrum_count_total += entry_count

      # write each entry to output_file
      for entry in xrange(entry_count):
        output_list = []
        output_list.append(line_begin)
        output_list.append(line_title)
        output_list.append("PEPMASS=" + precursor_mz_list[entry] + "\n")
        output_list.append("CHARGE=" + charge_list[entry] + "\n")
        output_list.append("SCANS=" + scan + ":" + str(entry) + "\n")
        output_list.append(line_rt)
        if deepnovo_config.FLAGS.header_seq:
          output_list.append("SEQ=" + raw_sequence_list[entry] + "\n")
        output_list += ion_list
        output_list.append(line_ion) # "END IONS"
        for line in output_list:
          output_handle.write(line)
        output_handle.write("\n")

  input_handle.close()
  output_handle.close()
  print("spectrum_count_dia = {0:d}".format(spectrum_count_dia))
  print("spectrum_count_split = {0:d}".format(spectrum_count_split))
  print("spectrum_count_total = {0:d}".format(spectrum_count_total))

#~ split_dia(input_file="data.training/dia.umpire.dec15/fraction_1.mgf.train.adjacent")
#~ split_dia(input_file="data.training/dia.umpire.dec15/fraction_1.mgf.valid.adjacent")
#~ split_dia(input_file="data.training/dia.umpire.dec15/fraction_1.mgf.test.adjacent")




# Partition a dataset into three random or adjacent sets train-valid-test
#     with a distribution, e.g. 90-5-5 percent.
#
def partition_train_valid_test_dup_mgf(input_file, prob, mode):
  print("partition_train_valid_test_dup_mgf()")

  print("input_file = ", os.path.join(input_file))
  print("prob = ", prob)
  print("mode = ", mode)

  output_file_train = input_file + ".train." + mode
  output_file_valid = input_file + ".valid." + mode
  output_file_test = input_file + ".test." + mode

  # if mode is adjacent, the set size is pre-determined
  with open(input_file, mode="r") as input_handle:
    total_number = sum(1 for line in input_handle if "BEGIN IONS" in line)
    train_number = int(round(total_number * prob[0]))
    valid_number = int(round(total_number * prob[1]))
    test_number = total_number - (train_number + valid_number)

  with open(input_file, mode="r") as input_handle:
    with open(output_file_train, mode="w") as output_handle_train:
      with open(output_file_valid, mode="w") as output_handle_valid:
        with open(output_file_test, mode="w") as output_handle_test:
          counter = 0
          counter_train = 0
          counter_valid = 0
          counter_test = 0
          line = input_handle.readline()
          while line:
            if "BEGIN IONS" in line: # a spectrum found
              counter += 1
              if mode == "random":
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
              else:
                if counter_train < train_number:
                  output_handle = output_handle_train
                  counter_train += 1
                elif counter_valid < valid_number:
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
  
  print("counter ", counter)
  print("counter_train ", counter_train)
  print("counter_valid ", counter_valid)
  print("counter_test ", counter_test)

#~ partition_train_valid_test_dup_mgf(input_file="data.training/yeast.low.heinemann_2015/peaks.db.frac.mgf",
                                   #~ prob=[0.90, 0.05, 0.05],
                                   #~ mode="random")







# Partition a dataset into three random sets train-valid-test
#     with a distribution, e.g. 90-5-5 percent.
#
# This version removes all duplicated peptides so that 
#     each peptide has only one spectrum (selected randomly).
#
def partition_train_valid_test_unique_mgf(input_file, prob):
  print("partition_train_valid_test_unique_mgf()")

  print("input_file = ", os.path.join(input_file))
  print("prob = ", prob)
  
  output_file_train = input_file + ".train" + ".unique"
  output_file_valid = input_file + ".valid" + ".unique"
  output_file_test = input_file + ".test" + ".unique"
  
  peptide_list = []
  
  with open(input_file, mode="r") as input_handle:
    with open(output_file_train, mode="w") as output_handle_train:
      with open(output_file_valid, mode="w") as output_handle_valid:
        with open(output_file_test, mode="w") as output_handle_test:
          counter = 0
          counter_train = 0
          counter_valid = 0
          counter_test = 0
          line = input_handle.readline()
          while line:
            if "BEGIN IONS" in line: # a spectrum found

              line_buffer = []
              line_buffer.append(line)
              
              # TITLE
              line = input_handle.readline()
              line_buffer.append(line)
        
              # PEPMASS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # CHARGE
              line = input_handle.readline()
              line_buffer.append(line)
              
              # SCANS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # RTINSECONDS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # SEQ
              line = input_handle.readline()
              line_buffer.append(line)
              peptide = re.split('=|\n|\r', line)[1]
              
              if not (peptide in peptide_list): # new peptide

                peptide_list.append(peptide)
                
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
                  
                for l in line_buffer:
                  output_handle.write(l)

                while line and not ("END IONS" in line):
                  line = input_handle.readline()
                  output_handle.write(line)

                output_handle.write("\n")
                  
            line = input_handle.readline()

  print("counter ", counter)
  print("counter_train ", counter_train)
  print("counter_valid ", counter_valid)
  print("counter_test ", counter_test)

#~ partition_train_valid_test_unique_mgf(input_file="data/human.PXD002179.sds/peaks.db.mgf",
                                      #~ prob=[1.0, 0.0, 0.0])
                                   
                                   




# Partition a dataset into three random sets train-valid-test
#     with a distribution, e.g. 90-5-5 percent.
#
# This version removes duplicated peptides so that 
#     each peptide has at most max_spectra_per_peptide (selected randomly).
#
def partition_train_valid_test_unique_control_mgf(input_file, prob, max_spectra_per_peptide):
  print("partition_train_valid_test_unique_control_mgf()")

  print("input_file = ", os.path.join(input_file))
  print("prob = ", prob)
  
  output_file_train = input_file + ".train" + ".unique" + str(max_spectra_per_peptide)
  output_file_valid = input_file + ".valid" + ".unique" + str(max_spectra_per_peptide)
  output_file_test = input_file + ".test" + ".unique" + str(max_spectra_per_peptide)
  
  peptide_list = []
  peptide_spectra_count = {}
  
  with open(input_file, mode="r") as input_handle:
    with open(output_file_train, mode="w") as output_handle_train:
      with open(output_file_valid, mode="w") as output_handle_valid:
        with open(output_file_test, mode="w") as output_handle_test:
          counter = 0
          counter_train = 0
          counter_valid = 0
          counter_test = 0
          line = input_handle.readline()
          while line:
            if "BEGIN IONS" in line: # a spectrum found

              line_buffer = []
              line_buffer.append(line)
              
              # TITLE
              line = input_handle.readline()
              line_buffer.append(line)
        
              # PEPMASS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # CHARGE
              line = input_handle.readline()
              line_buffer.append(line)
              
              # SCANS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # RTINSECONDS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # SEQ
              line = input_handle.readline()
              line_buffer.append(line)
              peptide = re.split('=|\n|\r', line)[1]
              
              if not (peptide in peptide_list): # new peptide
                peptide_list.append(peptide)
                peptide_spectra_count[peptide] = 0
                
              if (peptide_spectra_count[peptide] < max_spectra_per_peptide):
                
                peptide_spectra_count[peptide] += 1
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
                  
                for l in line_buffer:
                  output_handle.write(l)

                while line and not ("END IONS" in line):
                  line = input_handle.readline()
                  output_handle.write(line)

                output_handle.write("\n")
                  
            line = input_handle.readline()

  print("counter ", counter)
  print("counter_train ", counter_train)
  print("counter_valid ", counter_valid)
  print("counter_test ", counter_test)

#~ partition_train_valid_test_unique_control_mgf(input_file="data/human.cancer/peaks.db.frac_1_10.mgf",
                                              #~ prob=[1.0, 0.0, 0.0],
                                              #~ max_spectra_per_peptide=4)






# Partition a dataset into three random sets train-valid-test
#     with a distribution, e.g. 90-5-5 percent
#
# Each peptide may correspond to multiple different spectra,
#     but the three sets do not share common peptides.
#
def partition_train_valid_test_repeat_mgf(input_file,prob):
  print("partition_train_valid_test_repeat_mgf()")

  print("input_file = ", os.path.join(input_file))
  print("prob = ", prob)
  
  output_file_train = input_file + ".train" + ".repeat"
  output_file_valid = input_file + ".valid" + ".repeat"
  output_file_test = input_file + ".test" + ".repeat"
  
  peptide_train_list = []
  peptide_valid_list = []
  peptide_test_list = []
  
  with open(input_file, mode="r") as input_handle:
    with open(output_file_train, mode="w") as output_handle_train:
      with open(output_file_valid, mode="w") as output_handle_valid:
        with open(output_file_test, mode="w") as output_handle_test:
          counter = 0
          counter_train = 0
          counter_valid = 0
          counter_test = 0
          counter_unique = 0
          line = input_handle.readline()
          while line:
            if "BEGIN IONS" in line: # a spectrum found

              line_buffer = []
              line_buffer.append(line)
              
              # TITLE
              line = input_handle.readline()
              line_buffer.append(line)
        
              # PEPMASS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # CHARGE
              line = input_handle.readline()
              line_buffer.append(line)
              
              # SCANS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # RTINSECONDS
              line = input_handle.readline()
              line_buffer.append(line)
        
              # SEQ
              line = input_handle.readline()
              line_buffer.append(line)
              peptide = re.split('=|\n|\r', line)[1]

              # found a spectrum and a peptide
              counter += 1

              # check if the peptide already exists in any of the three lists
              # if yes, this new spectrum will be assigned to that list
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
                  
              for l in line_buffer:
                output_handle.write(l)

              while line and not ("END IONS" in line):
                line = input_handle.readline()
                output_handle.write(line)

              output_handle.write("\n")
                  
            line = input_handle.readline()

  print("counter ", counter)
  print("counter_train ", counter_train)
  print("counter_valid ", counter_valid)
  print("counter_test ", counter_test)
  print("counter_unique ", counter_unique)

#~ partition_train_valid_test_repeat_mgf(input_file="data.training/yeast.low.heinemann_2015/peaks.db.mgf",
                                      #~ prob=[0.90, 0.05, 0.05])







# filter spectra with 
#     MAX_LEN = 20 aa
#     MZ_MAX = 3000 Da 
#     PRECURSOR_MASS_PRECISION_INPUT_FILTER = 0.01 Da
#
# re-number scan id to unique because the input file was merged from 
#     different fractions
#
#
# and extract ground-truth peptide sequences from database-search
#
def prepare_test_file(input_file):
  print("prepare_test_file()")
  
  print("input_file = ", os.path.join(input_file))
  #
  #~ output_file = input_file + ".filter"
  #~ print("output_file = ", output_file)
  #
  dbseq_file = input_file + ".dbseq"
  print("dbseq_file = ", dbseq_file)

  counter = 0
  counter_skipped = 0
  counter_skipped_mod = 0
  counter_skipped_len = 0
  counter_skipped_mass = 0
  
  with open(input_file, mode="r") as input_handle:
    #~ with open(output_file, mode="w") as output_handle:
      with open(dbseq_file, mode="w") as dbseq_handle:
        
        print("scan \t target_seq \n", file=dbseq_handle, end="")

        line = input_handle.readline()
        while line:
  
          if "BEGIN IONS" in line: # a spectrum found
  
            line_buffer = []
            line_buffer.append(line)
        
            unknown_modification = False
            max_intensity = 0.0
  
            # header TITLE
            line = input_handle.readline()
            line_buffer.append(line)
  
            # header PEPMASS
            line = input_handle.readline()
            peptide_ion_mz = float(re.split('=|\n', line)[1])
            line_buffer.append(line)
  
            # header CHARGE
            line = input_handle.readline()
            charge = float(re.split('=|\+', line)[1])
            line_buffer.append(line)
        
            # header SCANS
            line = input_handle.readline()
            #~ scan = int(re.split('=', line)[1])
            scan = re.split('=|\n', line)[1]
            line_buffer.append(line)
  
            # header RTINSECONDS
            line = input_handle.readline()
            line_buffer.append(line)
  
            # header SEQ
            line = input_handle.readline()
            line_buffer.append(line)
            raw_sequence = re.split('=|\n|\r', line)[1]
            #~ peptide = peptide.translate(None, '+-.0123456789)') # modification signal "("
            #~ peptide = peptide.translate(None, '(+-.0123456789)') # ignore modifications
            raw_sequence_len = len(raw_sequence)
            peptide = []
            index = 0
            while (index < raw_sequence_len):
              if (raw_sequence[index]=="("):
                if (peptide[-1]=="C" and raw_sequence[index:index+8]=="(+57.02)"):
                #~ if (peptide[-1]=="C" and raw_sequence[index:index+8]=="(+58.01)"):
                  peptide[-1] = "C(Carbamidomethylation)"
                  index += 8
                elif (peptide[-1]=='M' and raw_sequence[index:index+8]=="(+15.99)"):
                  peptide[-1] = 'M(Oxidation)'
                  index += 8
                elif (peptide[-1]=='N' and raw_sequence[index:index+6]=="(+.98)"):
                  peptide[-1] = 'N(Deamidation)'
                  index += 6
                elif (peptide[-1]=='Q' and raw_sequence[index:index+6]=="(+.98)"):
                  peptide[-1] = 'Q(Deamidation)'
                  index += 6
                else: # unknown modification
                #~ elif ("".join(raw_sequence[index:index+8])=="(+42.01)"):
                  #~ print("ERROR: unknown modification!")
                  #~ print("raw_sequence = ", raw_sequence)
                  #~ sys.exit()
                  unknown_modification = True
                  break
              else:
                peptide.append(raw_sequence[index])
                index += 1
            #
            if (unknown_modification):
              counter_skipped += 1
              counter_skipped_mod += 1
              continue
              
            # neutral peptide_mass <= MZ_MAX(3000.0) # TEMP
            peptide_mass = peptide_ion_mz*charge - charge*deepnovo_config.mass_H
            if (peptide_mass > deepnovo_config.MZ_MAX):
              counter_skipped += 1
              counter_skipped_mass += 1
              continue
  
            # TRAINING-ONLY: consider peptide length <= MAX_LEN(20)
            peptide_len = len(peptide)
            if (peptide_len > deepnovo_config.MAX_LEN): 
              print("ERROR: peptide_len {0} exceeds MAX_LEN {1}".format(peptide_len, deepnovo_config.MAX_LEN))
              sys.exit()
              #~ counter_skipped += 1
              #~ counter_skipped_len += 1
              #~ continue
      
            # TRAINING-ONLY: testing peptide_mass & sequence_mass
            #~ sequence_mass = sum(deepnovo_config.mass_AA[aa] for aa in peptide)
            #~ sequence_mass += deepnovo_config.mass_N_terminus + deepnovo_config.mass_C_terminus
            #~ if (abs(peptide_mass-sequence_mass) > deepnovo_config.PRECURSOR_MASS_PRECISION_INPUT_FILTER):
              #
              #~ print("ERROR: peptide_mass and sequence_mass not matched")
              #~ print("peptide = ", peptide)
              #~ print("peptide_list_mod = ", peptide_list_mod)
              #~ print("peptide_list = ", peptide_list)
              #~ print("peptide_ion_mz = ",peptide_ion_mz)
              #~ print("charge = ", charge)
              #~ print("peptide_mass  ", peptide_mass)
              #~ print("sequence_mass ", sequence_mass)
              #~ sys.exit()
              #
              #~ counter_skipped += 1
              #~ continue
  
            # AN ENTRY FOUND
            counter += 1
            if counter % 10000 == 0:
              print("  reading peptide %d" % counter)
              
            # output ground-truth peptide sequence & re-number scan id
            print("%s\t%s\n" % (scan, ",".join(peptide)),
                  file=dbseq_handle,
                  end="")
      
            # output this data entry
            #~ for l in line_buffer:
              #~ output_handle.write(l)
            #
            while line and not ("END IONS" in line):
              line = input_handle.readline()
              #~ output_handle.write(line)
            #
            #~ output_handle.write("\n")
                
          line = input_handle.readline()
  
  print("  total peptide read %d" % counter)
  print("  total peptide skipped %d" % counter_skipped)
  print("  total peptide skipped by mod %d" % counter_skipped_mod)
  print("  total peptide skipped by len %d" % counter_skipped_len)
  print("  total peptide skipped by mass %d" % counter_skipped_mass)
  
#~ prepare_test_file(input_file="data.training/yeast.low.heinemann_2015/peaks.db.mgf")






# Partition a dbseq file into 2 sets: overlapping & nonoverlapping with the trainseq file
def partition_dbseq(dbseq_file, trainseq_file):
  
  print("partition_dbseq()")
  print("dbseq_file = ", dbseq_file)
  print("trainseq_file = ", trainseq_file)

  trainseq = []

  with open(trainseq_file, mode="r") as trainseq_handle:

    # header
    trainseq_handle.readline()

    for line in trainseq_handle:
      line_split = re.split('\t|\n', line)
      #~ scan = line_split[0]
      peptide = line_split[1]
      trainseq.append(peptide)
      
  overlap_file = dbseq_file + ".overlap"
  nonoverlap_file = dbseq_file + ".nonoverlap"
  count = 0
  count_overlap = 0
  count_nonoverlap = 0

  with open(dbseq_file, mode="r") as dbseq_handle:
    with open(overlap_file, mode="w") as overlap_handle:
      with open(nonoverlap_file, mode="w") as nonoverlap_handle:

        # header
        line = dbseq_handle.readline()
        overlap_handle.write(line)
        nonoverlap_handle.write(line)

        for line in dbseq_handle:
          line_split = re.split('\t|\n', line)
          #~ scan = line_split[0]
          peptide = line_split[1]
          
          if (peptide in trainseq):
            overlap_handle.write(line)
            count_overlap += 1
          else:
            nonoverlap_handle.write(line)
            count_nonoverlap += 1
          
          count += 1
  
  print("count = {0:d}".format(count))
  print("count_overlap = {0:d}".format(count_overlap))
  print("count_nonoverlap = {0:d}".format(count_nonoverlap))

#~ partition_dbseq(dbseq_file="data/human.cancer/peaks.db.frac_21_41.mgf.dbseq",
                #~ trainseq_file="data/human.cancer/peaks.db.frac_1_20.mgf.dbseq")
            
            

def read_dbseq(dbseq_file):
  
  print("read_dbseq()")
  print("dbseq_file = ", dbseq_file)

  dbseq = {}

  batch_len_AA = 0.0

  with open(dbseq_file, mode="r") as dbseq_handle:

    # header
    dbseq_handle.readline()

    for line in dbseq_handle:
      line_split = re.split('\t|\n', line)
      scan = line_split[0]
      peptide = re.split(',', line_split[1])
      #
      dbseq[scan] = [deepnovo_config.vocab[x] for x in peptide]
      batch_len_AA += len(peptide)
  
  batch_size = len(dbseq)

  print("batch_size = ", batch_size)
  print("batch_len_AA = ", batch_len_AA)
  
  return dbseq, batch_size, batch_len_AA


def read_novonet(novonet_file):

  print("read_novonet()")
  print("novonet_file = ", novonet_file)

  novonet = {}

  with open(novonet_file, mode="r") as novonet_handle:
    
    # header
    novonet_handle.readline()
    
    for line in novonet_handle:
      
      line_split = re.split('\t|\n', line)
      scan = line_split[0]
      if (line_split[2] == ""): # empty output
        novonet_seq_id = []
      else:
        novonet_seq = re.split(',', line_split[2])
        novonet_seq_id = [deepnovo_config.vocab[x] for x in novonet_seq]
      
      novonet[scan] = novonet_seq_id
      
  return novonet
  

def read_peaks(peaks_denovo_file, peaks_format, alc_threshold):

  print("read_peaks()")
  print("peaks_denovo_file = ", peaks_denovo_file)
  
  if (peaks_format == "old_7.5"):
    peptide_column = 1
    alc_column = 3
  elif (peaks_format == "new_8.0"):
    peptide_column = 3
    alc_column = 5
  else:
    print("ERROR: wrong PEAKS denovo format")
    sys.exit()

  peaks = {}
  peaks_raw = {}

  with open(peaks_denovo_file, mode="r") as peaks_handle:
    
    # header
    peaks_handle.readline()
    
    for line in peaks_handle:
      
      line_split = re.split(",", line)
      
      if (peaks_format == "old_7.5"):
        scan = line_split[0]
      elif (peaks_format == "new_8.0"):
        scan = "F" + line_split[0] + ":" + line_split[1]
        
      if (line_split[peptide_column] == ""): # empty output
        peaks_seq_id = []
      else:
        raw_sequence = line_split[peptide_column]
        raw_sequence_len = len(raw_sequence)
        peptide = []
        index = 0
        while (index < raw_sequence_len):
          if (raw_sequence[index]=="("):
            if (peptide[-1]=="C" and raw_sequence[index:index+8]=="(+57.02)"):
            #~ if (peptide[-1]=="C" and raw_sequence[index:index+8]=="(+58.01)"):
              peptide[-1] = "C(Carbamidomethylation)"
              index += 8
            elif (peptide[-1]=='M' and raw_sequence[index:index+8]=="(+15.99)"):
              peptide[-1] = 'M(Oxidation)'
              index += 8
            elif (peptide[-1]=='N' and raw_sequence[index:index+6]=="(+.98)"):
              peptide[-1] = 'N(Deamidation)'
              index += 6
            elif (peptide[-1]=='Q' and raw_sequence[index:index+6]=="(+.98)"):
              peptide[-1] = 'Q(Deamidation)'
              index += 6
            else: # unknown modification
            #~ elif ("".join(raw_sequence[index:index+8])=="(+42.01)"):
              #~ print("ERROR: unknown modification!")
              #~ print("raw_sequence = ", raw_sequence)
              #~ sys.exit()
              unknown_modification = True
              break
          else:
            peptide.append(raw_sequence[index])
            index += 1
        #
        peaks_seq_id = [deepnovo_config.vocab[x] for x in peptide]

      alc_score = float(line_split[alc_column])
      if (alc_score >= alc_threshold):
        peaks[scan] = peaks_seq_id
        peaks_raw[scan] = raw_sequence
      
  return peaks, peaks_raw
  

def get_peaks_denovo_spectra(output_spectra_file, raw_spectra_file, peaks_denovo_file, peaks_format, alc_threshold=0):
  
  print("get_peaks_denovo_spectra()")
  print("peaks_denovo_file = ", peaks_denovo_file)
  print("ALC cut-off = ", alc_threshold)
  print("raw_spectra_file = ", raw_spectra_file)
  
  _, peaks_denovo_peptides = read_peaks(peaks_denovo_file, peaks_format, alc_threshold)
  print("peaks_denovo_peptides: ", len(peaks_denovo_peptides))
  
  counter_spectra = 0
  
  with open(raw_spectra_file, mode="r") as input_handle:
    with open(output_spectra_file, mode="w") as output_handle:

      line = input_handle.readline()
      while line:

        if "BEGIN IONS" in line: # a spectrum found

          line_buffer = []
          line_buffer.append(line)
      
          # header TITLE
          line = input_handle.readline()
          line_buffer.append(line)

          # header PEPMASS
          line = input_handle.readline()
          line_buffer.append(line)

          # header CHARGE
          line = input_handle.readline()
          line_buffer.append(line)
      
          # header SCANS
          line = input_handle.readline()
          #~ scan = int(re.split('=', line)[1])
          scan = re.split('=|\n', line)[1]
          line_buffer.append(line)

          # lookup scan id
          if not (scan in peaks_denovo_peptides):
            continue
          else:
            counter_spectra += 1
            #
            for l in line_buffer:
              output_handle.write(l)
            #
            # RTINSECONDS
            line = input_handle.readline()
            output_handle.write(line)
            #
            # SEQ
            line = "SEQ=" + peaks_denovo_peptides[scan] + "\n"
            output_handle.write(line)
            #
            while line and not ("END IONS" in line):
              line = input_handle.readline()
              output_handle.write(line)
            #
            output_handle.write("\n")
              
        line = input_handle.readline()

  print("total spectra found %d" % counter_spectra)
  
  
def test_AA_match_novor(decoder_input, output):

  num_match = 0
  
  decoder_input_len = len(decoder_input)
  output_len = len(output)

  decoder_input_mass = [deepnovo_config.mass_ID[x] for x in decoder_input]
  decoder_input_mass_cum = np.cumsum(decoder_input_mass)

  output_mass = [deepnovo_config.mass_ID[x] for x in output]
  output_mass_cum = np.cumsum(output_mass)
  
  i = 0
  j = 0
  while (i < decoder_input_len and j < output_len):

    #~ # for testing
    #~ print(decoder_input_mass_cum[i])
    #~ print(output_mass_cum[j])

    if (abs(decoder_input_mass_cum[i] - output_mass_cum[j]) < 0.5):

      #~ # for testing
      #~ print(decoder_input_mass[i] )
      #~ print(output_mass[j])

      if (abs(decoder_input_mass[i] - output_mass[j]) < 0.1):
      #~ if  decoder_input[index_aa]==output[index_aa]:
        num_match += 1

      i += 1
      j += 1
    elif (decoder_input_mass_cum[i] < output_mass_cum[j]):
      i += 1
    else:
      j += 1

    #~ # for testing
    #~ print(num_match)
    
  return num_match


def test_accuracy(dbseq_file, denovo_file, tool, peaks_format=None, alc_threshold=None):
  
  print("test_accuracy()")
  
  batch_accuracy_AA = 0.0
  batch_len_decode = 0.0
  num_exact_match = 0.0
  num_len_match = 0.0

  dbseq, batch_size, batch_len_AA = read_dbseq(dbseq_file)
  
  if (tool == "novonet"):
    denovo = read_novonet(denovo_file)
  elif (tool == "peaks"):
    denovo, _ = read_peaks(denovo_file, peaks_format, alc_threshold)
    
  count_skipped = 0

  # for testing
  test_output = dict.fromkeys(dbseq.keys(),[])

  for scan, seq in denovo.iteritems():
    
    if (scan in dbseq):

      accuracy_AA = test_AA_match_novor(dbseq[scan], seq)
      len_AA = len(dbseq[scan])
      
      # for testing
      output_seq = [deepnovo_config.vocab_reverse[x] for x in seq]
      test_output[scan] = [output_seq, accuracy_AA]
      
      len_decode = len(seq) 
      batch_len_decode += len_decode
      
      batch_accuracy_AA += accuracy_AA
      #~ batch_accuracy_AA += accuracy_AA/len_AA
      #
      if (accuracy_AA==len_AA):
        num_exact_match += 1.0
      #
      if (len(seq)==len_AA):
        num_len_match += 1.0
    
    else:
      
      count_skipped += 1
    
  # for testing
  with open("test_accuracy.tab", "w") as file_handle:
    file_handle.write("scan \t target_seq \t target_len \t output_seq \t accuracy_AA \n")
    for scan, output in test_output.iteritems():
      target_seq = [deepnovo_config.vocab_reverse[x] for x in dbseq[scan]]
      target_len = len(target_seq)
      if (not output):
        file_handle.write("{0:s}\t{1:s}\t{2:d}\t{3:s}\t{4:d}\n".format(scan, target_seq, target_len, [], 0))
      else:
        file_handle.write("{0:s}\t{1:s}\t{2:d}\t{3:s}\t{4:d}\n".format(scan, target_seq, target_len, output[0], output[1]))

  print("  recall_AA %.4f" % (batch_accuracy_AA/batch_len_AA))
  #~ print("  accuracy_AA %.4f" % (batch_accuracy_AA/batch_size))
  print("  precision_AA %.4f" % (batch_accuracy_AA/batch_len_decode))
  #
  print("  recall_peptide %.4f" % (num_exact_match/batch_size))
  print("  recall_len %.4f" % (num_len_match/batch_size))
  print("  count_skipped (not in dbseq) %d" % (count_skipped))
  

# NovoNet
#~ test_accuracy(dbseq_file="data.training/yeast.low.coon_2013/peaks.db.mgf.test.dup.dbseq",
              #~ denovo_file = "train.example/decode_output.tab",
              #~ tool="novonet")

# PEAKS
#~ test_accuracy(dbseq_file="data.training/yeast.high.seidel_2016/peaks.db.mgf.dbseq",
              #~ denovo_file = "data.training/yeast.high.seidel_2016/peaks.denovo.csv",
              #~ tool="peaks",
              #~ #peaks_format="old_7.5",
              #~ peaks_format="new_8.0",
              #~ alc_threshold=0)





