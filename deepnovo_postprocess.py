# Copyright 2017 Hieu Tran. All Rights Reserved.
#
# DeepNovo is publicly available for non-commercial uses.
# ==============================================================================

"""TODO(nh2tran): docstring."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import random
import sys
import re

from Bio import SeqIO
from Bio.SeqIO import FastaIO
import Levenshtein

import csv
import numpy as np
import math
import deepnovo_config


def compute_peptide_mass(peptide):
  """TODO(nh2tran): docstring.
  """

  #~ print("".join(["="] * 80)) # section-separating line ===
  #~ print("WorkerDB: _compute_peptide_mass()")

  peptide_mass = (deepnovo_config.mass_N_terminus
                  + sum(deepnovo_config.mass_AA[aa] for aa in peptide)
                  + deepnovo_config.mass_C_terminus)

  return peptide_mass

# ~ peptide = 'TASSQRLR'
# ~ print(peptide)
# ~ print(compute_peptide_mass(peptide))


def read_feature_accuracy(input_file):

  feature_list = []
  with open(input_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter='\t')
    for row in csv_reader:
      feature = {}
      feature['feature_id'] = row['feature_id']
      feature['feature_area'] = math.log10(max(float(row['feature_area']), 1.0))
      feature['predicted_score'] = float(row['predicted_score'])
      feature['recall_AA'] = float(row['recall_AA'])
      feature['predicted_len'] = float(row['predicted_len'])
      feature_list.append(feature)
  return feature_list


def find_score_cutoff(accuracy_file, accuracy_cutoff):
  """TODO(nh2tran): docstring."""

  print("".join(["="] * 80)) # section-separating line
  print("find_score_cutoff()")

  print('accuracy_file =', accuracy_file)
  print('accuracy_cutoff =', accuracy_cutoff)

  feature_list = read_feature_accuracy(accuracy_file)
  feature_list_sorted = sorted(feature_list, key=lambda k: k['predicted_score'], reverse=True)
  recall_cumsum = np.cumsum([f['recall_AA'] for f in feature_list_sorted])
  predicted_len_cumsum = np.cumsum([f['predicted_len'] for f in feature_list_sorted])
  accuracy_cumsum = recall_cumsum / predicted_len_cumsum
  #cutoff_index = np.flatnonzero(accuracy_cumsum < accuracy_cutoff)[0]
  cutoff_index = np.flatnonzero(accuracy_cumsum >= accuracy_cutoff)[-1]
  cutoff_score = feature_list_sorted[cutoff_index]['predicted_score']
  print('cutoff_index = ', cutoff_index)
  print('cutoff_score = ', cutoff_score)
  print('cutoff_score = ', 100*math.exp(cutoff_score))

  return cutoff_score


def select_top_score(input_file, output_file, score_cutoff):
  """Select a threshold of de novo confidence scores to filter de novo results.
     The score threshold is calculated based on a 95% cutoff of the testing accuracy.

     Usage:
       accuracy_cutoff = 0.95
       accuracy_file = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo.accuracy"
       score_cutoff = find_score_cutoff(accuracy_file, accuracy_cutoff)
       input_file = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/feature.csv.mass_corrected.deepnovo_denovo"
       output_file = input_file + ".top95"
       select_top_score(input_file, output_file, score_cutoff)
  """

  print("".join(["="] * 80)) # section-separating line
  print("select_top_score()")

  print('input_file = ', input_file)
  print('output_file = ', output_file)
  print('score_cutoff = ', score_cutoff)

  total_feature = 0
  select_feature = 0
  with open(input_file, 'r') as input_handle:
    with open(output_file, 'w') as output_handle:
      csv_reader = csv.DictReader(input_handle, delimiter='\t')
      csv_writer = csv.DictWriter(output_handle, csv_reader.fieldnames, delimiter='\t')
      csv_writer.writeheader()
      for row in csv_reader:
        total_feature += 1
        predicted_score = float(row['predicted_score']) if row['predicted_score'] else -999
        if predicted_score >= score_cutoff:
          select_feature += 1
          csv_writer.writerow(row)
  print('total_feature = ', total_feature)
  print('select_feature = ', select_feature)
          

def convert_I_to_L(input_file, output_file):
  """Convert I (Isoleucine) to L (Leucine).

     Usage:
       input_file = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/feature.csv.mass_corrected.deepnovo_denovo.top95"
       output_file = input_file + ".I_to_L"
       convert_I_to_L(input_file, output_file)
  """

  print("".join(["="] * 80)) # section-separating line
  print("convert_I_to_L()")

  print('input_file = ', input_file)
  print('output_file = ', output_file)

  with open(input_file, 'r') as input_handle:
    with open(output_file, 'w') as output_handle:
      csv_reader = csv.DictReader(input_handle, delimiter='\t')
      csv_reader.fieldnames.append('before_I_to_L')
      csv_writer = csv.DictWriter(output_handle, csv_reader.fieldnames, delimiter='\t')
      csv_writer.writeheader()
      for row in csv_reader:
        predicted_sequence = row['predicted_sequence']
        row['before_I_to_L'] = predicted_sequence
        row['predicted_sequence'] = predicted_sequence.replace('I', 'L')
        csv_writer.writerow(row)
          

def compute_distance(predicted_sequence, consensus_sequence):
  """TODO(nh2tran): docstring.
  """

  #~ print("".join(["="] * 80)) # section-separating line ===
  #~ print("compute_distance()")

  # simplify the modifications
  modification_list = ['C(Carbamidomethylation)', 'M(Oxidation)', 'N(Deamidation)', 'Q(Deamidation)']
  simplified_list = ['c', 'm', 'n', 'q']
  for x in simplified_list:
    assert x not in deepnovo_config.vocab_reverse
  for x, y in zip(modification_list, simplified_list):
    predicted_sequence = [aa.replace(x, y) for aa in predicted_sequence]
    consensus_sequence = [aa.replace(x, y) for aa in consensus_sequence]
  predicted_sequence = ''.join(predicted_sequence)
  consensus_sequence = ''.join(consensus_sequence)

  return Levenshtein.distance(predicted_sequence, consensus_sequence)


def correct_by_consensus(input_file, output_file):
  """Correct de novo sequencing errors as following:
       group predicted sequences of the same mass together;
       vote the consensus sequence;
       replace the predicted by the consensus to correct errors: AB-BA, Q-AG, N-GG, etc.
  
     Usage:
       input_file = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L"
       output_file = input_file + ".consensus"
       correct_by_consensus(input_file, output_file)
  """

  print("".join(["="] * 80)) # section-separating line
  print("correct_by_consensus()")

  print('input_file = ', input_file)
  print('output_file = ', output_file)

  total_feature = 0
  empty_feature = 0
  mass_dict = {}
  with open(input_file, 'r') as input_handle:
    with open(output_file, 'w') as output_handle:
      csv_reader = csv.DictReader(input_handle, delimiter='\t')
      csv_reader.fieldnames.append('before_consensus')
      csv_writer = csv.DictWriter(output_handle, csv_reader.fieldnames, delimiter='\t')
      csv_writer.writeheader()

      # build the sequence mass dictionary
      # all sequences with the same mass are grouped together
      # (same mass up to resolution 1e4)
      for row in csv_reader:
        total_feature += 1
        predicted_sequence = row['predicted_sequence']
        # skip empty sequences that DeepNovo couldn't find a suitable candidate with the given mass
        if predicted_sequence == '':
          empty_feature += 1
          continue
        # save the original predicted sequence before correcting it later
        row['before_consensus'] = predicted_sequence

        predicted_sequence = predicted_sequence.split(',')
        predicted_score = float(row['predicted_score'])
        sequence_mass_index = int(round(compute_peptide_mass(predicted_sequence)
                                        * deepnovo_config.KNAPSACK_AA_RESOLUTION))
        feature = {'row': row,
                   'predicted_sequence': predicted_sequence,
                   'predicted_score': predicted_score}
        if sequence_mass_index in mass_dict:
          mass_dict[sequence_mass_index].append(feature)
        else:
          mass_dict[sequence_mass_index] = [feature]
      # check if all sequences have been assigned
      assigned_feature = sum([len(x) for x in mass_dict.values()])
      assert total_feature - empty_feature == assigned_feature

      # for each group of sequences of the same mass,
      # vote the consensus sequence;
      # calculate Levenshtein distance between each sequence and the consensus;
      # if 1 <= distance <= 2, replace the sequence by the consensus;
      # (distance = 2 examples: AB-BA, Q-AG, N-GG)
      # write to output.
      for group in mass_dict.values():
        if len(group) == 1:
          consensus_sequence = group[0]['predicted_sequence']
        else:
          # vote the consensus sequence
          # the easy way is to find the sequence with the highest score and frequency
          # (more complicated ways: De Bruijn graph, alignment)
          consensus_candidate = {}
          for feature in group:
            predicted_sequence = feature['predicted_sequence']
            predicted_score_prob = 100*math.exp(feature['predicted_score'])
            predicted_sequence = ','.join(predicted_sequence)
            if predicted_sequence in consensus_candidate:
              consensus_candidate[predicted_sequence] += predicted_score_prob
            else:
              consensus_candidate[predicted_sequence] = predicted_score_prob
          consensus_sequence = max(consensus_candidate.iterkeys(), key=(lambda key: consensus_candidate[key]))
          consensus_sequence = consensus_sequence.split(',')

        # calculate distance, correct sequence by the consensus, write to output
        for feature in group:
          distance = compute_distance(feature['predicted_sequence'], consensus_sequence)
          if 1 <= distance <= 2:
            feature['row']['predicted_sequence'] = ','.join(consensus_sequence)
          csv_writer.writerow(feature['row'])

      print('total_feature = ', total_feature)
      print('empty_feature = ', empty_feature)
      print('assigned_feature = ', assigned_feature)
          

def filter_by_minlen(input_file, output_file, minlen):
  """Filter out sequences of length less than minlen.
  
     Usage:
       minlen = 5
       input_file = "data.training/aa.hla.bassani.nature_2016.mel_16.class_1/feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus"
       output_file = input_file + ".minlen" + str(minlen)
       filter_by_minlen(input_file, output_file, minlen)
  """

  print("".join(["="] * 80)) # section-separating line
  print("filter_by_minlen()")
  print('input_file = ', input_file)
  print('output_file = ', output_file)
  print('minlen = ', minlen)

  total_feature = 0
  minlen_feature = 0
  removed_feature = 0
  with open(input_file, 'r') as input_handle:
    with open(output_file, 'w') as output_handle:
      csv_reader = csv.DictReader(input_handle, delimiter='\t')
      csv_writer = csv.DictWriter(output_handle, csv_reader.fieldnames, delimiter='\t')
      csv_writer.writeheader()
      for row in csv_reader:
        total_feature += 1
        predicted_sequence_len = len(re.split(',', row['predicted_sequence']))
        if predicted_sequence_len >= minlen:
          csv_writer.writerow(row)
          minlen_feature += 1
        else:
          removed_feature += 1
  print('total_feature = ', total_feature)
  print('minlen_feature = ', minlen_feature)
  print('removed_feature = ', removed_feature)
          

def database_lookup(input_fasta_file, input_denovo_file, output_file, split_char, col_sequence):

  print("".join(["="] * 80)) # section-separating line
  print("database_lookup()")

  print('input_fasta_file = ', input_fasta_file)
  print('input_denovo_file = ', input_denovo_file)
  print('output_file = ', output_file)

  with open(input_fasta_file, 'r') as input_fasta_handle:
    record_list = list(SeqIO.parse(input_fasta_handle, "fasta"))
    print("Number of protein sequences: ", len(record_list))

  total_count = 0 
  db_count = 0
  denovo_count = 0
  with open(input_denovo_file, 'r') as input_denovo_handle:
    with open(output_file, 'w') as output_handle:
      # header
      header_line = input_denovo_handle.readline()
      print(header_line, file=output_handle, end="")
      for line in input_denovo_handle:
        total_count += 1
        line_split = re.split(split_char, line)
        line_split = line_split[:-1] # exclude the last empty ""
        predicted_sequence = line_split[col_sequence]
        predicted_sequence = predicted_sequence.replace(',', '')
        predicted_sequence = predicted_sequence.replace('C(Carbamidomethylation)', 'C')
        indb = False
        for record in record_list:
          if predicted_sequence in record.seq:
            indb = True
            break
        if indb:
          db_count += 1
          line_split.append("db")
        else:
          denovo_count += 1
          line_split.append("denovo")
        print('\t'.join(line_split), file=output_handle, end="\n")
  print('total_count = ', total_count)
  print('db_count = ', db_count)
  print('denovo_count = ', denovo_count)

# ~ input_fasta_file = "data/uniprot_sprot.human.fasta"
# ~ input_denovo_file = "data.training/bassani.hla.2018_10_18.correct_mass_shift/unidentified_features.csv.deepnovo_denovo.top95"
# ~ output_file = input_denovo_file + ".lookup"
# ~ split_char = '\t|\n'
# ~ col_sequence = 2
# ~ database_lookup(input_fasta_file, input_denovo_file, output_file, split_char, col_sequence)


def select_top_k(input_file, output_file, top_k, split_char, col_score):
  """TODO(nh2tran): docstring."""

  print("".join(["="] * 80)) # section-separating line
  print("select_top_k()")

  print('input_file = ', input_file)
  print('output_file = ', output_file)
  print('top_k = ', top_k)

  with open(input_file, 'r') as input_handle:
    with open(output_file, 'w') as output_handle:
      # header
      header_line = input_handle.readline()
      print(header_line, file=output_handle, end="")
      predicted_list = []
      for line in input_handle:
        line_split = re.split(split_char, line)
        predicted = {}
        predicted["line"] = line
        predicted["score"] = float(line_split[col_score]) if line_split[col_score] else -999
        predicted_list.append(predicted)
      sorted_list = sorted(predicted_list, key=lambda k: k['score'], reverse=True) 
      for entry in sorted_list[:top_k]:
        print(entry["line"], file=output_handle, end="")
          
#~ top_k = 7673
#~ split_char = '\t|\n'
#~ col_score = deepnovo_config.pcol_score_max
#~ input_file = "data.training/dia.pecan.plasma.2018_03_29/testing.unlabeled.csv.deepnovo_denovo"
#~ output_file = input_file + ".topk"
#~ select_top_k(input_file, output_file, top_k, split_char, col_score)
#~ split_char = ',|\n'
#~ col_score = 5
#~ input_file = "data.training/dia.urine.2018_03_29/peaks.denovo.csv"


# filter features of single-feature (DDA-like) scan or multi-feature scan (DIA)
def filter_multifeature(input_file):
  """TODO(nh2tran): docstring."""

  print("".join(["="] * 80)) # section-separating line
  print("filter_multifeature()")

  print('input_file = ', input_file)
  output_file_1 = input_file + '.1fea'
  output_file_2 = input_file + '.2fea'
  print('output_file_1 = ', output_file_1)
  print('output_file_2 = ', output_file_2)

  # read feature and record feature_dict, scan_dict
  with open(input_file, 'r') as input_handle:
    # header
    header_line = input_handle.readline()
    col_feature_id = deepnovo_config.col_feature_id
    col_scan_list = deepnovo_config.col_scan_list
    feature_dict = {}
    scan_dict = {}
    # read feature and record feature_dict, scan_dict
    for line in input_handle:
      line_split = re.split(',|\n', line)
      feature_id = line_split[col_feature_id]
      scan_list = re.split(';', line_split[col_scan_list])
      feature_dict[feature_id] = {}
      feature_dict[feature_id]['line'] = line
      feature_dict[feature_id]['scan_list'] = scan_list
      for scan_id in scan_list:
        if scan_id in scan_dict:
          scan_dict[scan_id]['feature_list'].append(feature_id)
        else:
          scan_dict[scan_id] = {}
          scan_dict[scan_id]['feature_list'] = [feature_id]

  print('Total scan count = ', len(scan_dict))
  print('  Scan with single-feature = ',
        sum([1 if (len(scan['feature_list'])==1) else 0 for _, scan in scan_dict.iteritems()]))
  print('  Scan with multi-feature = ',
        sum([1 if (len(scan['feature_list'])>=2) else 0 for _, scan in scan_dict.iteritems()]))

  # write feature to separate files,
  # depending on its scan is single-feature (DDA-like) or multi-feature (DIA)
  single_feature_count = 0
  multi_feature_count = 0
  with open(output_file_1, 'w') as output_handle_1:
    with open(output_file_2, 'w') as output_handle_2:
      # header
      print(header_line, file=output_handle_1, end="")
      print(header_line, file=output_handle_2, end="")
      for feature_id, feature in feature_dict.iteritems():
        # assuming all scans are single-feature
        output_handle = output_handle_1
        single_feature_count += 1
        # at least 1 scan is multi-feature
        #~ for scan_id in feature['scan_list']:
          #~ if len(scan_dict[scan_id]['feature_list']) >= 2:
            #~ output_handle = output_handle_2
            #~ multi_feature_count += 1
            #~ single_feature_count -= 1
            #~ break
        # average feature count of scans
        feature_count = sum([len(scan_dict[scan_id]['feature_list']) for scan_id in feature['scan_list']])
        feature_count /= float(len(feature['scan_list']))
        if feature_count >= 2:
          output_handle = output_handle_2
          multi_feature_count += 1
          single_feature_count -= 1
        print(feature['line'], file=output_handle, end="")

  print('Total feature count = ', len(feature_dict))
  print('Feature with single-feature scans = ', single_feature_count)
  print('Feature with at least 1 multi-feature scans = ', multi_feature_count)

#~ input_file = "data.training/dia.urine.2018_03_29/testing_12.feature.csv"
#~ filter_multifeature(input_file)

