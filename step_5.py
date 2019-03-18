from __future__ import print_function

import csv
import re
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import Levenshtein
import multiprocessing
num_processes = 8
import time


# I/O files 
psm_file = "data.training/aa.hla.bassani.nature_2016.mel_15/step_4.DB search psm round_2_FDR_1%.csv"
netmhc_file = "data.training/aa.hla.bassani.nature_2016.mel_15/step_5.denovo_peptide.NetMHCpan.xls.csv"
WEAK_BINDING = 2.0 # NetMHC weak binding rank
STRONG_BINDING = 0.5 # NetMHC strong binding rank
db_fasta_file = "data/uniprot_sprot.human.plus_contaminants.fasta"
labeled_feature_file = "data.training/aa.hla.bassani.nature_2016.mel_15/feature.csv.labeled"
output_neoantigen_criteria = "data.training/aa.hla.bassani.nature_2016.mel_15/step_5.output_neoantigen_criteria.csv"


CODON_AA = { # dictionary {codon: aa}
  'TTT':'F',
  'TTC':'F',
  'TTA':'L',
  'TTG':'L',
  'TCT':'S',
  'TCC':'S',
  'TCA':'S',
  'TCG':'S',
  'TAT':'Y',
  'TAC':'Y',
  'TAA':'X',
  'TAG':'X',
  'TGT':'C',
  'TGC':'C',
  'TGA':'X',
  'TGG':'W',
  'CTT':'L',
  'CTC':'L',
  'CTA':'L',
  'CTG':'L',
  'CCT':'P',
  'CCC':'P',
  'CCA':'P',
  'CCG':'P',
  'CAT':'H',
  'CAC':'H',
  'CAA':'Q',
  'CAG':'Q',
  'CGT':'R',
  'CGC':'R',
  'CGA':'R',
  'CGG':'R',
  'ATT':'I',
  'ATC':'I',
  'ATA':'I',
  'ATG':'M',
  'ACT':'T',
  'ACC':'T',
  'ACA':'T',
  'ACG':'T',
  'AAT':'N',
  'AAC':'N',
  'AAA':'K',
  'AAG':'K',
  'AGT':'S',
  'AGC':'S',
  'AGA':'R',
  'AGG':'R',
  'GTT':'V',
  'GTC':'V',
  'GTA':'V',
  'GTG':'V',
  'GCT':'A',
  'GCC':'A',
  'GCA':'A',
  'GCG':'A',
  'GAT':'D',
  'GAC':'D',
  'GAA':'E',
  'GAG':'E',
  'GGT':'G',
  'GGC':'G',
  'GGA':'G',
  'GGG':'G'}

AA_CODON ={} # dictionary {aa: list of codons}
for codon, aa in CODON_AA.iteritems():
  if aa in AA_CODON:
    AA_CODON[aa].append(codon)
  else:
    AA_CODON[aa] = [codon]

AA_PAIRWISE_DISTANCE = {} # dictionary {(aa1, aa2): min_distance}
for aa1 in AA_CODON:
  for aa2 in AA_CODON:
    if (aa1, aa2) not in AA_PAIRWISE_DISTANCE:
      min_distance = 3
      for codon1 in AA_CODON[aa1]:
        for codon2 in AA_CODON[aa2]:
          distance = Levenshtein.hamming(codon1, codon2)
          assert distance <= 3, "Error: codon distance > 3"
          min_distance = min(min_distance, distance)
      AA_PAIRWISE_DISTANCE[(aa1, aa2)] = min_distance
      AA_PAIRWISE_DISTANCE[(aa2, aa1)] = min_distance

# a mutation pair (aa1, aa2) is missense if their codons are different by 1 nucleotide
AA_PAIR_MISSENSE = [(aa1, aa2) for (aa1, aa2), min_distance in AA_PAIRWISE_DISTANCE.iteritems()
                    if min_distance == 1]
# for now, remove N-D, Q-E because not sure mutations or modifications
AA_PAIR_MISSENSE.remove(('N', 'D'))
AA_PAIR_MISSENSE.remove(('D', 'N'))
AA_PAIR_MISSENSE.remove(('Q', 'E'))
AA_PAIR_MISSENSE.remove(('E', 'Q'))


def drop_mod_peaks(peptide):
  peptide = peptide.replace("M(+15.99)", "M")
  peptide = peptide.replace("N(+.98)", "N")
  peptide = peptide.replace("Q(+.98)", "Q")
  return peptide


def read_denovo_psm(psm_file):

  print("read_denovo_psm()")
  print("psm_file:", psm_file)

  # store PSM of denovo peptides in a dictionary 
  # {peptide: {'num_psm': , 'total_score': , 'total_abundance'}}
  denovo_peptide_psm = {}
  with open(psm_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      accession = drop_mod_peaks(row['Accession'])
      if accession == 'DENOVO':
        peptide = drop_mod_peaks(row['Peptide'])
        score = float(row['-10lgP'])
        abundance = float(row['Area']) if row['Area'] else 0
        if peptide not in denovo_peptide_psm:
          denovo_peptide_psm[peptide] = {'num_psm': 1,
                                         'total_score': score,
                                         'total_abundance': abundance}
        else:
          denovo_peptide_psm[peptide]['num_psm'] += 1
          denovo_peptide_psm[peptide]['total_score'] += score
          denovo_peptide_psm[peptide]['total_abundance'] += abundance

  print("Number of denovo peptides:", len(denovo_peptide_psm))
  num_psm_list = [x['num_psm'] for x in denovo_peptide_psm.values()]
  print("Number of denovo peptides with >= 1 psm: ", len([x for x in num_psm_list if x >= 1]))
  print("Number of denovo peptides with >= 2 psm: ", len([x for x in num_psm_list if x >= 2]))
  print("Number of denovo peptides with >= 3 psm: ", len([x for x in num_psm_list if x >= 3]))
  print("Number of denovo peptides with >= 4 psm: ", len([x for x in num_psm_list if x >= 4]))
  print()

  return denovo_peptide_psm


def read_denovo_netmhc(netmhc_file):

  print("read_denovo_netmhc()")
  print("netmhc_file:", netmhc_file)

  # store NetMHC predictions of denovo peptides in a dictionary 
  # {peptide: {'best_nM': , 'best_rank': , 'is_weak_binding': , 'is_strong_binding': }}
  denovo_peptide_netmhc = {}
  with open(netmhc_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      peptide = row['Peptide']
      if peptide not in denovo_peptide_netmhc:
        best_nM = min([float(row[x]) for x in ['nM1', 'nM2', 'nM3', 'nM4']])
        best_rank = min([float(row[x]) for x in ['Rank1', 'Rank2', 'Rank3', 'Rank4']])
        is_weak_binding = int(best_rank <= WEAK_BINDING)
        is_strong_binding = int(best_rank <= STRONG_BINDING)
        denovo_peptide_netmhc[peptide] = {
          'best_nM': best_nM,
          'best_rank': best_rank,
          'is_weak_binding': is_weak_binding,
          'is_strong_binding': is_strong_binding}
      else:
        print("Warning: duplicate peptide found in denovo_peptide_netmhc:", peptide)

  print("Number of denovo peptides:", len(denovo_peptide_netmhc))
  print("Number of denovo peptides with weak binding: ", sum([x['is_weak_binding'] for x in denovo_peptide_netmhc.values()]))
  print("Number of denovo peptides with strong binding: ", sum([x['is_strong_binding'] for x in denovo_peptide_netmhc.values()]))
  print()

  return denovo_peptide_netmhc


def hamming1_align((peptide, protein_list)):

  query = peptide.replace('I', 'L')
  query_length = len(query)
  match_list = []
  for protein in protein_list:
    subject = protein['seq_ItoL']
    subject_length = len(subject)

    # First, find candidate locations by pigeonhole principle:
    # if hamming distance is 1, the left or right half must be exact match
    # Then, calculate hamming distance at candidate locations and return those equal to 1
    query_left = query[:query_length/2]
    query_right = query[query_length/2:]
    left_index = [x.start() for x in re.finditer(query_left, subject)]
    right_index = [x.start() for x in re.finditer(query_right, subject)]
    right_index = [(x - query_length/2) for x in right_index]
    candidate_index = left_index + right_index
    candidate_index = [x for x in candidate_index if x >= 0 and (x + query_length) <= subject_length]
    hamming1_index = [x for x in candidate_index
                      if Levenshtein.hamming(query, subject[x : (x + query_length)]) == 1]

    if hamming1_index:
      match_list += [{'protein': protein['name'],
                      'wildtype': protein['seq'][x : (x + len(peptide))]}
                     for x in hamming1_index]

  return peptide, match_list


def find_denovo_mutation(denovo_peptide_list, db_fasta_file, labeled_feature_file):

  print("find_denovo_mutation()")
  print("db_fasta_file:", db_fasta_file)
  print("labeled_feature_file:", labeled_feature_file)

  print("Number of denovo peptides: ", len(denovo_peptide_list))

  with open(db_fasta_file, 'r') as input_fasta_handle:
    record_list = list(SeqIO.parse(input_fasta_handle, "fasta"))
    print("Number of protein sequences in the database: ", len(record_list))
  protein_list = [{'name': str(record.name), 
                   'seq': str(record.seq),
                   'seq_ItoL': str(record.seq).replace('I', 'L')}
                  for record in record_list]

  db_peptide_set = set()
  with open(labeled_feature_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      peptide = drop_mod_peaks(row['seq'])
      db_peptide_set.add(peptide)
  print("Number of db peptides identified at step 1: ", len(db_peptide_set))

  print("Align denovo peptides against protein sequences with 1 mismatch ...")
  start_time = time.time()
  pool = multiprocessing.Pool(processes=num_processes)
  search_list = [(peptide, protein_list) for peptide in denovo_peptide_list]
  result_list = pool.map(hamming1_align, search_list)
  print(time.time() - start_time, "seconds")

  print("Annotate denovo mutations and their wildtypes ...")
  denovo_peptide_mutation = {}
  for peptide, match_list in result_list:
    missense_list = []
    peptide_ItoL = peptide.replace('I', 'L')
    for match in match_list:
      wildtype = match['wildtype']
      wildtype_ItoL = wildtype.replace('I', 'L')
      assert len(peptide_ItoL) == len(wildtype_ItoL), "Error: peptide and wildtype lengths not matched"
      mutation_index = [x for x in range(len(peptide_ItoL)) if peptide_ItoL[x] != wildtype_ItoL[x]]
      assert len(mutation_index) >= 1, "Error: more than 1 mutation found"
      mutation_index = mutation_index[0]
      mutation_aa = peptide[mutation_index]
      mutation_wildtype = wildtype[mutation_index]
      match['mutation_pos'] = mutation_index + 1
      match['mutation_aa'] = '->'.join([mutation_wildtype, mutation_aa])
      match['is_missense'] = int((mutation_aa, mutation_wildtype) in AA_PAIR_MISSENSE)
      match['is_db'] = int(wildtype in db_peptide_set)
      match['is_missense_db'] = match['is_missense'] * match['is_db']

    num_hits = len(match_list)
    num_missense = len([x for x in match_list if x['is_missense'] == 1])
    num_db = len([x for x in match_list if x['is_db'] == 1])
    num_missense_db = len([x for x in match_list if x['is_missense_db'] == 1])
    denovo_peptide_mutation[peptide] = {'num_hits': num_hits,
                                        'num_missense': num_missense,
                                        'num_db': num_db,
                                        'num_missense_db': num_missense_db,
                                        'match_list': match_list}

  print("Number of denovo peptides with > 0 hits:",
        len([x for x in denovo_peptide_mutation.values() if x['num_hits'] > 0]))
  print("Number of denovo peptides with > 0 missense hits:",
        len([x for x in denovo_peptide_mutation.values() if x['num_missense'] > 0]))
  print("Number of denovo peptides with > 0 missense db hits:",
        len([x for x in denovo_peptide_mutation.values() if x['num_missense_db'] > 0]))
  print()

  return denovo_peptide_mutation
        

if __name__ == '__main__':

  denovo_peptide_psm = read_denovo_psm(psm_file)
  denovo_peptide_netmhc = read_denovo_netmhc(netmhc_file)
  denovo_peptide_list = denovo_peptide_psm.keys()
  denovo_peptide_mutation = find_denovo_mutation(denovo_peptide_list,
                                                 db_fasta_file,
                                                 labeled_feature_file)

  print("Write neoantigen criteria:")
  print("output_neoantigen_criteria:", output_neoantigen_criteria)
  print()
  with open(output_neoantigen_criteria, 'w') as output_handle:
    fieldnames = ['peptide',
                  'num_psm',
                  'total_score',
                  'total_abundance',
                  'best_nM',
                  'best_rank',
                  'is_weak_binding',
                  'is_strong_binding',
                  'num_hits',
                  'num_missense',
                  'num_db',
                  'num_missense_db',
                  'match_list']
    csv_writer = csv.DictWriter(output_handle, fieldnames=fieldnames, delimiter=',')
    csv_writer.writeheader()
    for peptide in denovo_peptide_list:
      row = {'peptide': peptide}
      row.update(denovo_peptide_psm[peptide])
      row.update(denovo_peptide_netmhc[peptide])
      row.update(denovo_peptide_mutation[peptide])
      csv_writer.writerow(row)






