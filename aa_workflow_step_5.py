from __future__ import print_function


import sys
import csv
import re
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import Levenshtein
import multiprocessing
num_processes = 8
import time


WEAK_BINDING = 2.0 # NetMHC weak binding rank
STRONG_BINDING = 0.5 # NetMHC strong binding rank

AA_3_to_1 = {
  'Ala':'A',
  'Arg':'R',
  'Asn':'N',
  'Asp':'D',
  'Cys':'C',
  'Glu':'E',
  'Gln':'Q',
  'Gly':'G',
  'His':'H',
  'Ile':'I',
  'Leu':'L',
  'Lys':'K',
  'Met':'M',
  'Phe':'F',
  'Pro':'P',
  'Ser':'S',
  'Thr':'T',
  'Trp':'W',
  'Tyr':'Y',
  'Val':'V'}

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
  print()

  return denovo_peptide_psm


def read_netmhc(netmhc_file):

  print("read_netmhc()")
  print("netmhc_file:", netmhc_file)

  # store NetMHC predictions of denovo peptides in a dictionary 
  # {peptide: {'best_nM': , 'best_rank': , 'is_weak_binding': , 'is_strong_binding': }}
  peptide_netmhc = {}
  with open(netmhc_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      peptide = row['Peptide']
      if peptide not in peptide_netmhc:
        best_nM = min([float(row[x]) for x in ['nM1', 'nM2', 'nM3', 'nM4']])
        best_rank = min([float(row[x]) for x in ['Rank1', 'Rank2', 'Rank3', 'Rank4']])
        is_weak_binding = int(best_rank <= WEAK_BINDING)
        is_strong_binding = int(best_rank <= STRONG_BINDING)
        peptide_netmhc[peptide] = {
          'best_nM': best_nM,
          'best_rank': best_rank,
          'is_weak_binding': is_weak_binding,
          'is_strong_binding': is_strong_binding}
      else:
        print("Warning: duplicate peptide found in peptide_netmhc:", peptide)

  print("Number of peptides:", len(peptide_netmhc))
  print("Number of peptides with weak binding: ", sum([x['is_weak_binding'] for x in peptide_netmhc.values()]))
  print("Number of peptides with strong binding: ", sum([x['is_strong_binding'] for x in peptide_netmhc.values()]))
  print()

  return peptide_netmhc


def read_fasta(fasta_file,
               get_uniprot_id=False,
               get_enst_id=False,
               get_gene_name=False):

  print("read_fasta()")
  print("fasta_file:", fasta_file)
  print("get_uniprot_id:", get_uniprot_id)
  print("get_enst_id:", get_enst_id)
  print("get_gene_name:", get_gene_name)

  with open(fasta_file, 'r') as file_handle:
    record_list = list(SeqIO.parse(file_handle, "fasta"))
    protein_list = []
    for record in record_list:
      uniprot_id = ''
      enst_id = ''
      gene_name = ''
      name = str(record.name)
      if get_uniprot_id:
        uniprot_id = name.split('|')[1]
      if get_enst_id:
        enst_id = name
      if get_gene_name:
        description_list = str(record.description).strip().split(' ')
        gene_name_list = [x for x in description_list if 'GN=' in x]
        if len(gene_name_list) == 1:
          gene_name = gene_name_list[0].split('=')[1]
      seq = str(record.seq)
      protein_list.append({'name': name,
                           'uniprot_id': uniprot_id,
                           'enst_id': enst_id,
                           'gene_name': gene_name,
                           'seq': seq})

  print("Number of protein sequences in the fasta file: ", len(protein_list))
  print()

  return protein_list


def read_db_peptide(labeled_feature_file):

  print("read_db_peptide()")
  print("labeled_feature_file:", labeled_feature_file)

  db_peptide_set = set()
  with open(labeled_feature_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      peptide = drop_mod_peaks(row['seq'])
      db_peptide_set.add(peptide)
  print("Number of db peptides identified at step 1: ", len(db_peptide_set))
  print()

  return db_peptide_set


def hamming1_align((peptide, protein_list)):

  # I and L are considered the same in this alignment
  query = peptide.replace('I', 'L')
  query_length = len(query)
  match_list = []
  for protein in protein_list:
    subject = protein['seq'].replace('I', 'L')
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
      match_list += [{'protein': protein, 'match_index': index}
                      for index in hamming1_index]

  return peptide, match_list


def find_mutation(peptide_list, protein_list):

  print("find_mutation()")

  print("Align peptides against protein sequences with 1 mismatch ...")
  print("Number of peptides: ", len(peptide_list))
  print("Number of protein sequences:", len(protein_list))
  print("I and L are considered the same in this alignment")
  start_time = time.time()
  pool = multiprocessing.Pool(processes=num_processes)
  search_list = [(peptide, protein_list) for peptide in peptide_list]
  result_list = pool.map(hamming1_align, search_list)
  print(time.time() - start_time, "seconds")
  print()

  peptide_mutation = {}
  protein_mutation = {}
  for peptide, match_list in result_list:
    missense_list = []
    peptide_length = len(peptide)
    peptide_ItoL = peptide.replace('I', 'L')
    for match in match_list:
      protein = match['protein']
      match_index = match['match_index']

      wildtype = protein['seq'][match_index : (match_index + peptide_length)]
      wildtype_ItoL = wildtype.replace('I', 'L')
      mutation_index = [x for x in range(len(peptide_ItoL)) if peptide_ItoL[x] != wildtype_ItoL[x]]
      assert len(mutation_index) == 1, "Error: not 1 mutation found"
      mutation_index = mutation_index[0]
      mutation_wildtype = wildtype[mutation_index]
      mutation_aa = peptide[mutation_index]
      match['wildtype'] = wildtype
      match['mutation_pos'] = mutation_index + 1
      match['mutation_wt'] = mutation_wildtype
      match['mutation_aa'] = mutation_aa
      match['is_missense'] = int((mutation_aa, mutation_wildtype) in AA_PAIR_MISSENSE)
      not_flanking = int(match['mutation_pos'] != 1 and match['mutation_pos'] != len(peptide))
      match['is_missense_not_flanking'] = match['is_missense'] * not_flanking

      if match['is_missense_not_flanking']:
        protein_mutation_entry = {'peptide': peptide, 'match_index': match['match_index']}
        if not protein['name'] in protein_mutation:
          protein_mutation[protein['name']] = [protein_mutation_entry]
        else:
          protein_mutation[protein['name']].append(protein_mutation_entry)

    num_hits = len(match_list)
    num_missense = len([x for x in match_list if x['is_missense'] == 1])
    num_missense_not_flanking = len([x for x in match_list if x['is_missense_not_flanking'] == 1])
    peptide_mutation[peptide] = {'num_hits': num_hits,
                                 'num_missense': num_missense,
                                 'num_missense_not_flanking': num_missense_not_flanking,
                                 'match_list': match_list}

  print("Number of denovo peptides with >= 1 hits:",
        len([x for x in peptide_mutation.values() if x['num_hits'] >= 1]))
  print("Number of denovo peptides with >= 1 missense hits:",
        len([x for x in peptide_mutation.values() if x['num_missense'] >= 1]))
  print("Number of denovo peptides with >= 1 missense, not flanking hits:",
        len([x for x in peptide_mutation.values() if x['num_missense_not_flanking'] >= 1]))
  print()

  return peptide_mutation, protein_mutation
        

def read_missense_snp(snp_file, snp_enst_fasta, snp_sample_id):

  print("read_missense_snp()")
  print("snp_file:", snp_file)
  print("snp_enst_fasta:", snp_enst_fasta)
  print("snp_sample_id:", snp_sample_id)

  # read SNP file
  snp_list = []
  with open(snp_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      mutation_type = row['Effect']
      if mutation_type == 'missense_variant' and snp_sample_id == row['Sample ID']:
        enst_id = row['ENSEMBL Transcript ID']
        mutation_change = row['Aa change']
        snp_list.append({'enst_id': enst_id, 'mutation_change': mutation_change})
  print("Number of SNPs:", len(snp_list))
  print()

  # cross-check snp_list and snp_enst_fasta for enst_id, location of mutated amino acid
  # because some transcripts were removed or updated, so their SNPs are no longer correct
  protein_list = read_fasta(snp_enst_fasta, get_enst_id=True)
  # clean letter 'X' from the 1st position of some enst protein sequences
  for protein in protein_list:
    if protein['seq'][0] == 'X':
      protein['seq'] = protein['seq'][1:]
  num_not_missense = 0
  num_protein_confirmed = 0
  snp_confirmed_list = []
  for snp in snp_list:
    # example: Pro575Leu; note that the location is 1-based, not 0-based
    aa_3letter_ref = snp['mutation_change'][:3]
    aa_loc = int(snp['mutation_change'][3:-3])
    aa_3letter_alt = snp['mutation_change'][-3:]
    aa_ref = AA_3_to_1[aa_3letter_ref]
    aa_alt = AA_3_to_1[aa_3letter_alt]
    if (aa_ref, aa_alt) not in AA_PAIR_MISSENSE:
      num_not_missense += 1
    protein_confirmed = False
    for protein in protein_list:
      if protein['enst_id'] == snp['enst_id']:
        num_protein_confirmed += 1
        if aa_loc-1 < len(protein['seq']) and aa_ref == protein['seq'][aa_loc-1]:
          snp_confirmed_list.append({'enst_id':snp['enst_id'],
                                     'aa_loc': aa_loc,
                                     'aa_ref': aa_ref,
                                     'aa_alt': aa_alt})

  print("len(snp_list):", len(snp_list))
  print("Warning: num_not_missense", num_not_missense)
  print("num_protein_confirmed:", num_protein_confirmed)
  print("len(snp_confirmed_list):", len(snp_confirmed_list))
  print()

  return snp_confirmed_list, protein_list


def match_peptide_snp(peptide_list, snp_file, snp_enst_fasta, snp_sample_id):

  print('match_peptide_snp()')

  snp_list, protein_list = read_missense_snp(snp_file, snp_enst_fasta, snp_sample_id)
  peptide_mutation, _ = find_mutation(peptide_list, protein_list)
  peptide_snp = {}
  for peptide, mutation in peptide_mutation.iteritems():
    peptide_snp[peptide] = {'snp_list': []}
    if mutation['num_hits'] > 0:
      for match in mutation['match_list']:
        enst_id = match['protein']['enst_id']
        match_index = match['match_index']
        for snp in snp_list:
          if (enst_id == snp['enst_id']
              and match_index + match['mutation_pos'] == snp['aa_loc']
              and match['mutation_wt'] == snp['aa_ref']
              and match['mutation_aa'].replace('I', 'L') == snp['aa_alt'].replace('I', 'L')):
            match_snp = snp
            match_snp.update({'wildtype': match['wildtype']})
            peptide_snp[peptide]['snp_list'].append(match_snp)

  num_peptide_snp = len([x for x in peptide_snp.values() if x['snp_list']])
  print('Number of peptide mutations match to SNPs:', num_peptide_snp)
  for peptide in peptide_snp:
    if peptide_snp[peptide]['snp_list']:
      print(peptide, peptide_snp[peptide]['snp_list'])
  print()

  return peptide_snp


def step_5(psm_file, netmhc_file, db_fasta_file, labeled_feature_file,
           snp_file, snp_enst_fasta, snp_sample_id,
           output_neoantigen_criteria, output_protein_mutation):

  print("".join(["="] * 80)) # section-separating line
  print("step_5()")

  denovo_psm = read_denovo_psm(psm_file)
  if netmhc_file:
    denovo_netmhc = read_netmhc(netmhc_file)
  else:
    denovo_netmhc = None
  denovo_peptide_list = denovo_psm.keys()

  print("Find denovo mutations with respect to the reference fasta:")
  protein_list = read_fasta(db_fasta_file)
  denovo_mutation, protein_mutation = find_mutation(denovo_peptide_list, protein_list)

  print("Write protein with missense and not flanking mutations:")
  print("output_protein_mutation:", output_protein_mutation)
  print()
  with open(output_protein_mutation, 'w') as output_handle:
    fieldnames = ['protein_name', 'num_peptide', 'peptide_list']
    csv_writer = csv.DictWriter(output_handle, fieldnames=fieldnames, delimiter=',')
    csv_writer.writeheader()
    for protein_name, peptide_list in protein_mutation.iteritems():
      row = {'protein_name': protein_name,
             'num_peptide': len(peptide_list),
             'peptide_list': peptide_list}
      csv_writer.writerow(row)

  print("Find wildtypes in identified db peptides")
  db_peptide_set = read_db_peptide(labeled_feature_file)
  for peptide in denovo_mutation:
    num_db = 0
    for match in denovo_mutation[peptide]['match_list']:
      match['is_db'] = int(match['wildtype'] in db_peptide_set)
      num_db += match['is_db']
    denovo_mutation[peptide]['num_db'] = num_db
  print("Number of denovo peptides with >= 1 wildtype hits:",
        len([x for x in denovo_mutation.values() if x['num_db'] >= 1]))
  print()

  print("Find denovo mutations match to SNPs:")
  denovo_snp = match_peptide_snp(denovo_peptide_list, snp_file, snp_enst_fasta, snp_sample_id)
  print()

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
                  'num_missense_not_flanking',
                  'num_db',
                  'match_list',
                  'snp_list']
    csv_writer = csv.DictWriter(output_handle, fieldnames=fieldnames, delimiter=',')
    csv_writer.writeheader()
    for peptide in denovo_peptide_list:
      row = {'peptide': peptide}
      row.update(denovo_psm[peptide])
      if denovo_netmhc is not None:
        row.update(denovo_netmhc[peptide])
      row.update(denovo_mutation[peptide])
      row.update(denovo_snp[peptide])
      for match in row['match_list']:
        match['protein'] = match['protein']['name']
      csv_writer.writerow(row)

  print("Selection criteria: >= 1 missense, not flanking hits AND >= 2 psm")
  num_selection = len([peptide for peptide in denovo_peptide_list
                       if denovo_mutation[peptide]['num_missense_not_flanking'] >= 1
                       and denovo_psm[peptide]['num_psm'] >= 2])
  print("num_selection :", num_selection)









