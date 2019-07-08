from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

from deepnovo_preprocess import *
from deepnovo_postprocess import *
import aa_workflow_step_4_2
import aa_workflow_step_5


data_fasta_dir = "data.fasta/"
patient_id = "MM15"
data_training_dir = "data.training/aa.hla.bassani.nature_2016.mel_15.class_1/train.exclude_mel_15/"
num_fractions = 10
model_dir = "train.mel_15.class_1" # create this empty folder at the same level as Python scripts.


# ================================================================================
# Workflow of neoantigen discovery by personalized de novo sequencing.
# ================================================================================

# Step-by-step instructions based on the following example dataset:

#   Patient Mel-16 (Bassani-Sternberg et al., Nature Communication, 2016)
#   HLA class 1: 12 raw files, 1 failed to run PEAKS

#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_1_A
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_1_B
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_2_A
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_2_B
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_3_A
#     20141212_QEp7_MiBa_SA_HLA-I-p_MM16_3_B
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_1_A_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_1_B_1, failed
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_2_A_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_2_B_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_3_A_1
#     20141213_QEp7_MiBa_SA_HLA-I-p_MM16_3_B_1




# ================================================================================
# Step 1: Build the immunopeptidome of the patient.
# ================================================================================

# This step 1 took about ?? hours on a laptop with 4 CPU cores i7, 16 GB memory

# ================================================================================
# Step 1.1: Run PEAKS X DB search on the raw files with the following parameters:
# ================================================================================

#     Enzyme: None
#     Instrument: Orbi-Orbi
#     Fragment: HCD
#     Acquisition: DDA

#     Parent Mass Error Tolerance: 15.0 ppm
#     Fragment Mass Error Tolerance: 0.05 Da
#     Precursor Mass Search Type: monoisotopic
#     Enzyme: None
#     Digest Mode: Unspecific
#     Max Missed Cleavages: 100
#     Variable Modifications:
#       Oxidation (M): 15.99
#       Deamidation (NQ): 0.98
#     Max Variable PTM Per Peptide: 3
#     Database: uniprot_sprot.human
#     Taxon: All
#     Contaminant Database: contaminants_maxquant
#     Searched Entry: 20488
#     FDR Estimation: Enabled
#     Merge Options: no merge
#     Precursor Options: corrected
#     Charge Options: no correction
#     Filter Options: no filter
#     Process: true
#     Associate chimera: no




# ================================================================================
# Step 1.2: Set FDR 1.0%.
# ================================================================================

# The number of MS/MS spectra is "694565", the number of peptide-spectrum matches (PSMs) is "207332", the number of peptide sequences is "26594".




# ================================================================================
# Step 1.3: Right-click on the DB search node "??", select "Deep Denovo Export".
# ================================================================================

# We will get the following 8 pairs of csv and mgf files in the PEAKS project folder:

#       export_0.csv, export_0.mgf
#       export_1.csv, export_1.mgf
#       export_2.csv, export_2.mgf
#       export_3.csv, export_3.mgf
#       export_4.csv, export_4.mgf
#       export_5.csv, export_5.mgf
#       export_6.csv, export_6.mgf
#       export_7.csv, export_7.mgf
#       export_8.csv, export_8.mgf
#       export_9.csv, export_9.mgf
#       export_10.csv, export_10.mgf




# ================================================================================
# Step 2: Train personalized DeepNovo model.
# ================================================================================

# This step 2 took about 12 hours on a server with GPU Titan X, 32 GB memory

# Note that you will need to specify the paths to your own data and model folders when you run the Python scripts. The following scripts just show examples of my data and model folders.

# ================================================================================
# Step 2.1: Prepare the training data.
# ================================================================================

# Run merge_mgf_file() and merge_feature_file()
# ======================= UNCOMMENT and RUN ======================================
# ~ folder_path = data_training_dir
# ~ fraction_list = range(0, num_fractions)
# ~ merge_mgf_file(
    # ~ input_file_list=[folder_path + "export_" + str(i) + ".mgf" for i in fraction_list],
    # ~ fraction_list=fraction_list,
    # ~ output_file=folder_path + "spectrum.mgf")
# ~ merge_feature_file(
    # ~ input_file_list=[folder_path + "export_" + str(i) + ".csv" for i in fraction_list],
    # ~ fraction_list=fraction_list,
    # ~ output_file=folder_path + "feature.csv")
# ================================================================================
# We will get two output files in the same folder: "spectrum.mgf" and "feature.csv".
# Both functions also report the number of entries that have been processed: "counter = 694565".
# That number should be the same as the total number of MS/MS spectra from the raw files.

# Run split_feature_unlabel()
# ======================= UNCOMMENT and RUN ======================================
# ~ input_feature_file = data_training_dir + "feature.csv"
# ~ split_feature_unlabel(input_feature_file)
# ================================================================================
# It will split the "feature.csv" into 2 files: "feature.csv.labeled" and "feature.csv.unlabeled".
# It also reports the number of labeled and unlabel features: "num_labeled = 207332" and "num_unlabeled = 487233".
# Note that "207332" is also the number of PSMs reported at FDR 1.0% in Step 1.

# Run calculate_mass_shift_ppm() and correct_mass_shift_ppm()
# ======================= UNCOMMENT and RUN ======================================
# ~ labeled_feature_file = data_training_dir + "feature.csv.labeled"
# ~ ppm = calculate_mass_shift_ppm(labeled_feature_file)
# ~ input_feature_file = data_training_dir + "feature.csv.labeled"
# ~ correct_mass_shift_ppm(input_feature_file, ppm)
# ~ input_feature_file = data_training_dir + "feature.csv"
# ~ correct_mass_shift_ppm(input_feature_file, ppm)
# ================================================================================
# The mass shift is calculated from "feature.csv.labeled".
# The mass shift ppm (part per million) is reported as: "mean_precursor_ppm = 7.07514819678".
# Then mass is corrected for 2 files: "feature.csv.labeled.mass_corrected" and "feature.csv.mass_corrected".

# Run split_feature_training_noshare()
# ======================= UNCOMMENT and RUN ======================================
# ~ input_feature_file = data_training_dir + "feature.csv.labeled.mass_corrected"
# ~ proportion = [0.90, 0.05, 0.05]
# ~ split_feature_training_noshare(input_feature_file, proportion)
# ================================================================================
# It will split "feature.csv.labeled.mass_corrected" into train/valid/test sets with "proportion = [0.9, 0.05, 0.05]".
# Those 3 sets do not share common peptides.
# Their sizes are also reported.
#   "num_total = 207332"
#   "num_unique = 26656"
#   "num_train = 185823"
#   "num_valid = 10900"
#   "num_test = 10609"




# ================================================================================
# Step 2.2: Training DeepNovo model.
# ================================================================================

# Run DeepNovo training
# The training will stop after 10 epoch. The model with best performance on the valid set, "ckpt-16200" is saved in the model folder "train.mel_16.class_1".
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --train"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--train_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--train_feature", data_training_dir + "feature.csv.labeled.mass_corrected.train.noshare"]
# ~ command += ["--valid_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--valid_feature", data_training_dir + "feature.csv.labeled.mass_corrected.valid.noshare"]
# ~ command += ["--reset_step"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================

# Run DeepNovo testing
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test_true_feeding"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--test_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--test_feature", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --search_denovo"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--denovo_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--denovo_feature", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test"]
# ~ command += ["--target_file", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command += ["--predicted_file", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================
# The testing accuracy at the amino acid (AA) and peptide levels will be reported as following:
#   "precision_AA_mass_db  = 0.8425"
#   "precision_peptide_mass_db  = 0.6430"




# ================================================================================
# Step 3: Perform personalized de novo sequencing with DeepNovo.
# ================================================================================

# This step 3 took about 5 hours on a server with GPU Titan X, 32 GB memory

# Run DeepNovo de novo sequencing on all features (label and unlabeled)
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --search_denovo"]
# ~ command += ["--train_dir", model_dir]
# ~ command += ["--denovo_spectrum", data_training_dir + "spectrum.mgf"]
# ~ command += ["--denovo_feature", data_training_dir + "feature.csv.mass_corrected"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================
# The de novo results will be written to the file "feature.csv.mass_corrected.deepnovo_denovo".
# The tool will also report the number of features that have been processed:
#   "Total spectra: 694565"
#     "read: 690354"
#     "skipped: 4211"
#       "by mass: 4211"




# ================================================================================
# Step 4: Quality control.
# ================================================================================

# ================================================================================
# Step 4.1: Post-process de novo results to improve their accuracy. 
# ================================================================================

# Run select_top_score()
# This script selects a threshold of de novo confidence scores and uses it to filter de novo results.
# The score threshold is calculated based on a 95% cutoff of the testing accuracy obtained at the end of Step 2 above.
# ======================= UNCOMMENT and RUN ======================================
# ~ accuracy_cutoff = 0.95
# ~ accuracy_file = data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo.accuracy"
# ~ score_cutoff = find_score_cutoff(accuracy_file, accuracy_cutoff)
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo"
# ~ output_file = input_file + ".top95"
# ~ select_top_score(input_file, output_file, score_cutoff)
# ================================================================================
# After this step we'll get the file "feature.csv.mass_corrected.deepnovo_denovo.top95".
# The score cutoff and the number of selected features will also be reported:
#   "score_cutoff =  -0.5"
#   "total_feature =  690354"
#   "select_feature =  233589"

# Run convert_I_to_L()
# This script converts I (Isoleucine) to L (Leucine) in all de novo peptides, because de novo sequencing is not able to distinguish them.
# ======================= UNCOMMENT and RUN ======================================
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95"
# ~ output_file = input_file + ".I_to_L"
# ~ convert_I_to_L(input_file, output_file)
# ================================================================================

# Run correct_by_consensus()
# This script corrects de novo sequencing errors by grouping predicted sequences of the same mass together and voting the consensus sequence.
# ======================= UNCOMMENT and RUN ======================================
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L"
# ~ output_file = input_file + ".consensus"
# ~ correct_by_consensus(input_file, output_file)
# ================================================================================

# Run filter_by_minlen()
# This script filters out sequences of length less than 5 amino acids.
# ======================= UNCOMMENT and RUN ======================================
# ~ minlen = 5
# ~ input_file = data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus"
# ~ output_file = input_file + ".minlen" + str(minlen)
# ~ filter_by_minlen(input_file, output_file, minlen)
# ================================================================================
# The numbers of features will be reported as:
#   "total_feature =  233589"
#   "minlen_feature =  223507"
#   "removed_feature =  10082"

# Up to this step, we get the following file: 
#   "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"
# We test its accuracy against the test set:
# Run DeepNovo testing
# ======================= UNCOMMENT and RUN ======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test"]
# ~ command += ["--target_file", data_training_dir + "feature.csv.labeled.mass_corrected.test.noshare"]
# ~ command += ["--predicted_file", data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================
# We get these results:
#   "precision_AA_mass_db  = 0.9530"
#   "precision_peptide_mass_db  = 0.8441"

# Repeat the same testing but now against all labeled features:
# Run DeepNovo testing
# ====================== UNCOMMENT and RUN =======================================
# ~ command = ["LD_PRELOAD=\"/usr/lib/libtcmalloc.so\" /usr/bin/time -v python deepnovo_main.py --test"]
# ~ command += ["--target_file", data_training_dir + "feature.csv.labeled.mass_corrected"]
# ~ command += ["--predicted_file", data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"]
# ~ command = " ".join(command)
# ~ print(command)
# ~ os.system(command)
# ================================================================================
# We get these results:
#   "precision_AA_mass_db  = 0.9797"
#   "precision_peptide_mass_db  = 0.9371"
# Note that these accuracy results look better than those against the test set because the test set was not used for training the model.
# The number of de novo only features is also reported as
#   "predicted_only: 68721"
# and they are written to the file 
#   "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"




# ================================================================================
# Step 4.2: Run second round of PEAKS X DB search against the list of database and de novo peptides. 
# ================================================================================

# Before running PEAKS, we need to combine database and de novo peptides into a list.
# This script will select unique de novo peptides, filter out those that belong to the human Swiss-Prot protein database, and combine the remaining de novo peptides and the database peptides identified from Step 1 into a fasta file.
# ======================= UNCOMMENT and RUN ======================================
# ~ aa_workflow_step_4_2.preprocess(
    # ~ denovo_file=data_training_dir + "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only",
    # ~ db_fasta_file=data_fasta_dir + "uniprot_sprot.human.plus_contaminants.fasta",
    # ~ labeled_feature_file=data_training_dir + "feature.csv.labeled.mass_corrected",
    # ~ peptide_list_fasta=data_training_dir + "aa_workflow.step_4.peptide_list.fasta")
# ================================================================================
# The numbers of de novo and database peptides are reported as following:
#   "Number of top-scoring denovo peptides: 17318"
#   "num_db_peptides = 25274"
#   "num_denovo_peptides = 6444" (not in database)

# Run PEAKS X DB search with as following:
#   Select the DENOVO node result from Step 1.1, and select PEAKS DB search;
#   Select option "No digestion" for "Digest mode";
#   Select the fasta file "aa_workflow.step_4.peptide_list.fasta" as the only database, no contaminant;
#   Leave other settings the same as in Step 1.1.
# Set FDR 1.0% and export the "DB search psm.csv" file, rename it to "aa_workflow.step_4.psm.csv".

# Extract de novo peptides from the PSMs of PEAKS X DB search round 2.
# ======================= UNCOMMENT and RUN ======================================
# ~ aa_workflow_step_4_2.postprocess(
    # ~ psm_file = data_training_dir + "aa_workflow.step_4.psm.csv",
    # ~ output_denovo_peptide_file = data_training_dir + "aa_workflow.step_4.output_peptide_list")
# ================================================================================
# The number of de novo peptides is reported as following:
#   "num_denovo_peptides = 1259"




# ================================================================================
# Step 5: Neoantigen selection. 
# ================================================================================
# ~ aa_workflow_step_5.step_5(
    # ~ psm_file=data_training_dir + "aa_workflow.step_4.psm.csv",
    # ~ netmhc_file=data_training_dir + "aa_workflow.step_5.netmhcpan.csv",
    # ~ db_fasta_file=data_fasta_dir + "uniprot_sprot.human.plus_contaminants.fasta",
    # ~ labeled_feature_file=data_training_dir + "feature.csv.labeled",
    # ~ snp_file=data_training_dir + "aa_workflow.step_5.supp_data5_snp.csv",
    # ~ snp_enst_fasta=data_training_dir + "aa_workflow.step_5.supp_data5_snp_enst.fasta",
    # ~ snp_sample_id=patient_id,
    # ~ output_neoantigen_criteria=data_training_dir + "aa_workflow.step_5.output_neoantigen_criteria.csv",
    # ~ output_protein_mutation=data_training_dir + "aa_workflow.step_5.protein_mutation.csv")














