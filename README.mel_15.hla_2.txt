================================================================================
Workflow of neoantigen discovery by personalized de novo sequencing.
================================================================================

Step-by-step instructions based on the following example dataset:

  Patient Mel-15 (Bassani-Sternberg et al., Nature Communication, 2016)
  HLA class 2: 8 raw files

    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_1.raw
    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_1_1.raw
    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_2.raw
    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_2_1.raw
    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_3.raw
    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_3_1.raw
    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_4.raw
    20141220_QEp7_MiBa_SA_HLA-II-p_MM15_4_1.raw




================================================================================
Step 1: Build the immunopeptidome of the patient
================================================================================

  This step 1 took about 7 hours on a laptop with 4 CPU cores i7, 16 GB memory

  Step 1.1: Run PEAKS X DB search on the raw files with the following parameters:

    Enzyme: None
    Instrument: Orbi-Orbi
    Fragment: HCD
    Acquisition: DDA

    Parent Mass Error Tolerance: 15.0 ppm
    Fragment Mass Error Tolerance: 0.05 Da
    Precursor Mass Search Type: monoisotopic
    Enzyme: None
    Digest Mode: Unspecific
    Max Missed Cleavages: 100
    Variable Modifications:
      Oxidation (M): 15.99
      Deamidation (NQ): 0.98
    Max Variable PTM Per Peptide: 3
    Database: uniprot_sprot.human
    Taxon: All
    Contaminant Database: contaminants_maxquant
    Searched Entry: 20488
    FDR Estimation: Enabled
    Merge Options: no merge
    Precursor Options: corrected
    Charge Options: no correction
    Filter Options: no filter
    Process: true
    Associate chimera: no

  Step 1.2: Set FDR 1.0%

    The corresponding number of peptide-spectrum matches (PSMs) will be about "67021", the number of unique peptides is be about "9664".

  Step 1.3: Right-click on the DB search node "PEAKS_10", select option "Deep Denovo Export".

    You will get the following 8 pairs of csv and mgf files in your PEAKS project folder:

      export_0.csv, export_0.mgf
      export_1.csv, export_1.mgf
      export_2.csv, export_2.mgf
      export_3.csv, export_3.mgf
      export_4.csv, export_4.mgf
      export_5.csv, export_5.mgf
      export_6.csv, export_6.mgf
      export_7.csv, export_7.mgf




================================================================================
Step 2: Train personalized DeepNovo model
================================================================================

  This step 2 took about 12 hours on a server with 32 CPU cores Xeon, 32 GB memory

  Note that you will need to specify the paths to your own data and model folders when you run the Python scripts. The provided scripts just show examples of my data and model folders.

  Step 2.1: Prepare the training data. The scripts are provided in "deepnovo_preprocess.py", uncomment and run them in the following order:

    Run merge_mgf_file() and merge_feature_file()

      You need to specify the path to your folder of csv and mgf files.
      You will get two output files: "spectrum.mgf" and "feature.csv".
      Both functions also report the number of entries that have been processed: "counter = 202511".
      That number should be the same as the total number of MS/MS spectra from the raw files.

    Run split_feature_unlabel()
  
      It will split the "feature.csv" into 2 files: "feature.csv.labeled" and "feature.csv.unlabeled".
      It also reports the number of labeled and unlabel features: "num_labeled = 67021" and "num_unlabeled = 135490".
      Note that "67021" is also the number of PSMs reported at FDR 1.0% in Step 1.

    Run calculate_mass_shift_ppm() and correct_mass_shift_ppm()

      The mass shift is calculated from labeled features: "mean_precursor_ppm = 2.098829".
      Then mass is corrected for all features: "feature.csv.labeled.mass_corrected" and "feature.csv.mass_corrected".

    Run split_feature_training_noshare()

      It will split the labeled features into train/valid/test sets with "proportion = [0.9, 0.05, 0.05]".
        "feature.csv.labeled.mass_corrected.train.noshare"
        "feature.csv.labeled.mass_corrected.valid.noshare"
        "feature.csv.labeled.mass_corrected.test.noshare"
      Those 3 sets do not share common peptides.
      Their numbers of features are also reported.

  Step 2.2: Training DeepNovo model.

    Copy a pre-trained model from "train.mel_15.hla_1" into a new folder "train.mel_15.hla_2.pretrained.mel_15.hla_1".

    (You can also train from scratch. In that case, just create an empty folder)

    Specify train/valid/test files in "deepnovo_config.py" in the following lines:

      input_spectrum_file_train = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/spectrum.mgf"
      input_feature_file_train = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected.train.noshare"
      input_spectrum_file_valid = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/spectrum.mgf"
      input_feature_file_valid = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected.valid.noshare"
      input_spectrum_file_test = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/spectrum.mgf"
      input_feature_file_test = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected.test.noshare"

    Run the following command to start re-training:

      python deepnovo_main.py --train_dir train.mel_15.hla_2.pretrained.mel_15.hla_1 --train --reset_step

    The re-training will stop after 10 epoch. The model with best performance on the valid set is recorded in the model folder "train.mel_15.hla_2.pretrained.mel_15.hla_1".

    After training, specify the following test files in "deepnovo_config.py" for testing:

      denovo_input_spectrum_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/spectrum.mgf"
      denovo_input_feature_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected.test.noshare"

      target_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected.test.noshare"
      predicted_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo"
    
    Run the following commands for testing:

      python deepnovo_main.py --train_dir train.mel_15.hla_2.pretrained.mel_15.hla_1 --search_denovo
      python deepnovo_main.py --test

    The testing accuracy at the amino acid (AA) and peptide levels will be reported in these two lines:

      "precision_AA_mass_db  = 0.6986"
      "precision_peptide_mass_db  = 0.6029"




================================================================================
Step 3: Perform personalized de novo sequencing with DeepNovo
================================================================================

  This step 3 took about 6 hours on a server with 32 CPU cores Xeon, 32 GB memory

  Specify these files in "deepnovo_config.py" to perform de novo sequencing on all features (label and unlabeled):

    denovo_input_spectrum_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/spectrum.mgf"
    denovo_input_feature_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.mass_corrected"

  Run the following command for de novo sequencing:

    python deepnovo_main.py --train_dir train.mel_15.hla_2.pretrained.mel_15.hla_1 --search_denovo

  The de novo results will be written to the file "feature.csv.mass_corrected.deepnovo_denovo"

  Note that the tool will report the total number of features "202511", the number of processed features "160698", and the number of skipped features "41813" (since their masses exceed 3000 Dalton) 




================================================================================
Step 4: Quality control
================================================================================

  Step 4.1: Post-process de novo results to improve their accuracy. The scripts are provided in "deepnovo_postprocess.py", uncomment and run them in the following order:

    Run select_top_score()

      This script selects a threshold of de novo confidence scores and uses it to filter de novo results.
      The score threshold is calculated based on a 95% cutoff of the testing accuracy obtained at the end of Step 2 above.
      After this step you'll get the file "feature.csv.mass_corrected.deepnovo_denovo.top95".
      The score cutoff and the number of selected features will also be reported:

        score_cutoff =  -0.54
        select_feature =  64832

    Run convert_I_to_L()

      This script converts I (Isoleucine) to L (Leucine) in all de novo peptides, because de novo sequencing is not able to distinguish them.

    Run correct_by_consensus()

      This script corrects some de novo errors by grouping predicted sequences of the same mass together and voting the consensus sequence.

    Run filter_by_minlen()

      This script filters out sequences of length less than 5 amino acids.

    Up to this step, you will get the following file: 

      "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"

    We test sequencing accuracy against the test set by specifying the following in "deepnovo_config.py":

      target_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected.test.noshare"
      predicted_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5"

    Then run

      python deepnovo_main.py --test

    and get these results

      "precision_AA_mass_db  = 0.9460"
      "precision_peptide_mass_db  = 0.9033"

    Repeat the same testing but now against all labeled features specifying the following in "deepnovo_config.py":

      target_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled.mass_corrected"

    You will get these results

      "precision_AA_mass_db  = 0.9746"
      "precision_peptide_mass_db  = 0.9624"

    Note that these accuracy results look better than those against the test set because the test set was not used for training the model.

    The number of de novo only features is also reported as

      "predicted_only: 8564"

    and they are written to the file 

      "feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"
      
    The number of unique peptides is 3913 (with modifications: Oxidation, Deamidation) and 3867 without.

  Step 4.2: Run second round of PEAKS X DB search against the list of database and de novo peptides.

    Specify the following input files in the script "step_4_2.py" and run it:

      denovo_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.mass_corrected.deepnovo_denovo.top95.I_to_L.consensus.minlen5.denovo_only"
      db_fasta_file = "data/uniprot_sprot.human.plus_contaminants.fasta"
      labeled_feature_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/feature.csv.labeled"
      peptide_list_fasta = "data.training/aa.hla.bassani.nature_2016.mel_15.class_2/step_4.peptide_list.fasta"

    The script will select unique de novo peptides, filter out those that belong to the human Swiss-Prot protein database, and combine the remaining de novo peptides and the database peptides identified from Step 1 into a fasta file "step_4.peptide_list.fasta".

      "num_db_peptides = 9664"
      "num_denovo_peptides = 2717"

    Run PEAKS X DB search with the following:

      Select the DENOVO node result from Step 1.1, and select PEAKS DB search;
      Select option "No digestion" for "Digest mode";
      Select the fasta file "step_4.peptide_list.fasta" as the only database, no contaminant;
      Leave other settings the same as in Step 1.1.

    Set FDR 1.0% and export the "DB search psm.csv" file, rename it to "step_4.DB search psm round_2_FDR_1%.csv".

    Run "step_4_2.postprocess.py" to extract de novo peptides into the file "step_4.output_peptide_list".
    

