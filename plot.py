import re
import math
import numpy as np
import pandas as pd
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from matplotlib_venn import venn2
from matplotlib_venn import venn3
matplotlib.rcParams.update({'font.size': 11})
from scipy import stats


# ~ # TEMP
# ~ file_path = "step_5.output_neoantigen_criteria.xlsx"
# ~ value_list = pd.read_excel(file_path, sheetname='5_targets_152_candidates')['total_abundance'].values
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([value_list], labels=['total_abundance'])
# ~ ax.set_yscale('log')
# ~ ax.set_ylabel('Total abundance of supporting PSMs')
# ~ ax.spines["top"].set_visible(False)
# ~ ax.spines["right"].set_visible(False)
# ~ # GRLAFFLKY
# ~ pyplot.plot([1], [134464000], color='red', marker='o', markersize=6)
# ~ pyplot.savefig("temp.png")


def read_netmhcpan_csv(input_file, num_allele):

  best_nM_list = []
  best_rank_list = []
  with open(input_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      best_nM = min([float(row['nM' + str(x)]) for x in range(1, num_allele+1)])
      best_rank = min([float(row['Rank' + str(x)]) for x in range(1, num_allele+1)])
      best_nM_list.append(best_nM)
      best_rank_list.append(best_rank)
  return best_nM_list, best_rank_list


def draw_figure2_boxplot_netmhcpan():

  num_allele = 4
  denovo_path = "deepnovo.aa.figure_2g.netmhcpan_denovo.csv"
  db_path = "deepnovo.aa.figure_2g.netmhcpan_db.csv"
  iedb_path = "deepnovo.aa.figure_2g.netmhcpan_iedb.csv"
  denovo_nM_list, denovo_rank_list = read_netmhcpan_csv(denovo_path, num_allele)
  db_nM_list, db_rank_list = read_netmhcpan_csv(db_path, num_allele)
  iedb_nM_list, iedb_rank_list = read_netmhcpan_csv(iedb_path, num_allele)

  # ~ colors = ['red', 'dodgerblue', 'lightgrey']
  # ~ nM_list = [denovo_nM_list, db_nM_list, iedb_nM_list]
  # ~ print([len(x) for x in nM_list])
  # ~ fig, ax = pyplot.subplots()
  # ~ nM_plot = pyplot.boxplot(nM_list, labels=['De novo', 'Database', 'IEDB'], patch_artist=True)
  # ~ for patch, color in zip(nM_plot['boxes'], colors):
    # ~ patch.set_facecolor(color)
  # ~ ax.set_yscale('log')
  # ~ ax.set_ylabel('Binding affinity (nM, log-scale)')
  # ~ ax.spines["top"].set_visible(False)
  # ~ ax.spines["right"].set_visible(False)
  # ~ # 500-nM threshold
  # ~ pyplot.plot([0, 6], [500, 500], color='black', linestyle='--', linewidth=1)
  # ~ pyplot.savefig("figure2.boxplot_nM.png")

  # ~ colors = ['red', 'dodgerblue']
  colors = ['red', 'dodgerblue', 'lightgrey']
  # ~ rank_list = [denovo_rank_list, db_rank_list]
  rank_list = [denovo_rank_list, db_rank_list, iedb_rank_list]
  print([len(x) for x in rank_list])
  fig, ax = pyplot.subplots()
  # ~ rank_plot = pyplot.boxplot(rank_list, labels=['De novo', 'Database'], patch_artist=True)
  rank_plot = pyplot.boxplot(rank_list, labels=['De novo', 'Database', 'IEDB'], patch_artist=True)
  for patch, color in zip(rank_plot['boxes'], colors):
    patch.set_facecolor(color)
  ax.set_yscale('log')
  ax.set_ylabel('Binding affinity rank (%, log-scale)')
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)
  # 2% and 0.5% threshold
  pyplot.plot([0, 6], [2, 2], color='black', linestyle='--', linewidth=1)
  pyplot.plot([0, 6], [0.5, 0.5], color='black', linestyle='--', linewidth=1)
  pyplot.savefig("figure2.boxplot_rank.png")

  print("np.log(np.median(denovo_rank_list)) =", np.log(np.median(denovo_rank_list)))
  print("np.log(np.median(db_rank_list)) =", np.log(np.median(db_rank_list)))
  print("np.log(np.median(iedb_rank_list))", np.log(np.median(iedb_rank_list)))
  mannwhitneyu, pvalue = stats.mannwhitneyu(denovo_rank_list, iedb_rank_list)
  print("mannwhitneyu =", mannwhitneyu)
  print("pvalue =", pvalue)

# ~ draw_figure2_boxplot_netmhcpan()


def read_immuno_csv(input_file):

  score_list = []
  with open(input_file, 'r') as input_handle:
    csv_reader = csv.DictReader(input_handle, delimiter=',')
    for row in csv_reader:
      score_list.append(float(row['score']))
  return score_list


def draw_figure2_boxplot_immuno():

  denovo_path = "deepnovo.aa.figure_S5.mel_16.immuno_denovo.csv"
  db_path = "deepnovo.aa.figure_S5.mel_16.immuno_db.csv"
  # ~ iedb_path = "deepnovo.aa.figure_2i.immuno_iedb.csv"
  # ~ model_path = "deepnovo.aa.figure_2i.immuno_model.csv"
  denovo_score_list = read_immuno_csv(denovo_path)
  db_score_list = read_immuno_csv(db_path)
  # ~ iedb_score_list = read_immuno_csv(iedb_path)
  # ~ model_score_list = read_immuno_csv(model_path)

  # ~ colors = ['red', 'dodgerblue', 'lightgrey', 'white']
  colors = ['red', 'dodgerblue']
  # ~ score_list = [denovo_score_list, db_score_list, iedb_score_list, model_score_list]
  score_list = [denovo_score_list, db_score_list]
  print([len(x) for x in score_list])
  fig, ax = pyplot.subplots()
  # ~ score_plot = pyplot.boxplot(score_list, labels=['De novo', 'Database', 'IEDB', 'Calis et al.'], patch_artist=True)
  score_plot = pyplot.boxplot(score_list, labels=['De novo', 'Database'], patch_artist=True)
  for patch, color in zip(score_plot['boxes'], colors):
    patch.set_facecolor(color)
  ax.set_ylabel('Immunogenicity')
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)
  pyplot.plot([0, 6], [0., 0.], color='black', linestyle='--', linewidth=1)
  pyplot.savefig("figure2.boxplot_immuno.png")

  print("np.median(denovo_score_list) =", np.median(denovo_score_list))
  print("np.median(db_score_list) =", np.median(db_score_list))
  # ~ print("np.median(iedb_score_list) =", np.median(iedb_score_list))
  # ~ print("np.median(model_score_list) =", np.median(model_score_list))
  mannwhitneyu, pvalue = stats.mannwhitneyu(denovo_score_list, db_score_list)
  print("mannwhitneyu =", mannwhitneyu)
  print("pvalue =", pvalue)

# ~ draw_figure2_boxplot_immuno()


def draw_figure2_venn():

  file_path = "temp.manuscript/deepnovo.aa.figure_2.step6.xlsx"
  denovo_set = set(pd.read_excel(file_path, sheet_name='denovo_peptide')['denovo_peptide'].values)
  db_set = set(pd.read_excel(file_path, sheet_name='db_peptide')['db_peptide'].values)
  iedb_set = set(pd.read_excel(file_path, sheet_name='iedb_peptide')['iedb_peptide'].values)
  set_labels = ('De novo', 'Database', 'IEDB')
  venn_plot = venn3(subsets=[denovo_set, db_set, iedb_set], set_labels=set_labels)
  venn_plot.get_patch_by_id('100').set_color('red')
  venn_plot.get_patch_by_id('100').set_alpha(0.75)
  venn_plot.get_patch_by_id('010').set_color('skyblue')
  venn_plot.get_patch_by_id('010').set_alpha(1.0)
  venn_plot.get_patch_by_id('001').set_color('grey')
  venn_plot.get_patch_by_id('001').set_alpha(0.5)
  pyplot.savefig("venn3.svg")

# ~ draw_figure2_venn()


def plot_spectrum_array(spectrum_array, figure_name):
  print("plot_spectrum_array()")

  figure = plt.figure(1)
  spectrum_count = spectrum_array.shape[0]
  for index in range(spectrum_count):
    plt.subplot(spectrum_count, 1, index+1)
    plt.plot(spectrum_array[index,:])
    plt.ylim((0.0,1.0))
  plt.show()
  figure.savefig(figure_name)
  plt.close()

#~ plot_spectrum_array(np.load("spectrum_original_forward.npy"), "spectrum_original_forward.pdf")


def read_feature_id(input_file, split_char):

  feature_set = set()
  with open(input_file, 'r') as handle:
    header_line = handle.readline()
    for line in handle:
      line = re.split(split_char, line)
      feature_id = line[0]
      feature_set.add(feature_id)
  return feature_set


# figure2.venn2/3.peaks_deepnovo.png
#~ matplotlib.rcParams.update({'font.size': 16})
#~ peaks_set = read_feature_id("data.training/dia.urine.2018_03_29/testing.feature.csv", ',|\r|\n')
#~ peaks_set2 = read_feature_id("data.training/dia.urine.2018_03_29/peaks.denovo.csv.top.feature_id", '\t|\r|\n')
#~ deepnovo_set = read_feature_id("data.training/dia.urine.2018_03_29/testing.unlabeled.csv.deepnovo_denovo.minlen_5.top", '\t|\r|\n')
#~ set_labels = ("PEAKS DB", "PEAKS Denovo", "DeepNovo-DIA")
#~ set_labels = ("", "", "")
#~ venn_plot = venn3(subsets=[peaks_set, peaks_set2, deepnovo_set], set_labels=set_labels)
#~ for text in venn_plot.set_labels:
  #~ text.set_fontsize(16)
#~ pyplot.savefig("figure2.venn2.peaks_18_20_deepnovo.png")


# figure2.bar.aa/peptide.png
def draw_figure2_bar(y_value, y_label, figure_file):

  fig, ax = pyplot.subplots()
  x_value = range(1, len(y_value)+1)
  bar_10k, bar_5k, bar_2k = pyplot.bar(x_value, y_value, width=0.4, align='center')
  bar_10k.set_facecolor('g')
  bar_10k.set_alpha(0.5)
  bar_5k.set_facecolor('lightskyblue')
  bar_2k.set_facecolor('blue')
  bar_2k.set_alpha(0.7)
  for index, value in zip(x_value, y_value):
    ax.text(index-0.1, value + 3, str(value), color='black')
  ax.set_xticks(x_value)
  ax.set_xticklabels(['Top 10k', 'Top 5k', 'Top 2k'])
  ax.set_xlim([0, 4])
  ax.set_ylim([0, 100])
  ax.set_ylabel(y_label)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')
  pyplot.savefig(figure_file)

#~ denovo_only = [41.7, 22.1, 5.5]
#~ draw_figure2_bar(denovo_only, 'Denovo only peptides on top of database (%)', 
                 #~ 'figure2.bar.denovo_only.png')
#~ aa_accuracy = [76.2, 83.6, 94.2]
#~ draw_figure2_bar(aa_accuracy, 'Amino acid accuracy (%)', 'figure2.bar.aa.png')
#~ peptide_accuracy = [41.4, 53.0, 79.9]
#~ draw_figure2_bar(peptide_accuracy, 'Peptide accuracy (%)', 'figure2.bar.peptide.png')


def read_feature_accuracy(input_file, split_char):

  feature_list = []
  with open(input_file, 'r') as handle:
    header_line = handle.readline()
    for line in handle:
      line = re.split(split_char, line)
      feature = {}
      feature["feature_id"] = line[0]
      feature["feature_area_log10"] = math.log10(max(float(line[1]), 1.0))
      feature["predicted_score"] = float(line[4])
      feature["recall_AA"] = float(line[5])
      feature["predicted_len"] = float(line[6])
      feature_list.append(feature)
  return feature_list


# figure2.accuracy.area.png
def draw_figure2_accuracy_area(accuracy_file):

  feature_list = read_feature_accuracy(accuracy_file, '\t|\r|\n')
  num_features = len(feature_list)
  feature_area_log10 = [f['feature_area_log10'] for f in feature_list]
  #~ x_max = int(max(feature_area_log10))
  x_range = np.arange(3, 11, 1.0)
  x_value = []
  y_accuracy = []
  y_proportion = []
  #~ for x in range(1, x_max+1):
  for x in x_range:
    feature_x = [f for f in feature_list if x-0.5 < f['feature_area_log10'] <= x+0.5]
    recall_AA = sum([f['recall_AA'] for f in feature_x])
    target_len = sum([f['predicted_len'] for f in feature_x])
    if target_len > 0:
      x_value.append(x)
      y_accuracy.append(100*recall_AA/target_len)
      y_proportion.append(100.0 * len(feature_x) / num_features)
  
  fig, left_ax = pyplot.subplots()
  pyplot.bar(x_value, y_proportion, width=0.6, align='center', color='salmon', alpha=0.75)
  left_ax.set_xlabel('Feature abundance (log10 scale)')
  left_ax.set_ylabel('Proportion of features (%)', color='salmon')
  left_ax.tick_params('y', colors='salmon')
  left_ax.set_xlim([0, 12])
  left_ax.set_ylim([0, 100])
  left_ax.spines['top'].set_visible(False)
  left_ax.xaxis.set_ticks_position('bottom')
  for index, value in zip(x_value, y_proportion):
    if value > 0:
      left_ax.text(index-0.2, value + 2, str(round(value,1)), fontsize=12, color='black')

  right_ax = left_ax.twinx()
  right_ax.plot(x_value, y_accuracy, '-o', linewidth=2.0, color='blue')
  right_ax.set_ylabel('Amino acid accuracy (%)', color='blue')
  right_ax.tick_params('y', colors='blue')
  right_ax.set_xlim([0, 12])
  right_ax.set_ylim([0, 100])
  for index, value in zip(x_value, y_accuracy):
    if value > 0:
      right_ax.text(index-0.2, value + 2, str(round(value,1)), fontsize=12, color='black')


  pyplot.savefig('figure2.accuracy.area.png')

#~ accuracy_file = "data.training/dia.pecan.plasma.2018_03_29/testing_plasma.feature.csv.deepnovo_denovo.accuracy"
#~ accuracy_file = "data.training/dia.hla.elife.jurkat_oxford/testing_jurkat_oxford.unlabeled.csv.deepnovo_denovo.accuracy"
#~ accuracy_file = "Supplementary Table S6.txt"
#~ draw_figure2_accuracy_area(accuracy_file)


def get_accuracy_score(accuracy_file):

  feature_list = read_feature_accuracy(accuracy_file, '\t|\r|\n')
  num_value = 10
  step = 100 // num_value
  x_value = [x*step for x in range(1, num_value+1)]
  y_value = []
  # find the accuracy for different cutoff
  for x in x_value:
    # ~ feature_x = [f for f in feature_list if x-(step//2) < 100*math.exp(f['predicted_score']) <= x+(step//2)]
    feature_x = [f for f in feature_list if x <= 100*math.exp(f['predicted_score'])]
    recall_AA = sum([f['recall_AA'] for f in feature_x])
    predicted_len = sum([f['predicted_len'] for f in feature_x])
    if predicted_len > 0:
      y_value.append(100*recall_AA/predicted_len)
    else:
      y_value.append(0)

  return x_value, y_value


# figure2.accuracy.score.png
def draw_figure2_accuracy_score():

  accuracy_file = "data.training/aa.hla.bassani.nature_2016.mel_15.class_1/feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo.accuracy"
  accuracy_file_generic = "data.training/aa.hla.bassani.nature_2016.mel_15.class_1/train.exclude_mel_15/feature.csv.labeled.mass_corrected.test.noshare.deepnovo_denovo.accuracy"
  x_test, y_test = get_accuracy_score(accuracy_file)
  x_test_generic, y_test_generic = get_accuracy_score(accuracy_file_generic)
  fig, ax = pyplot.subplots()
  plot_test, = pyplot.plot(x_test, y_test, '-s', linewidth=1.0, color='red', markeredgecolor='red', alpha=0.75)
  plot_test_generic, = pyplot.plot(x_test_generic, y_test_generic, '--o', linewidth=1.0, color='red', markeredgecolor='red', alpha=0.75)
  plot_cutoff_y, = pyplot.plot([0, 100], [95, 95], '--', linewidth=1.0, color='black', markeredgecolor='black', alpha=0.75)
  # ~ plot_cutoff_x, = pyplot.plot([59.5, 59.5], [0, 100], '--', linewidth=1.0, color='black', markeredgecolor='black', alpha=0.75)
  # ~ plot_cutoff_x, = pyplot.plot([61.9, 61.9], [0, 100], '--', linewidth=1.0, color='black', markeredgecolor='black', alpha=0.75)
  pyplot.legend([plot_test, plot_test_generic, plot_cutoff_y], ['Personalized model', 'Generic model', '95% cutoff'], loc='lower right')
  pyplot.yticks([80, 83.5, 86.7, 90, 95, 100], ['80', '83.5', '86.7', '90', '95', '100'])
  # ~ pyplot.title('DeepNovo confidence score for quality control')
  ax.set_xlabel('De novo confidence score')
  ax.set_xlim([0, 105])
  ax.set_ylim([80, 101])
  ax.set_ylabel('Amino acid accuracy (%)')
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')
  pyplot.savefig('figure2.accuracy.score.png')
  # ~ pyplot.savefig('figure2.accuracy.score.svg')

# ~ draw_figure2_accuracy_score()



#~ db_file = "data.training/dia.pecan.plasma.2018_03_28/testing.feature.csv.deepnovo_denovo"
#~ db_abundance = read_feature_abundance(db_file, '\t|\r|\n')
#~ db_abundance_log10 = np.log10(np.array(db_abundance))
#~ denovo_top_file = "data.training/dia.pecan.plasma.2018_03_28/testing.unlabeled.csv.deepnovo_denovo.len_5.7k"
#~ denovo_abundance = read_feature_abundance(denovo_top_file, '\t|\r|\n')
#~ denovo_abundance_log10 = np.log10(np.array(denovo_abundance))
#~ print(len(db_abundance_log10))
#~ print(len(denovo_abundance_log10))
#~ n, bins, patches = pyplot.hist(db_abundance_log10, 50, facecolor='salmon', alpha=0.5)
#~ n, bins, patches = pyplot.hist(denovo_abundance_log10, 50, facecolor='green', alpha=0.5)
#~ pyplot.xlabel('Feature abundance (log10 scale)')
#~ pyplot.ylabel('Number of features')
#~ pyplot.savefig('figure2.hist.denovo.png')
