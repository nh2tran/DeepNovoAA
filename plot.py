import re
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from matplotlib_venn import venn2
from matplotlib_venn import venn3
matplotlib.rcParams.update({'font.size': 15})


# ~ file_path = "testing_plasma.venn.xlsx"
# ~ #~ excel = pd.ExcelFile(file_path)
# ~ pecan_set = set(pd.read_excel(file_path, sheetname='pecan_unique')['pecan_unique'].values)
# ~ openswath_set = set(pd.read_excel(file_path, sheetname='openswath_unique')['openswath_unique'].values)
# ~ spectronaut_set = set(pd.read_excel(file_path, sheetname='spectronaut_unique')['spectronaut_unique'].values)
# ~ #~ inhouse_set = set(pd.read_excel(file_path, sheetname='inhouse_unique')['inhouse_unique'].values)
# ~ top90_set = set(pd.read_excel(file_path, sheetname='top90_final')['top90_final'].values)
# ~ print('len(pecan_set) =', len(pecan_set))
# ~ print('len(spectronaut_set) =', len(spectronaut_set))
# ~ #~ print('len(inhouse_set) =', len(inhouse_set))
# ~ print('len(top90_set) =', len(top90_set))

# ~ set_labels = ('DeepNovo', 'PECAN', 'Spectronaut')
# ~ venn_plot = venn3(subsets=[top90_set, pecan_set, spectronaut_set], set_labels=set_labels)
# ~ pyplot.savefig("venn3.svg")



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
    feature_x = [f for f in feature_list if x-(step//2) < 100*math.exp(f['predicted_score']) <= x+(step//2)]
    recall_AA = sum([f['recall_AA'] for f in feature_x])
    predicted_len = sum([f['predicted_len'] for f in feature_x])
    if predicted_len > 0:
      y_value.append(100*recall_AA/predicted_len)
    else:
      y_value.append(0)

  return x_value, y_value


# figure2.accuracy.score.png
def draw_figure2_accuracy_score():

  accuracy_file_oc = "data.training/bassani.hla.2018_10_18.correct_mass_shift/identified_features.csv.valid.nodup.deepnovo_denovo.accuracy"
  accuracy_file_uti = "data.training/bassani.hla.2018_10_18.correct_mass_shift/identified_features.csv.valid.nodup.deepnovo_denovo.accuracy"
  accuracy_file_plasma = "data.training/bassani.hla.2018_10_18.correct_mass_shift/identified_features.csv.valid.nodup.deepnovo_denovo.accuracy"
  # ~ accuracy_file_oc = "data.training/dia.hla.elife.jurkat_oxford/testing_jurkat_oxford.unlabeled.csv.deepnovo_denovo.accuracy"
  # ~ accuracy_file_uti = "data.training/dia.hla.elife.jurkat_oxford/testing_jurkat_oxford.unlabeled.csv.deepnovo_denovo.accuracy"
  # ~ accuracy_file_plasma = "data.training/dia.hla.elife.jurkat_oxford/testing_jurkat_oxford.unlabeled.csv.deepnovo_denovo.accuracy"
  x_oc, y_oc = get_accuracy_score(accuracy_file_oc)
  x_uti, y_uti = get_accuracy_score(accuracy_file_uti)
  x_plasma, y_plasma = get_accuracy_score(accuracy_file_plasma)
  fig, ax = pyplot.subplots()
  plot_oc, = pyplot.plot(x_oc, y_oc, '-o', linewidth=2.0, color='chocolate', markeredgecolor='chocolate', alpha=0.75)
  plot_uti, = pyplot.plot(x_uti, y_uti, '-s', linewidth=2.0, color='orange', markeredgecolor='orange', alpha=0.75)
  plot_plasma, = pyplot.plot(x_plasma, y_plasma, '-^', linewidth=2.0, color='salmon', markeredgecolor='salmon', alpha=0.75)
  pyplot.legend([plot_oc, plot_uti, plot_plasma], ['OC', 'UTI', 'Plasma'], loc='lower right')
  pyplot.title('DeepNovo confidence score for quality control')
  ax.set_xlabel('De novo confidence score')
  ax.set_xlim([0, 105])
  ax.set_ylim([0, 105])
  ax.set_ylabel('Amino acid accuracy (%)')
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')
  pyplot.savefig('figure2.accuracy.score.png')
  # ~ pyplot.savefig('figure2.accuracy.score.svg')

draw_figure2_accuracy_score()



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
