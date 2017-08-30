"""
This component of the Empirical Neighbourhood Calibration method conducts the
significance test on the specified data.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Evaluates the composition of different neighbourhoods.
from neighbourhood_evaluator import neighbourhood_evaluator
# Conducts the Mann-Whitney U-test for input neighbourhood composition
# dictionaries.
from MWU_test import mwu_test
# Conducts complex mathematical operations.
import math


# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\OneDrive\\Documents\\ENC_Py3_1\\"
# Set the case study
case_study = "Lisbon"
# Set the paths to the directories and relevant data
data_path = base_path + "EU_data\\"
output_path = base_path + "EU_output\\"
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Seaports",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 5

# Read in the map for time slice 1.
omap = read_map(map1_path)
# Read in the map for time slice 2.
amap = read_map(map2_path)
# Read in the masking map.
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
rows = map_dimensions[0]
cols = map_dimensions[1]

# Determine the distances that will be analysed using the module considered 
# distances.
temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered
cdl = temp[1]
# Determine the maximum neighbourhood size (unit) from considered distances
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
N = []
for c in range(0, max_distance):
    N.append(N_all[c])

# Conduct the significance testing using the module 'neighbourhood_evaluator.'
dummy = neighbourhood_evaluator(
    luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols
)
# Store the requisite output into specific dictionaries
all_cells_baseline = dummy[0]
no_new_cells_ci_baseline = dummy[1]
no_cells_ci_baseline = dummy[2]
transition_dictionary = dummy[3]
ef = dummy[4]
# Log scale the enrichment factor values.
log_ef = np.zeros(shape=(max_distance, luc, luc))
for p in range(0, luc):
    for q in range(0, luc):
        for c in range(0, max_distance):
            # If the enrichment factor is not evaluated a value of 0 is given.
            # Hence, a did not error value of -9999 is used.
            if ef[c, p, q] == 0:
                log_ef[c, p, q] = -9999
            else:
                log_ef[c, p, q] = math.log(ef[c, p, q], 10)
# Specify what the significance limit is (recommended value is 1.96 i.e. 95%).
z_limit = 1.96
# Conduct the MWU test using the module 'mwu_test.'
z_scores = mwu_test(
    max_distance, luc, transition_dictionary, all_cells_baseline, N
)

# Determine which rules have meaningful over-representation at nearby distances.
sig_distances = np.zeros(shape=(act, luc))
c = 1
for i in range(0, act):
    for j in range(0, luc):
        if abs(z_scores[c, i+pas, j]) > z_limit and log_ef[c, i+pas, j] > 0:
            sig_distances[i, j] = sig_distances[i, j] + 1
c = 2
for i in range(0, act):
    for j in range(0, luc):
        if abs(z_scores[c, i + pas, j]) > z_limit and log_ef[c, i+pas, j] > 0:
            sig_distances[i, j] = sig_distances[i, j] + 1
att_rules = np.zeros(shape=(act, luc))
for i in range(0, act):
    for j in range(0, luc):
        if i + pas == j:
            att_rules[i, j] = 1
        elif sig_distances[i, j] > 1:
            att_rules[i, j] = 1
att_rules = np.transpose(att_rules)

# Save the attraction rules to a file.
att_rules_file_path = output_path + case_study + "\\Rules\\att_rules.txt"
np.savetxt(att_rules_file_path, att_rules, fmt="%d")

# Completed!
