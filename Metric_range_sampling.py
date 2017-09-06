"""
This component of the Empirical Neighbourhood Calibration method conducts the
coarse calibration, setting neighbourhood rule parameters based on a
structured sampling evaluation.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Calculate the enrichment factor.
from enrichment_factor import ef
# Log scale the enrichment factor values (default is to base 10).
from log_scale_ef import log_scale_ef
# Generate a contingency table for the data.
from contingency_table import contingency_table
# Randomly sample values.
import random
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_lp_rule
# Set the random number seed in the Metronamica file.
from set_rand import set_rand
# Run the Metronamica model.
from run_metro import run_metro
# Import the Map Comparison Library for calculation of Fuzzy Kappa and FKS.
import mcl
# Calculate the area-weighted clumpiness error.
from area_weighted_clu import area_weighted_clu_error
# Interact with csv files.
import csv

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\OneDrive\\Documents\\ENC_Py3_1\\"
# Set the case study
case_study = "Madrid"
# Set the paths to the directories and relevant data
data_path = base_path + "EU_data\\"
output_path = base_path + "EU_output\\"
# Specify the original map (data at time slice 1).
omap_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
# Specify the actual map (data at time slice 2).
amap_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
# Specify the masking map.
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Specify the fuzzy weights for the calculation of fuzzy Kappa.
fuzzy_coefficients = data_path + "coeff14.txt"
# Specify the fuzzy transition weights for the calculation of FKS.
fuzzy_trans_coefficients = data_path + "coefficients14.txt"

# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Port area",
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
omap = read_map(omap_path)
# Read in the map for time slice 2.
amap = read_map(amap_path)
# Read in the masking map.
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
rows = map_dimensions[0]
cols = map_dimensions[1]

# Count the presence of each land-use class in the actual map. This is
# used in the calculation of area-weighted average clumpiness across the
# active classes.
luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[amap[i, j]] = luc_count[amap[i, j]] + 1

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

# Set the working directory, which contains the geoproject file.
working_directory = ("C:\\Geonamica\\Metronamica\\" + case_study)
# Set the project file path.
project_file = working_directory + "\\" + case_study + ".geoproj"
# Set the path to the command line version of Geonamica
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
# Set the path to the log file.
log_file = base_path + "LogSettings.xml"
# Set the path to the simulated output map
smap_path = (
    working_directory + "\\Log\\Land_use\\"
                        "Land use map_2000-Jan-01 00_00_00.rst"
)

# Specify the bands levels.
# The inertia band levels, for setting inertia values.
high_inertia_band = 0.95
mid_inertia_band = 0.90
# The conversion band levels, for setting conversion values.
high_conversion_band = 0.5
mid_conversion_band = 0.1
low_conversion_band = 0.025
# The enrichment factor band levels, for setting the tail values.
high_ef = 1.0
mid_ef = 0.5

# Specify high, medium and low inertia values.
# The default settings are a high inertia of 1000, med 500, low 250.
high_inertia = 1000.0
mid_inertia = 500.0
low_inertia = 250.0

# Generate the Enrichment Factor and contingency table.
data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)
log_data_ef = log_scale_ef(data_ef, 10, luc, act, pas, max_distance)
cont_table = contingency_table(omap, amap, mask, luc, rows, cols)

# Determine the rates of inertia and conversion between different land-uses.
ic_rates = np.zeros(shape=(luc, luc))
for i in range(0, luc):
    for j in range(0, luc):
        if i == j:
            if cont_table[i, luc] > 0:
                ic_rates[i, j] = cont_table[i, j] / cont_table[i, luc]
        else:
            conversions = abs(float(cont_table[j, j]) - float(cont_table[luc, j]))
            if conversions > 0:
                ic_rates[i, j] = float(cont_table[i, j]) / float(conversions)
# Load the attraction rules file
att_rule_file = output_path + case_study + "\\Rules\\att_rules.txt"
att_rules = np.loadtxt(att_rule_file)

# Initialise a dictionary to track the rules.
rules = {}
for i in range(0, act):
    for j in range(0, luc):
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        rules[key] = [0, 0, 0, 5]

# Specify the ranges of the different parameter values
theta_st_low = 0.0
theta_st_high = 0.250
theta_cp_low = 0.0
theta_cp_high = 0.250
theta_it_low = 0.0
theta_it_high = 0.050

# Initialise a tracking counter for different permutations
# of the randomly generated parameters.
iteration_counter = 0
max_iterations = 100
# Initialise a tracking array
tracker = np.zeros(shape=(max_iterations, 6))

while iteration_counter < max_iterations:
    # Generate random, uniformly distributed values of theta:
    # Inertia tail, conversion point, conversion tail.
    rand_st = random.uniform(theta_st_low, theta_st_high)
    rand_cp = random.uniform(theta_cp_low, theta_cp_high)
    rand_it = random.uniform(theta_it_low, theta_it_high)
    # Store these values in the tracking array
    tracker[iteration_counter, 0] = rand_st
    tracker[iteration_counter, 1] = rand_cp
    tracker[iteration_counter, 2] = rand_it
    # Calculate the relevant parameter values
    theta_st = rand_st
    theta_cp = rand_cp
    theta_it = rand_it
    # Generate the neighbourhood rules
    # Set the tested parameter values for the self-influence tails.
    d1_high_st_value = high_inertia * theta_st
    d1_mid_st_value = mid_inertia * theta_st
    d1_low_st_value = low_inertia * theta_st
    # Set influence values at distance 2 for self-influence tails.
    d2_high_st_value = high_inertia * theta_st * 0.1
    d2_mid_st_value = mid_inertia * theta_st * 0.1
    d2_low_st_value = low_inertia * theta_st * 0.1
    # Set the conversion parameter values.
    high_conversion = high_inertia * theta_cp
    mid_conversion = mid_inertia * theta_cp
    low_conversion = low_inertia * theta_cp
    # Set influence values at distance 1 for interaction tails.
    d1_high_it_value = high_inertia * theta_it
    d1_mid_it_value = mid_inertia * theta_it
    d1_low_it_value = low_inertia * theta_it
    # Set influence values at distance 2 for interaction tails.
    d2_high_it_value = high_inertia * theta_it * 0.1
    d2_mid_it_value = mid_inertia * theta_it * 0.1
    d2_low_it_value = low_inertia * theta_it * 0.1
    # Input the rules to be analysed.
    # Set the values for inertia and conversion.
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # If a self-influence rule, set the inertia point.
            if i + pas == j:
                if cont_table[i + pas, luc] > cont_table[luc, i + pas]:
                    rules[key][0] = low_inertia
                else:
                    inertia_rate = ic_rates[j, i + pas]
                    if inertia_rate > high_inertia_band:
                        rules[key][0] = high_inertia
                    elif inertia_rate > mid_inertia_band:
                        rules[key][0] = mid_inertia
                    else:
                        rules[key][0] = low_inertia
            # If an interactive rule, set the conversion point.
            else:
                conversion_rate = ic_rates[j, i + pas]
                if conversion_rate > high_conversion_band:
                    rules[key][0] = high_conversion
                elif conversion_rate > mid_conversion_band:
                    rules[key][0] = mid_conversion
                elif conversion_rate > low_conversion_band:
                    rules[key][0] = low_conversion
    # Set the values for self-influence and conversion tails.
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # If a self-influence rule, set the self-influence attraction values.
            if i + pas == j:
                for c in range(1, 3):
                    if c == 1:
                        if att_rules[j, i] == 1:
                            if log_data_ef[c, j, i] > high_ef:
                                rules[key][c] = d1_high_st_value
                            elif log_data_ef[c, j, i] > mid_ef:
                                rules[key][c] = d1_mid_st_value
                            else:
                                rules[key][c] = d1_low_st_value
                    elif c == 2:
                        if att_rules[j, i] == 1:
                            if log_data_ef[c, j, i] > high_ef:
                                rules[key][c] = d2_high_st_value
                            elif log_data_ef[c, j, i] > mid_ef:
                                rules[key][c] = d2_mid_st_value
                            else:
                                rules[key][c] = d2_low_st_value
            # If a conversion rule, set the interactive attraction values.
            else:
                if (
                                    att_rules[j, i] == 1 and log_data_ef[1, j, i] > 0
                        and log_data_ef[2, j, i] > 0
                ):
                    for c in range(1, 3):
                        if c == 1:
                            if log_data_ef[c, j, i] > high_ef:
                                rules[key][c] = d1_high_it_value
                            elif log_data_ef[c, j, i] > mid_ef:
                                rules[key][c] = d1_mid_it_value
                            elif log_data_ef[c, j, i] > 0:
                                rules[key][c] = d1_low_it_value
                        elif c == 2:
                            if log_data_ef[c, j, i] > high_ef:
                                rules[key][c] = d2_high_it_value
                            elif log_data_ef[c, j, i] > mid_ef:
                                rules[key][c] = d2_mid_it_value
                            elif log_data_ef[c, j, i] > 0:
                                rules[key][c] = d2_low_it_value
    # Set the end-points of each attraction rule
    for i in range(0, act):
        for j in range(0, luc):
            if att_rules[j, i] == 0:
                pass
            else:
                # Specify the neighbourhood rule key.
                key = "from " + luc_names[j] + " to " + luc_names[i + pas]
                # Iterate through to find end point
                for c in range(2, 5):
                    if att_rules[j, i] == 1 and log_data_ef[c, j, i] > 0:
                        rules[key][3] = c + 1
    # Input the rules into the model.
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            fu_elem = j
            lu_elem = i
            y0 = rules[key][0]
            y1 = rules[key][1]
            y2 = rules[key][2]
            xe = rules[key][3]
            set_lp_rule(project_file, fu_elem, lu_elem, y0, y1, y2, xe)
    # Generate the simulated output and track the results.
    it_fk = 0
    it_fks = 0
    it_clu = 0
    rseed = 1000
    set_rand(project_file, rseed)
    run_metro(project_file, log_file, working_directory, geo_cmd)
    # Create the analysis id for the Fuzzy Kappa.
    analysis_id_fk = mcl.createAnalysis()
    # Create the analysis id for the FKS.
    analysis_id_fks = mcl.createAnalysis()
    # Load the original map for the analysis of FK.
    mcl.loadOriginalMap(analysis_id_fks, omap_path)
    # Load the actual map for the analysis of FK and FKS.
    mcl.loadMapActual(analysis_id_fk, amap_path)
    mcl.loadMapActual(analysis_id_fks, amap_path)
    # Load the fuzzy weights for Fuzzy Kappa.
    mcl.loadFuzzyWeights(analysis_id_fk, fuzzy_coefficients)
    # Load the fuzzy transition weights for FKS.
    mcl.loadTransitionFuzzyWeights(analysis_id_fks, fuzzy_trans_coefficients)
    # Read in the simulated map.
    smap = read_map(smap_path)
    # Read in the map for analysis of fuzzy kappa.
    mcl.loadMapSimulated(analysis_id_fk, smap_path)
    # Read in the map for analysis of fuzzy Kappa Simulation.
    mcl.loadMapSimulated(analysis_id_fks, smap_path)
    # Calculate the corresponding metrics
    it_fk = mcl.getFuzzyKappa(analysis_id_fk)
    it_fks = mcl.getFuzzyKappaSim(analysis_id_fks)
    it_clu = area_weighted_clu_error(amap, smap, mask, luc, pas, act,
                                     luc_count)
    # Clear the analysis for the fuzzy kappa (prevents errors).
    mcl.clear(analysis_id_fk)
    # Clear the analysis for the FKS (prevents errors).
    mcl.clear(analysis_id_fks)
    # Store the values corresponding to the metric values generated.
    tracker[iteration_counter, 3] = it_fk
    tracker[iteration_counter, 4] = it_fks
    tracker[iteration_counter, 5] = it_clu
    # Add one to the iterations counter
    iteration_counter = iteration_counter + 1
    # Provide user feedback
    print("Interactions completed: " + str(iteration_counter))

# Write the output
output_file = (output_path + case_study + "\\Metric_ranges\\metric_ranges.csv")
store = [0] * 6
with open (output_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    values = ["Theta_st", "Theta_cp", "Theta_it", "FK", "FKS", "CLU"]
    writer.writerow(values)
    for i in range(0, max_iterations):
        store[0] = tracker[i, 0]
        store[1] = tracker[i, 1]
        store[2] = tracker[i, 2]
        store[3] = tracker[i, 3]
        store[4] = tracker[i, 4]
        store[5] = tracker[i, 5]
        writer.writerow(store)

# Finished!
