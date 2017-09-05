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
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_lp_rule
# Set the random number seed in the Metronamica file.
from set_rand import set_rand
# Run the Metronamica model.
from run_metro import run_metro
# Calculate the Kappa metric
from kappa import kappa
# Calculate the Kappa Simulation metric
from kappa import ksim
# Calculate the area-weighted clumpiness error.
from area_weighted_clu import area_weighted_clu_error
# Interact with csv files.
import csv


# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.

# base_path = "C:\\Users\\charl\OneDrive\\Documents\\ENC_Py3_1\\"
base_path = "C:\\Users\\a1210607\\ENC_Py3_1\\"

# Set the case study
case_study = "Berlin"
# Set the paths to the directories and relevant data
data_path = base_path + "EU_data\\"
output_path = base_path + "EU_output\\"
# Specify the original map (data at time slice 1).
omap_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
# Specify the actual map (data at time slice 2).
amap_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
# Specify the masking map.
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

# working_directory = ("C:\\Geonamica\\Metronamica\\" + case_study)
working_directory = "C:\\Users\\a1210607\\Geonamica\\Metronamica\\" + case_study

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

# Set the base random number seed.
base_seed = 1000
# Set the number of simulation runs per iteration.
max_runs = 1
# Specify the minimum, maximum, and number of intervals to be tested for the
# meta-parameter values (theta).
min_theta_st = 0.075
max_theta_st = 0.100
min_theta_cp = 0.075
max_theta_cp = 0.100
min_theta_it = 0.0125
max_theta_it = 0.025
intervals = 5

# Generate a list of meta-parameter values to test.
theta_st_values = []
theta_cp_values = []
theta_it_values = []
for i in range(0, intervals):
    st_range = i * (max_theta_st - min_theta_st) / (intervals - 1)
    theta_st_values.append(min_theta_st + st_range)
    cp_range = i * (max_theta_cp - min_theta_cp) / (intervals - 1)
    theta_cp_values.append(min_theta_cp + cp_range)
    it_range = i * (max_theta_it - min_theta_it) / (intervals - 1)
    theta_it_values.append(min_theta_it + it_range)

# Initialise a dictionary to store metric values.
coarse_metrics = {}

# Initialise the map comparison library for the calculation of Fuzzy Kappa.
"""
fuzzy_coefficients = ("C:\\Users\\charl\\OneDrive\\Documents\\ENC_Py3_1\\"
                      "coeff.txt")
analysis_id = mcl.createAnalysis()
mcl.loadMapActual(analysis_id, amap_path)
mcl.loadMaskingMap(analysis_id, mask_path)
mcl.loadFuzzyWeights(analysis_id, fuzzy_coefficients)
"""

# Perform the iterative testing.
for p in range(0, intervals):
    # Allocate the corresponding meta-parameter values.
    theta_st = theta_st_values[p]
    print("Theta self-influence tail value: " + str(theta_st))
    for q in range(0, intervals):
        theta_cp = theta_cp_values[q]
        print("Theta conversion point value: " + str(theta_cp))
        for r in range(0, intervals):
            theta_it = theta_it_values[r]
            print("Theta interactive tail value: " + str(theta_it))
            # Specify a tuple as the coarse_metrics key
            coarse_metrics_key = (theta_st, theta_cp, theta_it)
            coarse_metrics[coarse_metrics_key] = []
            # Set influence values at distance 1 for self-influence tails.
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
            # Set the values for inertia and conversion.
            for i in range(0, act):
                for j in range(0, luc):
                    # Specify the neighbourhood rule key.
                    key = (
                        "from " + luc_names[j] + " to " + luc_names[i + pas]
                    )
                    # If a self-influence rule, set the inertia value.
                    if i + pas == j:
                        if (cont_table[i + pas, luc] >
                                cont_table[luc, i + pas]):
                            rules[key][0] = low_inertia
                        else:
                            inertia_rate = ic_rates[j, i + pas]
                            if inertia_rate > high_inertia_band:
                                rules[key][0] = high_inertia
                            elif inertia_rate > mid_inertia_band:
                                rules[key][0] = mid_inertia
                            else:
                                rules[key][0] = low_inertia
                    # If an interactive rule, set the conversion rule.
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
            # Generate the simulated output and record the results.
            kappa_log = [0] * max_runs
            kappa_sim_log = [0] * max_runs
            clu_log = [0] * max_runs
            # Reset the run count to zero.
            run_count = 0
            while run_count < max_runs:
                print("Run: " + str(run_count))
                # Generate seed, run model to generate output.
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                run_metro(project_file, log_file, working_directory,
                          geo_cmd)
                # Read in the map.
                smap = read_map(smap_path)
                # Calculate the corresponding metrics
                kappa_log[run_count] = kappa(amap, smap, mask)
                kappa_sim_log[run_count] = ksim(omap, amap, smap, mask)
                clu_log[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                             luc, pas, act,
                                                             luc_count)
                # Add 1 to iterator to avoid infinite loop.
                run_count = run_count + 1
            # Log the output metrics in the dictionary.
            coarse_metrics[coarse_metrics_key].append(sum(kappa_log) /
                                                      len(kappa_log))
            coarse_metrics[coarse_metrics_key].append(sum(kappa_sim_log) /
                                                      len(kappa_sim_log))
            coarse_metrics[coarse_metrics_key].append(sum(clu_log) /
                                                      len(clu_log))

# Write the output metrics to a csv file.
metrics_output_file = (output_path + case_study + 
    "\\Meta_cal_output\\fine_cal_output.csv")
# Generate an empty list to store metric values.
store = [0]*6
# Write to csv file.
with open (metrics_output_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    values = ["theta_st", "theta_cp", "theta_it", "FK", "FKS", "CLU"]
    writer.writerow(values)
    for p in range(0, intervals):
        for q in range(0, intervals):
            for r in range(0, intervals):
                theta_st = theta_st_values[p]
                theta_cp = theta_cp_values[q]
                theta_it = theta_it_values[r]
                store[0] = theta_st
                store[1] = theta_cp
                store[2] = theta_it
                # Specify a tuple as the coarse_metrics key
                coarse_metrics_key = (theta_st, theta_cp, theta_it)
                store[3] = coarse_metrics[coarse_metrics_key][0]
                store[4] = coarse_metrics[coarse_metrics_key][1]
                store[5] = coarse_metrics[coarse_metrics_key][2]
                writer.writerow(store)

# Indicate completion with a beep
import winsound
Freq = 2500
Dur = 1000
winsound.Beep(Freq, Dur)
# Finished!
