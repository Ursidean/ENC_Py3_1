"""
This component of the Empirical Neighbourhood Calibration method conducts the
fine calibration, setting neighbourhood rule parameters based on an
adapted implementation of the Van Vliet et al. (2013) method.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Seed the model with a set of neighbourhood rules.
from seed_rules import seed_rules
# Calculate the enrichment factor.
from enrichment_factor import ef
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_lp_rule
# Set the random number seed in the Metronamica file.
from set_rand import set_rand
# Run the Metronamica model to generate output.
from run_metro import run_metro
# Module for calculation of Fuzzy Kappa.
from fuzzy_kappa import fuzzy_kappa
# Module for calculation of Fuzzy Kappa Simulation.
from fuzzy_kappa import fks
# Module for calculation of Area Weighted Absolute Avg. Clumpiness Error (AAWCE)
from area_weigthed_clu import area_weighted_clu_error
# Import the time module to track the duration of the calibration method.
import time
# Interaction with csv file format.
import csv

# Initialisation.
# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\OneDrive\\Documents\\ENC_Py3_1\\"
# Set the case study.
case_study = "Berlin"
# Set the paths to the directories and relevant data.
data_path = base_path + "EU_data\\"
output_path = base_path + "EU_output\\"
# Specify the original map (data at time slice 1).
omap_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
# Specify the actual map (data at time slice 2).
amap_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
# Specify the masking map.
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Specify the fuzzy weights for the calculation of fuzzy Kappa.
fuzzy_coefficients = data_path + "coeff13.txt"
# Specify the fuzzy transition weights for the calculation of FKS.
fuzzy_trans_coefficients = data_path + "coefficients13.txt"

# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Port area",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the socio-economic group for each land-use class. This is used to set the
# order of calibration.
socio_eco_group = ["Recreational", "Other", "Other", "Other", "Other", 
                   "Residential", "Work", "Recreational", "Recreational", 
                   "Other", "Work", "Work", "Work", "Other", "Other"]
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

# Determine the enrichment factor values for the data.
data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)

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
# Load the attraction rules file.
att_rule_file = output_path + case_study + "\\Rules\\att_rules.txt"
att_rules = np.loadtxt(att_rule_file)

# Generate the rules in one of two ways. Rules can either be generated
# by seeding the model with the specified meta-parameters, or read in
# the rules from a csv file.
seed = True
if seed is True:
    theta_st = 0.075
    theta_cp = 0.050
    theta_it = 0.020
    rules = seed_rules(omap, amap, mask, max_distance, luc_names, luc, act, pas,
                       att_rules, theta_st, theta_cp, theta_it, project_file)
else:
    # If false read the rules from an input file.
    initial_rule_file = (output_path + "\\" + case_study + 
                         "\\Rules\\initial_rules.csv")
    # Initialise a dictionary for storing rule values.
    rules = {}
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            rules[key] = [0] * 4
    # Read inputs from csv file
    with open(initial_rule_file, 'r', newline='') as f:
        readCSV = csv.reader(f)
        next(f)  # This skips the header line
        for row in readCSV:
            i = row[0]
            j = row[1]
            key = "from " + i + " to " + j
            rules[key][0] = float(row[2])
            rules[key][1] = float(row[3])
            rules[key][2] = float(row[4])
            rules[key][3] = float(row[5])

# Initialise an array to track the conversion points and distances.
con_rules = np.zeros(shape=(luc, act))
lim_rules = np.zeros(shape=(luc, act))
# Input the rules into the model, analyse the input points for inclusion.
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        fu_elem = j
        lu_elem = i
        y0 = rules[key][0]
        # Conversion points for feature classes are meaningless, as these 
        # cannot change. 
        if i > (act + pas - 1):
            pass
        elif y0 > 0:
            con_rules[i, j] = 1
        y1 = rules[key][1]
        y2 = rules[key][2]
        xe = rules[key][3]
        if att_rules[i, j] == 1:
            lim_rules[i, j] = xe
        set_lp_rule(project_file, fu_elem, lu_elem, y0, y1, y2, xe)

# Fine tuning.
# Initialisation of variables.
# Initialise a counter to track the number of model simulations performed.
iterations_counter = 0
# Set the base random number seed & max number of runs.
base_seed = 1000
max_runs = 10
# Specify the maximum and minimum bounding values for different neighbourhood
# components, and the golden section search tolerance:
# Inertia point
max_ip = 1000
min_ip = 250
gss_ip_tol = (max_ip - min_ip)/100
# The self-influence tail influence value at distance 1.
max_si = 100
min_si = 0
gss_si_tol = (max_si - min_si)/100
# The conversion point
max_cp = 100
min_cp = 0
gss_cp_tol = (max_cp - min_cp)/100
# The interaction tail at distance 1.
max_ct = 100
min_ct = 0
gss_ct_tol = (max_ct - min_ct)/100
# Initialise two lists to track the rule adjusted by name.
rule_from_tracker = []
rule_to_tracker = []
# Initialise a list to track the point that is fine tuned,
# and the final value taken.
pt_tracker = []
pt_value = []
# Initialisation of metrics.
# Set the weighting values for weighted sum calculations (Fuzzy Kappa, Fuzzy 
# Kappa Simulation, AAWCE).
w1 = 1/3
w2 = 1/3
w3 = 1/3
# Set the transformation metric ranges.
r1 = [0.600, 1.000]
r2 = [0.000, 0.500]
r3 = [0.000, 0.100]
# Initialise a list to track the individual metrics averaged over the number of 
# runs that are performed.
clu_log = []
fk_log = []
fks_log = []
# Initialise a set of temp variables for storing metric values.
clu_temp = [0] * max_runs
fk_temp = [0] * max_runs
fks_temp = [0] * max_runs
# Initialise a list to track the composite metric.
comp_log = []
# Determine the starting point metrics.
# Reset the run count.
run_count = 0
while run_count < max_runs:
    # Set the random seed.
    rseed = base_seed + run_count
    set_rand(project_file, rseed)
    # Run Metronamica to generate output.
    run_metro(project_file, log_file, working_directory, geo_cmd)    
    # Read in the simulated map.
    smap = read_map(smap_path)
    # Calculate the AAWCE.
    clu_temp[run_count] = area_weighted_clu_error(amap, smap, mask, luc, pas,
                                                  act, luc_count)
    # Calculate the run Fuzzy Kappa.
    fk_temp[run_count] = fuzzy_kappa(amap_path, smap_path, fuzzy_coefficients)
    # Calculate the run Fuzzy Kappa Simulation.
    fks_temp[run_count] = fks(omap_path, amap_path, smap_path,
                              fuzzy_trans_coefficients)
    # Add one to iterator (prevents an infinite loop!)
    run_count = run_count + 1
# Find the average over the number of runs performed.
clu_avg = sum(clu_temp) / len(clu_temp)
fk_avg = sum(fk_temp) / len(fk_temp)
fks_avg = sum(fks_temp) / len(fks_temp)
# Now determine the weighted sum.




"""
# Line-search calibration.
# CODE
# Initialise a matrix to track rule adjustment
rule_tracker = np.zeros(shape=(total_iterations, 6))
# Track the start of the calibration.
start = time.time()
"""















"""
# Track the end of the calibration.
end = time.time()
# Determine the duration of the calibration method.
duration = end - start
duration_h = duration/3600

# Record the duration of calibration.
output_duration_file = output_path + "\\" + case_study + "\\duration.txt"
f = open(output_duration_file, "w")
f.write(str(duration_h))
f.close()

# Write the output rules.
output_rules_file = output_path + case_study + "\\Rules\\final_rules.csv"
store = [0]*6
with open(output_rules_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line to the file
    values = ["from", "to", "y0", "y1", "y2", "xe"]
    writer.writerow(values)
    # Now write the neighbourhood rules in the form from ... to ...
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            store[0] = luc_names[i]
            store[1] = luc_names[j + pas]
            store[2] = rules[key][0]
            store[3] = rules[key][1]
            store[4] = rules[key][2]
            store[5] = rules[key][3]
            writer.writerow(store)
# Write the output log.
log_file = (output_path + "\\" + case_study + "\\" + case_study + 
            "_fine_tuning_output.csv")
store = [0]*6
with open(log_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line to the file
    values = ["from", "to", "distance", "new value", "max deviation", "KSIM"]
    writer.writerow(values)
    # Now write the output
    for i in range(0, total_iterations):
        store[0] = luc_names[int(rule_tracker[i, 0])]
        store[1] = luc_names[int(rule_tracker[i, 1] + pas)]
        store[2] = rule_tracker[i, 2]
        store[3] = rule_tracker[i, 3]
        store[4] = rule_tracker[i, 4]
        store[5] = rule_tracker[i, 5]
        writer.writerow(store)


# At conclusion print final Kappa Simulation value.
final_ksim = ksim(omap, amap, smap, mask)
print("Final map kappa simulation value: " + str(final_ksim))

# Indicate completion with a beep.
import winsound
Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 1000 # Set Duration To 1000 ms == 1 second
winsound.Beep(Freq,Dur)

# Completed!
"""