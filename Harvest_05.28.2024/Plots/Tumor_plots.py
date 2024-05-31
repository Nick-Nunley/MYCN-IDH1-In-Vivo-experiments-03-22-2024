import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import TumorMass as tm

def make_condition_helper(input_value, string_to_strip = '-M\d+-\d+'):
    value_to_return = re.split(string_to_strip, input_value)[0]
    return value_to_return

def calculate_volume_helper(input_value):
    components = str(input_value).split('x')
    # Convert the components to floats
    float_components = [float(component) for component in components]
    # Calculate the product of the components
    product = 1.0
    for number in float_components:
        product *= number
    return product

raw_data = pd.read_csv(filepath_or_buffer = 'Harvest_05.28.2024/Whole_tumor_weights_and_volumes.csv')

# Removing NaN values
raw_data = raw_data.dropna(subset = ['Weight (g)'])

# Making condition columns
raw_data['Condition'] = [make_condition_helper(value) for value in list(raw_data.loc[:, 'Sample-ID'])]
raw_data['Mouse_condition'] = [make_condition_helper(value, '-\d+$') for value in list(raw_data.loc[:, 'Sample-ID'])]

# Converting the caliper measurement to volume
raw_data['Volume (mm^3)'] = [calculate_volume_helper(value) for value in list(raw_data.loc[:, 'Caliper dims (mm^3)'])]

tumor_plots = tm.TumorMass(raw_data)

# Making a plot for each genotypic condition
tumor_plots.plot_data(
    xvals_key = 'Condition',
    yvals_key = 'Weight (g)',
    y_label = 'Mass (g)',
    color = 'black', 
    s = 10
    )

# Making a plot for different mouse samples
tumor_plots.plot_data(
    xvals_key = 'Mouse_condition',
    yvals_key = 'Weight (g)',
    y_label = 'Mass (g)',
    color = 'black',
    s = 10
    )

# Removing NaN volumes
tumor_plots.tumor_data = tumor_plots.tumor_data.dropna()

# Making a plot for each genotypic condition and volume
tumor_plots.plot_data(
    xvals_key = 'Condition',
    yvals_key = 'Volume (mm^3)',
    y_label = r'$\mathbf{Volume (mm^{3})}$',
    color = 'black',
    s = 10
    )

# Making a plot for each different mouse and volume
tumor_plots.plot_data(
    xvals_key = 'Mouse_condition',
    yvals_key = 'Volume (mm^3)',
    y_label = r'$\mathbf{Volume (mm^{3})}$',
    color = 'black',
    s = 10
    )
