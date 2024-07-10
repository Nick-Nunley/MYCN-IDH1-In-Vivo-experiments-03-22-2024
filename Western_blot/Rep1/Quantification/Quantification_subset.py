import re
import argparse
import numpy as np
import pandas as pd
import Quantification as qf

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f',
        '--filepath',
        default = 'Western_blot/Rep1/Quantification/ImageJ_data_subset.csv',
        help = 'Path to raw BCA data'
        )

    parser.add_argument(
        '-n',
        '--normalization_key',
        default = 'Vinculin',
        help = 'The column index for the normalization values'
        )

    parser.add_argument(
        '-b',
        '--blank_key',
        default = 'Blank',
        help = 'The row index value for the blank values'
        )
    
    parser.add_argument(
        '-c',
        '--control_indecies',
        default = [4, 5, 6, 7],
        help = 'An int or list of integers with the control indecies'
        )

    args = parser.parse_args()

    if not (isinstance(args.control_indecies, int) or isinstance(args.control_indecies, list)):
        raise TypeError(f'{args.control_indecies} is not type \'list\' or \'int\'')
    if isinstance(args.control_indecies, list) and not all(isinstance(x, int) for x in args.control_indecies):
        raise TypeError(f'{args.control_indecies} must contain integers')
    
    np.random.seed(1717)

    input_data = pd.read_csv(
        filepath_or_buffer = args.filepath,
        sep = ','
        )
    input_data = input_data.set_index(keys = ['Sample'])

    ij_object = qf.IJ_data(
        data = input_data,
        normalization_data = list(input_data.loc[input_data.index != args.blank_key, args.normalization_key]),
        blank_data = list(input_data.loc[args.blank_key, :]),
        norm_blank_index = input_data.columns.get_loc(args.normalization_key),
        control_indecies = args.control_indecies
        )
    
    ij_object_relative = ij_object.main_method(blank_key = args.blank_key)
    ij_object_relative['Condition'] = [re.sub('-M[0-9]-[0-9]-(E|F)', '', re.sub(' \(1/2\)', '', x)) for x in list(ij_object_relative.index)]
    ij_object_relative['Color'] = (['orange'] * 4) + (['purple'] * 4)

    ij_object.feature_labels = [x for x in list(ij_object_relative.columns) if x != 'Condition' and x != 'Color' and x != args.normalization_key]

    for feature in ij_object.feature_labels:
        qf.barplot_data(
            input_data = ij_object_relative,
            xvals_key = 'Condition',
            yvals_key = feature,
            color_index = 'Color',
            y_label = f'{feature} Intesity / {args.normalization_key} Intensity',
            ylimits = [0, 4]
            )
