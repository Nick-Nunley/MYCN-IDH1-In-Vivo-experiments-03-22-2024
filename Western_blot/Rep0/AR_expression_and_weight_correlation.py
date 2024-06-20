import re
import argparse
import hashlib
from datetime import date
from typing import Union
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import Quantification as qf


def obtain_series_average(
    input_data: pd.DataFrame, indexer: Union[str, int],
    unique_val: Union[str, int], average_key: Union[str, int]
    ) -> float:
        input_data = list(input_data.loc[input_data[indexer] == unique_val, average_key])
        return np.mean(input_data)

def sample_collapse(
    input_data: pd.DataFrame, condition_key: Union[str, int],
    average_key: Union[str, int], sample_regex: str = '-(E|F)'
    ) -> pd.DataFrame:
        temp_index = hashlib.md5(str(date.today()).encode()).hexdigest()
        input_data[temp_index] = [re.sub(sample_regex, '', x) for x in input_data.loc[:, condition_key]]
        input_data = input_data.set_index(keys = [condition_key])
        input_data = input_data.loc[:, [average_key, temp_index]]
        unique_conditions = list(set(list(input_data.loc[:, temp_index])))
        average_vals = [obtain_series_average(
            input_data = input_data,
            indexer = temp_index,
            unique_val = unique_condition,
            average_key = average_key
            ) for unique_condition in unique_conditions]
        input_data = pd.DataFrame({
            str(condition_key): unique_conditions,
            str(average_key): average_vals
            })
        return input_data

def combine_dataframe_objects(*data_frame_objs, index_col: Union[str, int] = 'Sample-ID') -> pd.DataFrame:
    # Type check for multiple args
    for df in data_frame_objs:
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f'Expected pandas.DataFrame, got {type(df).__name__}')
    # Obtaining intersecting conditions
    common_ids = set(data_frame_objs[0].loc[:, index_col])
    for df in data_frame_objs[1:]:
        common_ids &= set(df.loc[:, index_col])
    common_ids = list(common_ids)
    # Filter each DataFrame to include only common index values
    filtered_dfs = []
    for df in data_frame_objs:
        filtered_df = df[df[index_col].isin(common_ids)]
        filtered_df = filtered_df.sort_values(by = index_col)
        filtered_df = filtered_df.reset_index()
        filtered_dfs.append(filtered_df)
    # Merge DataFrames on common index values
    combined_df = pd.concat(filtered_dfs, axis = 1, join = 'inner')
    return combined_df

def main_plot(
    input_data: pd.DataFrame, xvals_key: Union[str, int], yvals_key: Union[str, int],
    x_label: str = None, y_label: str = None,
    xlimits: Union[tuple[Union[float, int]], list[Union[float, int]]] = None,
    ylimits: Union[tuple[Union[float, int]], list[Union[float, int]]] = None,
    output_stats: bool = True, **kwargs
    ):
        xvalues = list(input_data.loc[:, xvals_key])
        yvalues = list(input_data.loc[:, yvals_key])
        plt.figure(figsize = (5, 4))
        ax = plt.subplot(1, 1, 1)
        ax.scatter(
            x = xvalues,
            y = yvalues,
            **kwargs
            )
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.xaxis.set_tick_params(width = 1.5)
        ax.yaxis.set_tick_params(width = 1.5)
        if x_label is not None:
            plt.xlabel(
                xlabel = x_label,
                fontsize = 'large',
                fontweight = 'bold'
                )
        if y_label is not None:
            plt.ylabel(
                ylabel = y_label,
                fontsize = 'large',
                fontweight = 'bold'
                )
        if xlimits is not None:
            if len(xlimits) != 2:
                raise ValueError(f'{xlimits} should be of length: 2')
            plt.xlim(xlimits)
        if ylimits is not None:
            if len(ylimits) != 2:
                raise ValueError(f'{ylimits} should be of length: 2')
            plt.ylim(ylimits)
        plt.show()
        if output_stats is True:
            pearson_test = pearsonr(
                x = xvalues,
                y = yvalues
                )
            spearman_test = spearmanr(
                a = xvalues,
                b = yvalues
                )
            print(f'Pearson test result: {pearson_test}')
            print(f'Spearman test: {spearman_test}')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f',
        '--filepath',
        default = 'Western_blot/Rep0/ImageJ_data.csv',
        help = 'Path to raw BCA data'
        )
    
    parser.add_argument(
        '-w',
        '--weights_data',
        default = 'Harvest_05.28.2024/Whole_tumor_weights_and_volumes.csv',
        help = 'Path to whole tumor weights and volumes .csv file'
        )

    parser.add_argument(
        '-n',
        '--normalization_key',
        default = 'Actin',
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
    ij_object_relative['Condition'] = [re.sub(' \(1/2\)', '', x) for x in list(ij_object_relative.index)]

    # Collapsing down same '-E' and '-F' samples
    ij_object_relative = sample_collapse(
        input_data = ij_object_relative,
        condition_key = 'Condition',
        average_key = 'AR'
        )

    # Renaming 'Condition' column to 'Sample-ID'
    ij_object_relative = ij_object_relative.rename(columns = {'Condition': 'Sample-ID'})
    
    weight_and_vol = pd.read_csv(args.weights_data)

    plotting_data = combine_dataframe_objects(
        ij_object_relative,
        weight_and_vol
        )

    main_plot(
        input_data = plotting_data,
        xvals_key = 'Weight (g)',
        yvals_key = 'AR'
        )
