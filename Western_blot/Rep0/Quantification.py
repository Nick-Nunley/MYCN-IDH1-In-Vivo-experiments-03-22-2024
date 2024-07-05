import re
import argparse
from typing import Union
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def subtract_blank(sample: Union[float , int], blank: Union[float, int]) -> Union[float, int]:
    transformed_value = sample - blank
    return transformed_value

def norm_division(sample: float, normalization_factor: float) -> float:
    transformed_value = sample / normalization_factor
    return transformed_value

def run_anova(input_data, xvals_key, yvals_key, post_hoc_test = True):
    unique_x_vals = list(set(list(input_data.loc[:, xvals_key])))
    grouped_vals = [list(input_data.loc[input_data[xvals_key] == unique_value, yvals_key]) for unique_value in unique_x_vals]
    anova_result = stats.f_oneway(*grouped_vals)
    if post_hoc_test is True:
        tukey_result = pairwise_tukeyhsd(
            endog = list(input_data.loc[:, yvals_key]),
            groups = list(input_data.loc[:, xvals_key]),
            alpha = 0.05
            )
        return anova_result, tukey_result
    else:
        return anova_result

def barplot_data(
    input_data, xvals_key, yvals_key, color_index = None,
    ylimits = None, xlimits = None, x_label = None, y_label = None,
    output_stats = True, post_hoc_test = True, jitter_strength = 0.1, **kwargs
    ):
    xvalues = list(input_data.loc[:, xvals_key])
    grouped_data = input_data.groupby(xvals_key)[yvals_key]
    unique_xvalues = np.unique(xvalues)
    x_mapping = {val: i for i, val in enumerate(unique_xvalues)}
    # Generate random jitter for each point within each category
    jitter = {}
    for key in enumerate(unique_xvalues):
        jitter[key[1]] = np.random.uniform(-jitter_strength, jitter_strength, size = len(grouped_data.get_group(key[1])))
    if color_index is None:
        colors = {key: 'black' for key in xvalues}
    else:
        colors = {key: value for key in xvalues for value in list(input_data.loc[input_data[xvals_key] == key, color_index])}
    plt.figure(figsize = (5, 4))
    ax = plt.subplot(1, 1, 1)
    for x, y in grouped_data:
        # Scatter plot with jittered x-values
        ax.scatter(
            x = [x_mapping[x] + jitter[x][i] for i in range(len(y))],
            y = y.values,
            c = ([colors[x]] * len(y)),
            **kwargs
            )
        # Calculating and plotting SEM for each group
        sem = np.std(y) / np.sqrt(len(y))
        ax.bar(
            x,
            np.mean(y),
            color = 'none',
            edgecolor = colors[x]
            )
        ax.errorbar(
            x,
            np.mean(y),
            yerr = sem,
            fmt = '_',
            capsize = 5,
            color = 'black',
            label = 'SEM'
            )
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.xaxis.set_tick_params(width = 1.5)
    ax.yaxis.set_tick_params(width = 1.5)
    # Making ticks bold
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontweight('bold')
    if x_label is not None:
        if not isinstance(x_label, str):
            raise ValueError(f'{x_label} should be a string')
        plt.xlabel(
            xlabel = x_label,
            fontsize = 'large',
            fontweight = 'bold'
            )
    if y_label is not None:
        if not isinstance(y_label, str):
            raise ValueError(f'{y_label} should be a string')
        plt.ylabel(
            ylabel = y_label,
            fontsize = 'large',
            fontweight = 'bold'
            )
    if xlimits is not None:
        if not isinstance(xlimits, list):
            raise ValueError(f'{xlimits} should be a list dtype')
        if len(xlimits) != 2:
            raise ValueError(f'{xlimits} should be of length: 2')
        plt.xlim(xlimits)
    if ylimits is not None:
        if not isinstance(ylimits, list):
            raise ValueError(f'{ylimits} should be a list dtype')
        if len(ylimits) != 2:
            raise ValueError(f'{ylimits} should be of length: 2')
        plt.ylim(ylimits)
    plt.show()
    if output_stats is True:
        anov_output, tuckey_output = run_anova(
            input_data = input_data,
            xvals_key = xvals_key,
            yvals_key = yvals_key,
            post_hoc_test = post_hoc_test
            )
        print(anov_output)
        print(tuckey_output)


class IJ_data:
    
    def __init__(
        self,
        data: pd.DataFrame,
        normalization_data: list,
        blank_data: list = None,
        norm_blank_index: int = 0,
        control_indecies: Union[int, list[int]] = None,
        sample_labels: list = None,
        feature_labels: list = None
        ):
            self.data = data
            self.normalization_data = normalization_data
            self.blank_data = blank_data
            self.norm_blank_index = norm_blank_index
            self.control_indecies = control_indecies
            self.sample_labels = sample_labels
            self.feature_labels = feature_labels

    def zero(self, blank_key: Union[str, int] = 'Blank'):
        upper_range = len(np.transpose(self.data.values))
        for i in range(0, upper_range):
            self.data.iloc[:, i] = self.data.iloc[:, i].apply(
                func = subtract_blank,
                blank = self.blank_data[i]
                )
        self.normalization_data = list(pd.DataFrame({'protein': self.normalization_data}).apply(
            func = subtract_blank,
            axis = 1,
            blank = self.blank_data[self.norm_blank_index]
            )['protein'])
        # Remove the blank row from self.data
        self.data = self.data.loc[self.data.index != blank_key]

    def normalize(self):
        for i in range(0, len(self.data)):
            self.data.iloc[i, :] = self.data.iloc[i, :].apply(
                func = norm_division,
                normalization_factor = self.normalization_data[i]
                )

    def relative_expressions(self):
        upper_range = len(np.transpose(self.data.values))
        relative_values = self.data
        for i in range(0, upper_range):
            if isinstance(self.control_indecies, list):
                norm_value = np.mean(self.data.iloc[self.control_indecies, i])
            else:
                norm_value = self.data.iloc[self.control_indecies, i]
            relative_values.iloc[:, i] = self.data.iloc[:, i].apply(
                func = norm_division,
                normalization_factor = norm_value
                )
        return relative_values
    
    def main_method(self, blank_key = 'Blank') -> pd.DataFrame:
        if self.control_indecies is None:
            self.control_indecies = 0
        self.zero(blank_key = blank_key)
        self.normalize()
        ij_relative = self.relative_expressions()
        return ij_relative


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f',
        '--filepath',
        default = 'Western_blot/Rep0/ImageJ_data.csv',
        help = 'Path to raw BCA data'
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

    ij_object = IJ_data(
        data = input_data,
        normalization_data = list(input_data.loc[input_data.index != args.blank_key, args.normalization_key]),
        blank_data = list(input_data.loc[args.blank_key, :]),
        norm_blank_index = input_data.columns.get_loc(args.normalization_key),
        control_indecies = args.control_indecies
        )
    
    ij_object_relative = ij_object.main_method(blank_key = args.blank_key)
    ij_object_relative['Condition'] = [re.sub('-M[0-9]-[0-9]-(E|F)', '', re.sub(' \(1/2\)', '', x)) for x in list(ij_object_relative.index)]
    ij_object_relative['Color'] = (['orange'] * 4) + (['black'] * 4) + (['purple'] * 6)

    barplot_data(
        input_data = ij_object_relative,
        xvals_key = 'Condition',
        yvals_key = 'AR',
        color_index = 'Color',
        y_label = 'AR Intesity / Actin Intensity'
        )
