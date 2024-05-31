import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

class TumorMass:

    def __init__(self, tumor_data = None):
        self.tumor_data = tumor_data

    def read_in_data(self, data_path, **kwargs):
        self.tumor_data = pd.read_csv(data_path, **kwargs)

    def run_anova(self, xvals_key, yvals_key, post_hoc_test = True):
        unique_x_vals = list(set(list(self.tumor_data.loc[:, xvals_key])))
        grouped_vals = [list(self.tumor_data.loc[self.tumor_data[xvals_key] == unique_value, yvals_key]) for unique_value in unique_x_vals]
        anova_result = stats.f_oneway(*grouped_vals)
        if post_hoc_test is True:
            tukey_result = pairwise_tukeyhsd(
                endog = list(self.tumor_data.loc[:, yvals_key]),
                groups = list(self.tumor_data.loc[:, xvals_key]),
                alpha = 0.05
                )
            return anova_result, tukey_result
        else:
            return anova_result

    def plot_data(
            self, xvals_key, yvals_key,
            ylimits = None, xlimits = None, x_label = None, y_label = None,
            output_stats = True, post_hoc_test = True, **kwargs
            ):
            xvalues = list(self.tumor_data.loc[:, xvals_key])
            yvalues = list(self.tumor_data.loc[:, yvals_key])
            grouped_data = self.tumor_data.groupby(xvals_key)[yvals_key]
            plt.figure(figsize = (5, 4))
            ax = plt.subplot(1, 1, 1)
            ax.scatter(
                x = xvalues,
                y = yvalues,
                **kwargs
                )
            # Calculating and plotting SEM for each group
            for x, y in grouped_data:
                sem = np.std(y) / np.sqrt(len(y))
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
                anov_output, tuckey_output = self.run_anova(
                    xvals_key = xvals_key,
                    yvals_key = yvals_key,
                    post_hoc_test = post_hoc_test
                    )
                print(anov_output)
                print(tuckey_output)
        
