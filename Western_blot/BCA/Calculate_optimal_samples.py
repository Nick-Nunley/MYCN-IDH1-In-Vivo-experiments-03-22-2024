import argparse
import re
import numpy as np
import pandas as pd

class BCAdata:

    def __init__(self, numeric_dataframe, sample_id_dataframe):
        self.numeric_dataframe = numeric_dataframe
        self.sample_id_dataframe = sample_id_dataframe

    def obtain_max_value(self, max_keys = None):
        if max_keys is None:
            max_keys = [['A', '9'], ['B', '9']]
        row_entries = [max_keys[0][0], max_keys[1][0]]
        col_entries = [max_keys[0][1], max_keys[1][1]]
        max_value = np.mean(self.numeric_dataframe.loc[row_entries, col_entries].to_numpy())
        return max_value

    def return_optimal_sample_id(self, sample_regex, max_val):
        mask = self.sample_id_dataframe.applymap(lambda x: bool(sample_regex.match(x)))
        matching_indices = mask[mask].stack().index.tolist()
        if len(matching_indices) == 0:
            return
        subset_data = []
        for index in matching_indices:
            row_index, col_index = index
            subset_data.append(self.numeric_dataframe.loc[row_index, col_index])
        subset_numeric_data = pd.DataFrame(subset_data, index = matching_indices, columns = ['Value'])
        optimal_row, optimal_col = subset_numeric_data.loc[subset_numeric_data['Value'] <= max_val, 'Value'].idxmax()
        optimal_sample_id = self.sample_id_dataframe.loc[optimal_row, optimal_col]
        return optimal_sample_id
    
    def main_calculations(self, sample_regexs, max_keys = None):
        max_absorbance = self.obtain_max_value(max_keys = max_keys)
        optimal_samples = [self.return_optimal_sample_id(sample_regex = reg_val, max_val = max_absorbance) for reg_val in sample_regexs]
        return optimal_samples
    
    def generate_unknowns_data(self, optimum_samples = None, sample_regexs = None, max_keys = None):
        if optimum_samples is None:
            if sample_regexs is None:
                raise ValueError(f'{sample_regexs} must be specified to do automatic sample optimization')
            optimum_samples = self.main_calculations(sample_regexs = sample_regexs, max_keys = max_keys)
            optimum_samples = [x for x in optimum_samples if x is not None]
        mask = self.sample_id_dataframe.isin(optimum_samples)
        index_values = self.sample_id_dataframe[mask].stack().index.tolist()
        subset_data = {}
        for val in index_values:
            row_index, col_index = val
            subset_data[self.sample_id_dataframe.loc[row_index, col_index]] = [self.numeric_dataframe.loc[row_index, col_index]]
        return subset_data
    
def main(absorbance_filepath, sample_filepath, sample_regexs, max_keys = None, output_path = 'unknowns_data.xlsx'):
    raw_data = pd.read_csv(absorbance_filepath)
    raw_data = raw_data.set_index(keys = '0')
    sample_keys = pd.read_csv(sample_filepath)
    sample_keys = sample_keys.set_index(keys = '0')
    sample_keys = sample_keys.fillna('') 
    bca_test = BCAdata(raw_data, sample_keys)
    optimal_samples = bca_test.main_calculations(sample_regexs = sample_regexs, max_keys = max_keys)
    optimal_samples = [x for x in optimal_samples if x is not None]
    print(f'Optimal samples are as follows: {optimal_samples}')
    unknowns_data = bca_test.generate_unknowns_data(optimum_samples = optimal_samples, max_keys = max_keys) 
    ordered_unknowns_data = {'Replicate': [int(1)]}
    for key in optimal_samples:
        if key in unknowns_data:
            ordered_unknowns_data[key] = unknowns_data[key]
    ordered_unknowns_data = pd.DataFrame(ordered_unknowns_data)
    ordered_unknowns_data.to_excel(output_path, index = False)
        

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f',
        '--filepath',
        default = 'Western_blot/BCA/BCA_data.csv',
        help = 'Path to raw BCA data'
        )
    
    parser.add_argument(
        '-s',
        '--sample_key',
        default = 'Western_blot/BCA/Plating_scheme.csv',
        help = 'Path to plating scheme for samples'
        )
    
    parser.add_argument(
        '-m',
        '--max_keys',
        default = [['A', '9'], ['B', '9']],
        type = list,
        help = '2D list containing max standard key entries where inner list is of length 2 containing row and column keys'
        )
    
    parser.add_argument(
        '-o',
        '--output_path',
        default = 'Western_blot/BCA/unknown_data.xlsx',
        help = 'Output path for uknowns_data.xlsx'
        )

    args = parser.parse_args()

    SAMPLE_REGEXS = [
        r'FUCRW-M1-1-E',
        r'FUCRW-M1-1-F',
        r'FUCRW-M1-2-E',
        r'FUCRW-M1-2-F',
        r'FUCRW-M2-1-E',
        r'FUCRW-M2-1-F',
        r'FUCRW-M2-2-E',
        r'FUCRW-M2-2-F',
        r'MYCN-M1-1-E',
        r'MYCN-M1-1-F',
        r'MYCN-M1-2-E',
        r'MYCN-M1-2-F',
        r'MYCN-M2-1-E',
        r'MYCN-M2-1-F',
        r'MYCN-M2-2-E',
        r'MYCN-M2-2-F',
        r'IDH1-M1-1-E',
        r'IDH1-M1-1-F',
        r'IDH1-M1-2-E',
        r'IDH1-M1-2-F',
        r'IDH1-M2-1-E',
        r'IDH1-M2-1-F',
        r'IDH1-M2-2-E',
        r'IDH1-M2-2-F',
        r'MYCN;IDH1-M1-1-E',
        r'MYCN;IDH1-M1-1-F',
        r'MYCN;IDH1-M1-2-E',
        r'MYCN;IDH1-M1-2-F',
        r'MYCN;IDH1-M2-1-E',
        r'MYCN;IDH1-M2-1-F',
        r'MYCN;IDH1-M2-2-E',
        r'MYCN;IDH1-M2-2-F'
        ]
    
    sample_regexs = [re.compile(val) for val in SAMPLE_REGEXS]

    main(
        absorbance_filepath = args.filepath,
        sample_filepath = args.sample_key,
        sample_regexs = sample_regexs,
        max_keys = args.max_keys,
        output_path = args.output_path
        )
