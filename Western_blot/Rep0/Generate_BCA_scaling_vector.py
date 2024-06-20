import argparse
from typing import Union
import pandas as pd


def str_or_int(value):
    try:
        return int(value)
    except ValueError:
        return value

def main_element_transformation(data_value: float, bool_value: bool) -> float:
    if bool_value is True:
        return data_value
    else:
        return 1.0

def main_transformation(input_data: pd.DataFrame, main_col: Union[str, int], bool_col: Union[str, int]) -> pd.DataFrame:
    input_data = input_data.loc[:, [main_col, bool_col]]
    input_data = input_data.apply(
        func = lambda row: main_element_transformation(data_value = row[main_col], bool_value = row[bool_col]),
        axis = 1
        )
    return input_data


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i',
        '--input',
        default = 'Western_blot/Rep0/ImageJ_data_for_bca_scaling.csv',
        type = str,
        help = 'Path to raw ImageJ data'
        )
    
    parser.add_argument(
        '-c',
        '--column_key',
        default = 'Actin_relative_to_FUCRW-M1-1-E',
        type = str_or_int,
        help = 'Main column key for transformation'
        )
    
    parser.add_argument(
        '-b',
        '--bool_key',
        default = 'Use_for_adjustment',
        type = str_or_int,
        help = 'Key indicating column for boolean vector'
        )

    parser.add_argument(
        '-o',
        '--output',
        default = 'Western_blot/Rep0/BCA_scaling_vector.csv',
        type = str,
        help = 'Output path'
        )

    args = parser.parse_args()

    input_data = pd.read_csv(args.input)

    input_data = main_transformation(
        input_data = input_data,
        main_col = args.column_key,
        bool_col = args.bool_key
        )

    input_data.to_csv(
        path_or_buf = args.output,
        index = False,
        header = False
        )
