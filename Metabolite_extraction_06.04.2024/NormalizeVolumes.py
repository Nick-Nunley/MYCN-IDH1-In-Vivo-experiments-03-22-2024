import pandas as pd
import argparse

def calculate_volumes(t_mass: float, init_vol: float = 1.0, desired_mass: float = 4.0, mass_in_g: bool = True):
    if mass_in_g is True:
        t_mass = t_mass * 1000
    scaling_percent = desired_mass / t_mass
    new_vol = init_vol * scaling_percent
    return new_vol

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f',
        '--file_path',
        default = './Harvest_05.28.2024/Metabolite_sample_weights.csv',
        help = 'Path sample weights'
        )
    
    parser.add_argument(
        '-w',
        '--weight_key',
        default = 'Mass (g)',
        help = 'Key for weight in input.csv'
        )
    
    parser.add_argument(
        '-o',
        '--output_path',
        default = './Metabolite_extraction_06.04.2024/Metabolite_sample_weights.csv'
        )

    args = parser.parse_args()

    input_data = pd.read_csv(args.file_path)

    input_data['Volume of methanol to extract (mL)'] = input_data.loc[:, args.weight_key].apply(
        func = calculate_volumes
        )
    
    input_data.to_csv(args.output_path, index = False)
