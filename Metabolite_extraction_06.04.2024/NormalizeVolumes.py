import pandas as pd
import argparse

def calculate_volumes(t_mass: float, init_vol: float = 1.0, desired_mass: float = 4.0, mass_in_g: bool = True):
    if mass_in_g is True:
        t_mass = t_mass * 1000
    scaling_percent = desired_mass / t_mass
    new_vol = init_vol * scaling_percent
    return new_vol

def calculate_solv_volume(volume_val: float, init_vol: float = 1.0):
    return init_vol - volume_val

def main_volume_calculation(
    input_file: str,
    output_path: str,
    weight_key: str = 'Mass (g)',
    init_vol: float = 1.0,
    desired_mass: float = 4.0,
    mass_in_g: bool = True
    ):
        input_data = pd.read_csv(input_file)
        input_data['Volume of methanol to extract (mL)'] = input_data.loc[:, weight_key].apply(
            func = calculate_volumes,
            init_vol = init_vol,
            desired_mass = desired_mass,
            mass_in_g = mass_in_g
            )
        input_data['Volume of methanol to add (mL)'] = input_data.loc[:, 'Volume of methanol to extract (mL)'].apply(
            func = calculate_solv_volume,
            init_vol = init_vol
            )
        input_data.insert(0, 'Vial-ID', range(1, len(input_data) + 1))
        input_data.to_csv(output_path, index = False)
        return

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

    main_volume_calculation(
        input_file = args.file_path,
        output_path = args.output_path,
        weight_key = args.weight_key
        )
