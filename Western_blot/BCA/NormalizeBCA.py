import argparse
import csv
import re
import numpy as np

class NormalizeBca:

    def __init__(self, total_volume: float = None, scaling_factor: float = None):
        self.total_volume = total_volume
        self.scaling_factor = scaling_factor

    def load_init_data(self, filepath: str, sep: str = ',', lines_to_skip: int = 3):
        """Method to convert a csv file to a dictionary"""
        data_dict = {}
        with open(filepath, encoding = 'UTF-8') as file:
            # Skip the first 3 lines
            for _ in range(lines_to_skip):
                next(file)
            keys = file.readline().strip('\n').replace('"', '').split(sep)
            for key in keys:
                data_dict[key] = []
            csv_file = csv.reader(file, delimiter = sep)
            for line in csv_file:
                if len(line) == 0:
                    break
                else:
                    for key in enumerate(list(data_dict.keys())):
                        if key[1] == 'Sample':
                            data_dict[key[1]].append(str(line[key[0]]))
                        else:
                            data_dict[key[1]].append(float(str(line[key[0]])))
            return data_dict

    def load_scaling_vector(self, input_path: str):
        """Method to load in the scaling vector"""
        data_vec = []
        with open(input_path, 'r', encoding = 'UTF-8') as file:
            reader = csv.reader(file)
            for line in reader:
                data_vec.append(float(line[0]))
        return data_vec

    def obtain_scaling_factor(self, filepath: str):
        """Method to obtain scaling factor from initial input file"""
        with open(filepath, 'r', encoding = 'UTF-8') as file:
            file_header = file.readline().strip()
            float_match = re.search(r'(\d+\.\d+)', file_header)
            if float_match:
                self.scaling_factor = float(float_match.group(1))
            else:
                raise ValueError('Could not find a scaling factor in input file')

    def scale_protein_volumes(self, input_dict: dict, scaling_vector: list):
        """Scale and shift 'Protein' and 'DI Water' fields of input_dict accordingly"""
        scaling_vector = np.array(scaling_vector)
        old_proteins = np.array(input_dict['Protein'])
        new_proteins = old_proteins / scaling_vector
        volume_deltas = new_proteins - old_proteins
        new_waters = np.array(input_dict['DI Water']) - volume_deltas
        # Shifting everything to make min water value new 0 and recalculating 'Reducing agent' and 'Sample Buffer' volumes
        min_water = np.min(new_waters)
        new_proteins = new_proteins - min_water
        new_waters = new_waters - min_water
        total_volume = (20.0 / 13.0) * max(new_proteins)
        input_dict['Reducing Agent'] = [total_volume / 10.0] * len(input_dict['Sample'])
        input_dict['Sample Buffer'] = [total_volume / 4.0] * len(input_dict['Sample'])
        # Overwriting old data
        input_dict['Protein'] = list(new_proteins)
        input_dict['DI Water'] = list(new_waters)
        return input_dict

    def calculate_total_volume(self, input_dict: dict):
        """Method to calculate total volume"""
        if self.scaling_factor is None:
            raise ValueError(f'{self.scaling_factor} must be a floating point value. Try calling self.obtain_scaling_factor()')
        # Assume that all total volumes are the same for efficiency and calculate first entries
        total_vol = []
        for key in enumerate(input_dict):
            if key[1] == 'Sample':
                continue
            else:
                total_vol.append(input_dict[key[1]][0])
        self.total_volume = np.sum(np.array(total_vol)) / self.scaling_factor

    def write_to_csv(
        self,
        dictionary: dict,
        filename: str,
        lines: list,
        data_start: int = 4,
        data_end: int = 26,
        columns: list = None
        ):
        """Method to write a dictionary object to a .csv file along with header and footer"""
        if columns is None:
            columns = list(dictionary.keys())
        with open(filename, 'w', encoding = 'UTF-8') as csv_file:
            # Write the header and any lines before the data section
            for line in lines[:data_start]:
                csv_file.write(line)
            # Check if columns headers already exist in the lines before data_start
            if lines[data_start-1].strip() != ','.join(columns):
                csv_file.write(','.join(columns) + '\n')
            # Write the data
            for index in range(len(dictionary[columns[0]])):
                row = [str(dictionary[col][index]) for col in columns]
                csv_file.write(','.join(row) + '\n')
            # Write any lines after the data section
            for line in lines[data_end:]:
                # Replace the Total Volume line with the updated total volume
                if line.startswith('Total Volume:'):
                    line = f'Total Volume:,{self.total_volume}\n'
                csv_file.write(line)

    def main_method(self, input_filepath: str, norm_path: str, output_path: str):
        """Main method to wrap everything into"""
        # Read the entire file and split into sections
        with open(input_filepath, 'r', encoding='UTF-8') as file:
            lines = file.readlines()
        # Extract header, data, and footer sections
        data_start = 4  # Assuming the data starts from the 4th line
        data_end = len(lines)
        for i in range(4, len(lines)):
            if lines[i].strip() == '':
                data_end = i
                break
        # Load initial data into a dictionary
        data_dict = self.load_init_data(filepath = input_filepath)
        # Load the scaling vector
        scaling_vector = self.load_scaling_vector(input_path = norm_path)
        self.obtain_scaling_factor(filepath = input_filepath)
        # Scale protein volumes
        data_dict = self.scale_protein_volumes(
            input_dict = data_dict,
            scaling_vector = scaling_vector
            )
        # Calculate total volume
        self.calculate_total_volume(input_dict = data_dict)
        # Write the updated content back to a CSV file
        self.write_to_csv(
            dictionary = data_dict,
            filename = output_path,
            lines = lines,
            data_start = data_start,
            data_end = data_end
            )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i',
        '--initial_input',
        default = 'Western_blot/BCA/bca_output.csv',
        help = 'Path to initial bca_output.csv'
        )

    parser.add_argument(
        '-n',
        '--normalize_vector',
        default = 'Western_blot/BCA/test_vector.csv',
        help = 'Path to list of scaling values'
        )

    parser.add_argument(
        '-o',
        '--output_path',
        default = 'Western_blot/BCA/test_bca_output.csv',
        help = 'Output path for scaled bca_output.csv'
        )

    args = parser.parse_args()

    NormBCA = NormalizeBca()

    NormBCA.main_method(
        input_filepath = args.initial_input,
        norm_path = args.normalize_vector,
        output_path = args.output_path
        )
