import argparse
import csv
import re
import numpy as np
            
class NormalizeBca:

    def __init__(self, total_volume: float = None, scaling_factor: float = None):
        self.total_volume = total_volume
        self.scaling_factor = scaling_factor

    def load_init_data(self, filepath: str, sep: str = ',', lines_to_skip: int = 3):
        """Function to convert a csv file to a dictionary"""
        data_dict = {}
        with open(filepath, encoding = 'UTF-8') as file:
            # Skip the first 3 lines
            for _ in range(lines_to_skip):
                next(file)
            keys = file.readline().strip('\n').replace('"', '').split(sep)
            for key in keys:
                data_dict[key] = []
            csv_file = csv.reader(
                file,
                delimiter = sep
                )
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
        """Function to load in the scaling vector"""
        data_vec = []
        with open(input_path, 'r', encoding = 'UTF-8') as file:
            reader = csv.reader(file)
            for line in reader:
                data_vec.append(float(line[0]))
        return data_vec
    
    def obtain_scaling_factor(self, filepath: str):
        with open(filepath, 'r', encoding = 'UTF-8') as file:
            file_header = file.readline().strip()
            float_match = re.search(r'(\d+\.\d+)', file_header)
            if float_match:
                return float(float_match.group(1))
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

    def write_to_csv(self, dictionary: dict, filename: str, columns: list = None):
        """Function to write a dictionary object to a .csv file"""
        if columns is None:
            columns = list(dictionary.keys())
        with open(filename, 'w', encoding = 'UTF-8') as csv_file:
            for index in range(len(dictionary[columns[0]]) + 1):
                for col_t in enumerate(columns):
                    if index == 0:
                        if col_t[0] == (len(columns) - 1):
                            csv_file.write(str(columns[col_t[0]]))
                        else:
                            csv_file.write(str(columns[col_t[0]]) + ',')
                    else:
                        if col_t[0] == (len(columns) - 1):
                            csv_file.write(str(dictionary[columns[col_t[0]]][index - 1]))
                        else:
                            csv_file.write(str(dictionary[columns[col_t[0]]][index - 1]) + ',')
                csv_file.write('\n')

    def main_function(self, input_filepath: str, norm_path: str, output_path: str):
        """Main function to wrap everything into"""
        data_dict = self.load_init_data(filepath = input_filepath)
        scaling_vector = self.load_scaling_vector(input_path = norm_path)
        data_dict = self.scale_protein_volumes(
            input_dict = data_dict,
            scaling_vector = scaling_vector
            )
        self.write_to_csv(
            dictionary = data_dict,
            filename = output_path
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

    NormBCA.main_function(
        args.initial_input,
        args.normalize_vector,
        args.output_path
        )
