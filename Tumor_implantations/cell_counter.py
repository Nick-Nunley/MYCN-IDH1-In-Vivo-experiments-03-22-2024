import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '-i',
    '--input_path',
    default = 'Tumor_implantations/cell_count_input.csv',
    help = 'Path to input.csv file'
    )

parser.add_argument(
    '-o',
    '--output_path',
    default = 'Tumor_implantations/cell_count_output.csv',
    help = 'Path to output .csv'
    )

args = parser.parse_args()

class cell_counter:

    def __init__(
        self, raw_counts, groups = [], volumes = [],
        total_cells = [], cell_concentrations = []
        ):
            self.raw_counts = raw_counts
            self.groups = groups
            self.volumes = volumes
            self.total_cells = total_cells
            self.cell_concentrations = cell_concentrations

    def calculate_cell_concentration_wrapper(self, raw_count):
        cell_conc = raw_count * (10 ** 4) / 2.0
        return cell_conc

    def calculate_cell_concentrations(self):
        array_of_counts = np.array(self.raw_counts)
        concentrations = np.apply_along_axis(
            func1d = self.calculate_cell_concentration_wrapper,
            axis = 0,
            arr = array_of_counts
            )
        self.cell_concentrations = list(concentrations)

    def calculate_total_cells(self):
        for i in range(len(self.raw_counts)):
            cell_count = self.cell_concentrations[i] * self.volumes[i]
            self.total_cells.append(cell_count)

    def calculate(self):
        self.calculate_cell_concentrations()
        self.calculate_total_cells()

    def master_mix_wrapper(self, num_cells_to_plate, volume_to_plate_at, cell_concentration, scaling_factor):
        vol_of_cells = num_cells_to_plate * scaling_factor / cell_concentration
        vol_of_media = (volume_to_plate_at * scaling_factor) - vol_of_cells
        return [vol_of_cells, vol_of_media]

    def calculate_master_mixes(self, num_of_cells_to_plate, volumes_to_plate_at, scaling_factors = 2.0):
        if isinstance(scaling_factors, float) or isinstance(scaling_factors, int):
            scaling_factors = [scaling_factors] * len(self.raw_counts)
        if isinstance(num_of_cells_to_plate, float) or isinstance(num_of_cells_to_plate, int):
            num_of_cells_to_plate = [num_of_cells_to_plate] * len(self.raw_counts)
        if isinstance(volumes_to_plate_at, float) or isinstance(volumes_to_plate_at, int):
            volumes_to_plate_at = [volumes_to_plate_at] * len(self.raw_counts)
        for i in list(self.__dict__.keys()):
            if len(self.__dict__[i]) != len(self.raw_counts):
                new_value = ['NA'] * len(self.raw_counts)
                self.__dict__[i] = new_value 
        master_mix = {
            'Group': self.groups,
            'Raw_count': self.raw_counts,
            'Total_cells': self.total_cells,
            'Concentration': self.cell_concentrations,
            'Volume_of_cells_to_add': [],
            'Volume_of_media_to_add': [],
            'Total_volume': []
            }
        for i in range(len(self.raw_counts)):
            values_to_append = self.master_mix_wrapper(
                num_cells_to_plate = num_of_cells_to_plate[i],
                volume_to_plate_at = volumes_to_plate_at[i],
                cell_concentration = self.cell_concentrations[i],
                scaling_factor = scaling_factors[i]
                )
            master_mix['Volume_of_cells_to_add'].append(values_to_append[0])
            master_mix['Volume_of_media_to_add'].append(values_to_append[1])
            master_mix['Total_volume'].append(volumes_to_plate_at[i] * scaling_factors[i])
        return master_mix
    
    def write_to_csv(self, master_mix, filename):
        columns = list(master_mix.keys())
        csv_file = open(filename, 'w')
        for i in range(len(master_mix[columns[0]]) + 1):
            for t in range(len(columns)):
                if i == 0:
                    if t == (len(columns) - 1):
                        csv_file.write(str(columns[t]))
                    else:
                        csv_file.write(str(columns[t]) + ',')
                else:
                    if t == (len(columns) - 1):
                        csv_file.write(str(master_mix[columns[t]][i - 1]))
                    else:
                        csv_file.write(str(master_mix[columns[t]][i - 1]) + ',')
            csv_file.write('\n')


input_data = pd.read_csv(args.input_path)

cell_counter_class = cell_counter(
    raw_counts = list(input_data['count']),
    groups = list(input_data['Group']),
    volumes = list(input_data['volume'])
    )

cell_counter_class.calculate()

master_mix = cell_counter_class.calculate_master_mixes(
    num_of_cells_to_plate = list(input_data['cells_to_plate']),
    volumes_to_plate_at = list(input_data['volume_to_plate']),
    scaling_factors = list(input_data['scaling_factor'])
    )

cell_counter_class.write_to_csv(
    master_mix = master_mix,
    filename = args.output_path
    )
