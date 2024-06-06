import argparse

class csv_mapper:

    def __init__(self, dictionary):
        self.dictionary = dictionary

    def write_to_csv(self, filename):
        columns = list(self.dictionary.keys())
        tsv_file = open(filename, 'w')
        for i in range(len(self.dictionary[columns[0]]) + 1):
            for t in range(len(columns)):
                if i == 0:
                    if t == (len(columns) - 1):
                        tsv_file.write(str(columns[t]))
                    else:
                        tsv_file.write(str(columns[t]) + ',')
                else:
                    if t == (len(columns) - 1):
                        tsv_file.write(self.dictionary[columns[t]][i - 1])
                    else:
                        tsv_file.write(self.dictionary[columns[t]][i - 1] + ',')
            tsv_file.write('\n')

class solution:
    
    def __init__(
        self,
        concentration,
        volume
        ):
            self.concentration = concentration
            self.volume = volume
            
def solution_list_maker(concentrations):
    lst_of_solutions = []
    for i in range(0, len(concentrations)):
        lst_of_solutions.append(solution(
            concentrations[i],
            0
            ))
    return lst_of_solutions
            
def volume_calculator(
    conc_1,
    conc_2,
    final_volume
    ):
        vol_to_add = (conc_2 * final_volume) / conc_1
        return vol_to_add
    
def volume_of_solvent_calculator(
    solution_1,
    solution_2
    ):
        volume_of_solvent = solution_2.volume * (1 - (solution_2.concentration / solution_1.concentration))
        return volume_of_solvent

def dilution_calculator(
    lst_of_solutions,
    desired_volume = 100
    ):
        lst_of_solutions = sorted(
            lst_of_solutions,
            key = lambda x: x.concentration
            )
        for i in range(0, len(lst_of_solutions)):
            if i == 0:
                lst_of_solutions[i].volume = desired_volume
            else:
                volume_of_concentrate = volume_calculator(
                    lst_of_solutions[i].concentration,
                    lst_of_solutions[i - 1].concentration,
                    lst_of_solutions[i - 1].volume
                    )
                lst_of_solutions[i].volume = desired_volume + volume_of_concentrate
        return lst_of_solutions
        
def solution_display(
    lst_of_sols,
    stock_concentration
    ):
        for i in range(0, len(lst_of_sols)):
            if i < len(lst_of_sols) - 1:
                print(
                    'Concentration: ' + str(lst_of_sols[i].concentration) + 
                    '    Total volume: ' + str(lst_of_sols[i].volume) + 
                    '    Volume of previous stock: ' + str(volume_calculator(
                        lst_of_sols[i + 1].concentration,
                        lst_of_sols[i].concentration,
                        lst_of_sols[i].volume
                        )) +
                    '    Volume of solvent: ' + str(volume_of_solvent_calculator(
                        lst_of_sols[i + 1],
                        lst_of_sols[i]
                        ))
                    )
            else:
                print('Concentration: ' + str(lst_of_sols[i].concentration) + 
                    '    Total volume: ' + str(lst_of_sols[i].volume) +
                    '    Volume of previous stock: ' + str(volume_calculator(
                        stock_concentration,
                        lst_of_sols[i].concentration,
                        lst_of_sols[i].volume
                        )) +
                    '    Volume of solvent: ' + str(volume_of_solvent_calculator(
                        solution(
                            stock_concentration,
                            0
                            ),
                        lst_of_sols[i]
                        ))
                    )
                
def output_table_to_csv(
    lst_of_solutions,
    stock_concentration,
    output_file = 'Table_of_solutions_output.csv'
    ):
        lst_of_solutions.reverse()
        output_dict = {}
        for i in range(len(lst_of_solutions)):
            if i == 0:
                output_dict['Concentration'] = [str(lst_of_solutions[i].concentration)]
                output_dict['Volume of previous stock'] = [str(volume_calculator(
                    conc_1 = stock_concentration,
                    conc_2 = lst_of_solutions[i].concentration,
                    final_volume = lst_of_solutions[i].volume
                    ))]
                output_dict['Volume of solvent'] = [str(volume_of_solvent_calculator(
                    solution(
                        concentration = stock_concentration,
                        volume = 0
                        ),
                    lst_of_solutions[i]
                    ))]
                output_dict['Total volume'] = [str(lst_of_solutions[i].volume)]
            else:
                output_dict['Concentration'].append(str(lst_of_solutions[i].concentration))
                output_dict['Volume of previous stock'].append(str(volume_calculator(
                    conc_1 = lst_of_solutions[i - 1].concentration,
                    conc_2 = lst_of_solutions[i].concentration,
                    final_volume = lst_of_solutions[i].volume
                    )))
                output_dict['Volume of solvent'].append(str(volume_of_solvent_calculator(
                    lst_of_solutions[i - 1],
                    lst_of_solutions[i]
                    )))
                output_dict['Total volume'].append(str(lst_of_solutions[i].volume))
        csv_class = csv_mapper(output_dict)
        csv_class.write_to_csv(output_file)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-l',
        '--limiting_vol',
        default = 100.0,
        type = float,
        help = 'Floating point value of the limiting volume which determines the final volume of solutions to be made'
        )
    
    parser.add_argument(
        '-d',
        '--dilutions',
        default = [0.5, 0.2],
        type = list,
        help = 'List of dilutions desired'
        )
    
    parser.add_argument(
        '-o',
        '--output_path',
        default = './Western_blot/BCA/Iteration_0/Calculate_dilutions/Serial_dilutions.csv'
        )

    args = parser.parse_args()

    dilutions = args.dilutions

    for dilution in dilutions:
        if not isinstance(dilution, float):
            raise TypeError(f'{dilution} must be of dtype float')

    lst_of_sols = solution_list_maker(dilutions)

    # Calculate necessary volumes for serial dilutions
    lst_of_sols = dilution_calculator(
        lst_of_sols,
        desired_volume = args.limiting_vol
        )

    # Writing to .csv file
    output_table_to_csv(
        lst_of_solutions = lst_of_sols,
        stock_concentration = 1.0,
        output_file = args.output_path
        )
