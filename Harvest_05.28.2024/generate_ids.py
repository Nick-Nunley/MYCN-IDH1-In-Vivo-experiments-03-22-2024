import json

MOUSE_PHENOS = ['FUCRW', 'MYCN', 'IDH1', 'MYCN;IDH1']
NUMBER_OF_MICE = 2
TUMOR_NUM = 2
TECH_REPS = 6


def generate_sample_ids(mouse_phenos, number_of_mice, tumor_num, tech_reps):
    sample_json = {}
    for pheno in mouse_phenos:
        for mouse_num in range(number_of_mice):
            for tumor in range(tumor_num):
                for rep in range(tech_reps):
                    sample_json_key = f'{pheno}-M{mouse_num + 1}-{tumor + 1}'
                    if rep == 0:
                        sample_json_key = f'{sample_json_key}-A'
                    elif rep == 1:
                        sample_json_key = f'{sample_json_key}-B'
                    elif rep == 2:
                        sample_json_key = f'{sample_json_key}-C'
                    elif rep == 3:
                        sample_json_key = f'{sample_json_key}-D'
                    elif rep == 4:
                        sample_json_key = f'{sample_json_key}-E'
                    elif rep == 5:
                        sample_json_key = f'{sample_json_key}-F'
                    else:
                        raise IndexError(f'More technical reps were specified than can be handled')
                    if '-A' in sample_json_key or '-B' in sample_json_key or '-C' in sample_json_key:
                        sample_json[sample_json_key] = 'Metabolomics'
                    elif '-D' in sample_json_key:
                        sample_json[sample_json_key] = 'Tissue fixation'
                    elif '-E' in sample_json_key or '-F' in sample_json_key:
                        sample_json[sample_json_key] = 'Protein lysate'
                    else:
                        raise ValueError(f'{sample_json_key} is missing a technical ID')

    with open('Harvest_05.28.2024/sample_ids.json', 'w', encoding = 'UTF-8') as json_file:
        json.dump(sample_json, json_file, indent = 4)

if __name__ == '__main__':
    generate_sample_ids(
        mouse_phenos = MOUSE_PHENOS,
        number_of_mice = NUMBER_OF_MICE,
        tumor_num = TUMOR_NUM,
        tech_reps = TECH_REPS
        )

