import json

mouse_phenos = ['FUCRW', 'MYCN', 'IDH1', 'MYCN;IDH1']
number_of_mice = 2
tumor_num = 2
tech_reps = 6

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
                sample_json[sample_json_key] = ''

with open('Harvest_05.28.2024/sample_ids.json', 'w', encoding = 'UTF-8') as json_file:
    json.dump(sample_json, json_file, indent = 4)

