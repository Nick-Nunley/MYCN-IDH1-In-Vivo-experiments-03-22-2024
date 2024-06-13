#!/bin/bash

python Submodules/BCA-Assay/BCA_assay.py \
-s Western_blot/BCA/standards_data.xlsx \
-u Western_blot/BCA/unknown_data_subset.xlsx \
-p 10.0 \
-a \
-f 4.5 \
-o Western_blot/BCA/bca_output_normalized_to_actin.csv \
-S \
-n Western_blot/Rep0/BCA_scaling_vector.csv
