#!/bin/bash

python Submodules/BCA-Assay/BCA_assay.py \
-s Western_blot/BCA/standards_data.xlsx \
-u Western_blot/BCA/unknown_data.xlsx \
-p 10.0 \
-a \
-f 1.25 \
-o Western_blot/BCA/bca_output_loading_control_only.csv \
