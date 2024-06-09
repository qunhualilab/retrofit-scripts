#!/bin/bash

source activate stereoscope
#Mouse brain
stereoscope run --sc_cnt CerebellumPuck_sc_count.tsv --sc_labels CerebellumPuck_sc_label.tsv \
--st_cnt CerebellumPuck_counts.tsv -o mouse_res_record --gpu -lr 0.01
#Simulation
stereoscope run --sc_cnt Cerebellum_sc_count_K=10.tsv --sc_labels Cerebellum_sc_label_K=10.tsv \
--st_cnt N=10,M=3_loc_X.tsv -o N=10,M=3_res --gpu -lr 0.01
