#!/bin/bash

source activate stereoscope

/usr/bin/time -v stereoscope run --sc_cnt CerebellumPuck_sc_count.tsv --sc_labels CerebellumPuck_sc_label.tsv \
--st_cnt CerebellumPuck_counts.tsv -o mouse_res_record --gpu -lr 0.01
