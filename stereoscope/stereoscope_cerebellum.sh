#!/bin/bash

conda init bash
source ~/.bashrc 
source activate stereoscope

/usr/bin/time -v stereoscope run --sc_cnt /storage/group/qul12/default/retrofit/Mouse/CerebellumPuck_sc_count.tsv --sc_labels /storage/group/qul12/default/retrofit/Mouse/CerebellumPuck_sc_label.tsv \
--st_cnt /storage/group/qul12/default/retrofit/Mouse/CerebellumPuck_counts.tsv -o mouse_res_record --gpu -lr 0.01