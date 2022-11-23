#!/usr/bin/env bash

date
python ../../../arboreto/scripts/run_arboretum.py  --grnboost2 \
-i ./data/net30tf.states.tsv \
-tf ./data/transcription_factors.tsv \
-o output.tsv \
> run_log.txt
date
