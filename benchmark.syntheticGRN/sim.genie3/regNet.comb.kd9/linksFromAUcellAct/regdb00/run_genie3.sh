#!/usr/bin/env bash

python ../../../../arboreto/scripts/run_arboretum.py  --genie3 \
-i ./data/net30tf.states.tsv \
-tf ./data/transcription_factors.tsv \
-o output.tsv \
> run_log.txt
