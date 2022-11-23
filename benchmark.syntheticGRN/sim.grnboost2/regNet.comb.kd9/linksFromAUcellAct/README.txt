
Run the programs in the following format
create.input.R
    format input data into grnboost2 format 

run_grnboost2.sh
    predict interactions  
    output file: output.tsv

evalPerform-TF2targets.R
    input: predicted network - output.tsv 
    calculate recall and other metric to evaluate prediction performace 
    output: ./data/precision_recall.tsv 
            ./fig.TF2targets/stats.pdf

evalPerform-TF2TF.R
