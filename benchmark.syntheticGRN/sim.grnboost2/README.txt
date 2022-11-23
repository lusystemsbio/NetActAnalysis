create.input.R
    format input data into grnboost2 format 
    output: ./net30tf.states.tsv, transcription_factors.tsv 

run_grnboost2.sh
    predict interactions  

evalPerform-TF2targets.R
    calculate recall and other metric to evaluate prediction performace 
    output: 
        ./data/precision-recall.tsv 
        ./figs.TF2targets/stats.pdf 



evalPerform-TF2TF.R
