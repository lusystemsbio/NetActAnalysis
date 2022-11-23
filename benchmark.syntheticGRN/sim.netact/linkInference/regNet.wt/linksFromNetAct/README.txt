
Links inferred and precision-recall calculated for regulon DBs
at different perturbation levelsi 0% (0), 25% (5), 50% (10), 
and 75% (15): regdb00/, regdb05/, regdb10/, regdb15/

Each folder has the following R scripts: 
functions.R
    functions defined in this file are used by other programs
create.input.R
    input: 
         (1) racipe simulated expressions 
         (2) NetAct calculated activities for TF - 
            for different perturbation level appropriate activity file 
            is used. 
    action: TF activities and target expressions are combined to create 
            activities for both TF and targets. 
    output: combined activity file - ./data/net30tf.states.tsv

inferLinks.R
    input: combined activity file - ./data/net30tf.states.tsv
    links are inferred from the activities and saved
    output: inferredLinks.tsv

calPRvalues.R
    input: inferred links - inferredLinks.tsv
    Precision-Recall values are calculated and saved
    output: ./data/prValues.tsv, ./figs/stats.pdf 

calPR4newlinks.R - only for regdb05, regdb10, and regdb15
    Calcualtes NetAct's performance in discovering new links. 
    For this, ground truth are the links in regdb00 that were replaced while 
    creating regdb05, rebdb10, and regdb15 and the counts for those links are 
    150, 300, and 450 respectively.  

    input: inferred links - inferredLinks.tsv
    Calculate Precision-Recall values for the links replaced while 
        creating perturbed regulon 
    output: ./data/prValues.tsv, ./figs/stats.pdf 

eval.inferredLinks.R
    Calculates precision, recall, and roc curve values (TPR and FPR) 
    Histogram binning needs to calibrated to make it suitable for both precision-recall 
    and roc curve values. Currently, it is tuned only for precision-recall. 
