
benchmark 
 benchmarking TF activity and regulatory link inference using a synthetic GRN 
 tfActivity/ - benchmark results for TF activity inference 
 linkInference/ - benchmark results for regulatory link inference 

The following folders have the data and simulation of different methods used in the benchmark.
regDBs/ - folder for regulons with and without perturbations  
regDBsAUCell/ - folder for regulons modifed according to AUCell protocol
regNet/ - topology of the synthetic regulatory network and its visualization 

sim.racipe/ - folder for racipe simulations of the regulatory network 
    sim.racipe.wt - racipe simulation without perturbation 
    sim.racipe.kd9 - racipe simulation with knockdown of TF9
    sim.racipe.comb.kd9 - combined racipe simulation:with 40 WT models and 43 knockdown models
sim.netact/
    activity and link inference using NetAct 
sim.nca/
    activity calculations using NCA
sim.viper/
    activity calculations using viper
sim.aucell.1/
    activity and link inference using AUCell with AUCell regulon protocol
sim.aucell.2/
    activity and link inference using AUCell with original regulons 
sim.genie3
    link inference using GENIE3
sim.grnboost2
    link inference using GRNBOOST2
sim.ppcor
    link inference using PPCOR
arboreto
    GENIE3 and GRNBOOST2 implementation 
