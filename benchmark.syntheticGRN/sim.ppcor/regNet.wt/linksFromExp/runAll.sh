#!/usr/bin/env bash 

date 
R CMD BATCH create.input.R
R CMD BATCH inferLinks.R
R CMD BATCH calPRvalues.R
date
