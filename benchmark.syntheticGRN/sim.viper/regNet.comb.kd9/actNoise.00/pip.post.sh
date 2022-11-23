#!/usr/bin/env bash

R CMD BATCH calAvgCor.R
R CMD BATCH plotAvgCor.R

# the following scripts are not needed any more
#R CMD BATCH cal.cor.across.TFs.R
#R CMD BATCH cal.dist.across.models.R
