/*
This file produces Table 2. Time Series Validation.

Outline:
0. Set workspace
1. Set variable lists
2. Adjust simulated data
3. Table 2
*/
*******************************************************

*0. Set workspace
clear all
set more off
set mem 400m
use "panel.dta", clear

*******************************************************

*1. Set variable lists 
local controls_interact "t*numHH t*sav t*shg t*fracGM t*fractionL" 

*******************************************************

*2. Adjust simulated data
sort village t
local original_variable dynamicMF_simulated 
g dynamicMF_simulated_adjust = `original_variable' 

local adjust dynamicMF_simulated_adjust // This is the simulated takeup variable we adjust.  

*Determine value of empirical takeup at time t == 1 in each village
g empiricaltakeup_t1 = dynamicMF_empirical if t == 1 
bys village: egen temp = max(empiricaltakeup_t1)
replace empiricaltakeup_t1 = temp
drop temp

*Determine when simulated take up is greater than empirical takeup at time t == 1
gen greater = (`adjust' >= empiricaltakeup_t1) 
bys village: egen simulatedtakeup_t1 = min(t) if greater == 1 // This gives the first period we have simulated takeup greater than empirical takeup. 
bys village: egen temp = max(simulatedtakeup_t1)
replace simulatedtakeup_t1 = temp
drop temp 

*Adjust data 
g scale = simulatedtakeup_t1 - 1
bys village: replace `adjust' = `adjust'[_n+scale] if t != 0 

*******************************************************

*3. Table 2 

xi3: areg dynamicMF_empirical dynamicMF_simulated_adjust i.t if t>0, absorb(village) clust(village)
outreg2 dynamicMF_empirical dynamicMF_simulated_adjust using "ScienceRegressionTable2", excel paren(se) bdec(3) noaster replace

xi3: areg dynamicMF_empirical dynamicMF_simulated_adjust `controls_interact' i.t if t>0, absorb(village) clust(village)
outreg2 dynamicMF_empirical dynamicMF_simulated_adjust using "ScienceRegressionTable2", excel paren(se) bdec(3) noaster append 

*******************************************************
