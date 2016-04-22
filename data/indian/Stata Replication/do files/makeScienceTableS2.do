/*
This file produces Table S2: Explaining the Average Centrality of Leaders. 

Outline:
0. Set workspace
1. Set variable lists
2. Table S2

*/
*******************************************************

*0. Set workspace
clear all
set more off
set mem 400m
use "cross_sectional.dta", clear

*******************************************************

*1. Set variable lists 

#delimit ; 
local leader_list1 "degree_leader eigenvector_centrality_leader between_centrality_leader bonacich_centrality_leader 
	decay_centrality_leader closeness_centrality_leader diffusion_centrality_leader communication_centrality_leader"; // full list
local controls "numHH sav shg fracGM fractionLeader"; 
#delimit cr 

*******************************************************

*2. Table S2
foreach measure of varlist `leader_list1'{
	foreach var of varlist `controls'{
		reg `measure' `var', robust
		outreg2 `measure' `var' using "TableS2_`measure'", excel paren(se) bdec(3) noaster append
	}
	reg `measure' `controls', robust  
	outreg2 `measure' `controls' using "TableS2_`measure'", excel paren(se) bdec(3) noaster append
}

*******************************************************

