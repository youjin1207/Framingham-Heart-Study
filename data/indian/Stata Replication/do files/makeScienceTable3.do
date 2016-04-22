/*
This file produces Table 3: Microfinance Take-Up versus Centralities of Leaders. 

Outline:
0. Set workspace
1. Set variable lists
2. Table 3
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

local leader_list1 "communication_centrality_leader diffusion_centrality_leader degree_leader eigenvector_centrality_leader
	between_centrality_leader bonacich_centrality_leader decay_centrality_leader closeness_centrality_leader"; // full list

local leader_list2 "communication_centrality_leader degree_leader eigenvector_centrality_leader 
	between_centrality_leader bonacich_centrality_leader decay_centrality_leader closeness_centrality_leader"; // remove diffusion_centrality_leader

local leader_list3 "diffusion_centrality_leader degree_leader eigenvector_centrality_leader
	between_centrality_leader bonacich_centrality_leader decay_centrality_leader closeness_centrality_leader"; // remove communication_centrality_leader

local short_list "communication_centrality_leader diffusion_centrality_leader";

local controls "numHH sav shg fracGM fractionLeader"; 

#delimit cr

*******************************************************

*2. Table 3 
foreach var of varlist `short_list'{
	reg mf `var' `controls', robust 
	outreg2 mf `var' using "ScienceRegressionTable3", excel paren(se) bdec(3) noaster append
}

forvalues i = 2(1)3{
	reg mf `leader_list`i'' `controls', robust
	outreg2 mf `leader_list`i'' using "ScienceRegressionTable3", excel paren(se) bdec(3) noaster append
}

*******************************************************
