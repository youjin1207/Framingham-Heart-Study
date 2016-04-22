/*
This file produces Table S3: Microfinance Take-Up versus Centralities of Leaders. 

Outline:
0. Set workspace
1. Set variable lists
2. Table S3
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

local controls "numHH sav shg fracGM fractionLeader"; 
#delimit cr 
*******************************************************

*2. Table S3
*Panel A: Basic Regressions
foreach var of varlist `leader_list1'{
	reg mf `var', robust
	outreg2 mf `var' using "TableS3_PanelA", excel paren(se) bdec(3) noaster append
}
forvalues i = 2(1)3{
	reg mf `leader_list`i'', robust
	outreg2 mf `leader_list`i'' using "TableS3_PanelA", excel paren(se) bdec(3) noaster append
}

*Panel B: Regressions with numHH
foreach var of varlist `leader_list1'{
	reg mf `var' numHH, robust
	outreg2 mf `var' numHH using "TableS3_PanelB", excel paren(se) bdec(3) noaster append
}
forvalues i = 2(1)3{
	reg mf `leader_list`i'' numHH, robust
	outreg2 mf `leader_list`i'' numHH using "TableS3_PanelB", excel paren(se) bdec(3) noaster append
}

*Panel C: Regressions with full controls
foreach var of varlist `leader_list1'{
	reg mf `var' `controls', robust
	outreg2 mf `var' `controls' using "TableS3_PanelC", excel paren(se) bdec(3) noaster append

}
forvalues i = 2(1)3{
	reg mf `leader_list`i'' `controls', robust 
	outreg2 mf `leader_list`i'' `controls' using "TableS3_PanelC", excel paren(se) noaster bdec(3) append
}

*******************************************************

