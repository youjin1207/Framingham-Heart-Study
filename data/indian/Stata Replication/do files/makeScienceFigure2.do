/*
This file produces Figure 2: Microfinance participation versus measures of leader centrality. 

Outline:
0. Set workspace
1. Figure 2
*/
*******************************************************

*0. Set workspace
clear all
set more off
set mem 400m
use "cross_sectional.dta", clear
set scheme s1mono

*******************************************************

*1. Figure 2

*Panel A: Degree of leader
twoway (lfitci mf degree_leader, clpattern(solid) lcolor(blue) ciplot(rscatter) mcolor(brown) msymbol(p) msize(small)) (scatter mf degree_leader, msymbol(smplus) mcolor(sienna)), ytitle(Microfinance take-up rate) xtitle(Degree of leaders) legend(off) graphregion(color(white)) subtitle(A., pos(10))
graph set eps fontface Helvetica;
graph export graph_Deg_MF.eps, replace

*Panel B: Communication centrality of leader 
twoway (lfitci mf communication_centrality_leader, clpattern(solid) lcolor(blue) ciplot(rscatter) mcolor(brown) msymbol(p) msize(small)) (scatter mf communication_centrality_leader, msymbol(smplus) mcolor(sienna)), ytitle(Microfinance take-up rate) xtitle(Communication centrality of leaders) legend(off) graphregion(color(white)) subtitle(B., pos(10))
graph set eps fontface Helvetica;
graph export graph_Comm_MF.eps, replace

*Panel C: Diffusion centrality of leader 
twoway (lfitci mf diffusion_centrality_leader, clpattern(solid) lcolor(blue) ciplot(rscatter) mcolor(brown) msymbol(p) msize(small)) (scatter mf diffusion_centrality_leader, msymbol(smplus) mcolor(sienna)), ytitle(Microfinance take-up rate) xtitle(Diffusion centrality of leaders) legend(off) graphregion(color(white)) subtitle(C., pos(10))
graph set eps fontface Helvetica;
graph export graph_Diff_MF.eps, replace

*******************************************************




