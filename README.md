# Network Dependence Can Lead to Spurious Associations and Invalid Inference

Youjin Lee and Elizabeth L. Ogburn

## Abstract

Researchers across the health and social sciences generally assume that observations are independent, even while relying on convenience samples that draw subjects from one or a small number of communities, schools, hospitals, etc. A paradigmatic example of this is the Framingham Heart Study (FHS). Many of the limitations of such samples are well-known, but the issue of statistical dependence due to social network ties has not previously been addressed. We show that, along with anticonservative variance estimation, this can result in spurious associations due to network dependence. Using a statistical test that we adapted from one developed for spatial autocorrelation, we test for network dependence in several of the thousands of influential papers that have been published using FHS data. Results suggest that some of the many decades of research on coronary heart
disease, other health outcomes, and peer influence using FHS data may suffer from spurious associations, error-prone point estimates, and anticonservative inference due to unacknowledged network dependence. These issues are not unique to the FHS; as researchers in psychology, medicine, and beyond grapple with replication failures, this unacknowledged source of invalid statistical inference should be part of the conversation.


## Instructions for Use

All table and figures from our [paper](https://arxiv.org/pdf/1908.00520.pdf) can be reproducible by implementing following `R` code.


Code                      Data                      Figure        Table
----------------------    ----------------          ----------    -----------------------------
code/ReadPeerConti.R      data/PeerConti.RData      Figure 2      Table 1         
code/ReadLatentConti.R    data/LatentConti.RData    Figure A1                
code/sim2.R                                         Figure 3      Table 2
code/sim3.R                                         Figure 1
code/FHSsim.R                                       Figure 4      Table 3, Table A1
code/FHSanalysis.R                                  Figure A3     Table 4, Table A2, Table A3       
code/HRVexample.R                                                 Table 5, Table A4, Table A5
code/FHSstrokestudy.R                                             Table A6
code/FHSCHDstudy.R                                                Table A7      
code/FHSCVDstudy.R                                                Table A8
code/Mendelian.R                                                  Section A4.3                  
-------------------------------------------------------------------------------------------------

We also created the [`netdep`](https://github.com/youjin1207/netdep) `R` package available at [CRAN](https://cran.r-project.org/web/packages/netdep/index.html). 


## Data

Our work is primarily based on Framingham Heart Study (FHS) so as to examine whether the health-related outcomes are dependent on their underlying social network. As health-related data and social network data are separately saved for each subject, the analysis requires two data sources, [phs000153.v9.p8 (Framingham SHARe Social Network)](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000153.v9.p8) and  [phs000007.v29.p10 (Framingham Cohort)](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000007.v29.p10). 

Researchers may request access to the data through the [dbGaP website](https://www.ncbi.nlm.nih.gov/gap/).
dbGap also provides [tips for preparing a data access request](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetPdf.cgi?document_name=GeneralAAInstructions.pdf), where principal investigators (PIs) are required to receive an Institutional Review Board (IRB) approval from their local institutions and to submit a Research Use Statement.

For the analyses described in the main manuscript, there are two structures for the data – (1) social network data and (2) phenotype data. Variables associated with social network data are listed [here](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/dataset.cgi?study_id=phs000153.v9.p8&phv=73855&phd=1560&pha=&pht=836&phvf=1&phdf=&phaf=&phtf=1&dssp=1&consent=&temp=1). For our analyses, we constructed the adjacency matrix for each network regardless of `ALTERTYPE`, i.e., $A_{ij}$=1 if and only if $i$ and $j$ have any kind of tie, unless otherwise mentioned (e.g., sibling relationships in Section 7.2). We used `CAUSEINIT` and `CAUSESEVERED` to determine the existence of edges in replicating the studies in five published papers; for instance, if we are only interested in the underlying social network of the subjects at the time of the Original Cohort Exam 16 (1979 - 1982), we only considered (approximately) the edges that were initiated no later than 1982 and were not severed by 1978.  
A list of variables can be found at [dbGaP website for each exam](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?study_id=phs000007.v30.p11&phv=159482&phd=4398&pha=4313&pht=3099&phvf=1&phdf=&phaf=&phtf=&dssp=1&consent=&temp=1). Indices from the adjacency matrix described above were matched with phenotype data through a variable named “shareid” (SHARe ID).  We used variables from the publicly available data set without changing their names, so which phenotype variables were used for each analysis can be found in the code as described in the above table.


