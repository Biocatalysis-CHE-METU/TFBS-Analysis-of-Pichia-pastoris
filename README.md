# Transcription-Factor-Analysis-of-Pichia-pastoris
Original code prepared for the article "Functional Annotation of cis-acting DNA Elements Interconnecting Regulatory Networks of the Central Pathways in Pichia pastoris Through the Constructed Yeast Curation Pipeline"

Jupyter files Test_Factors and Test_Results are used to check for the accuracy for different transcription factor consensus matrices and allowed matches per 10000bp in secondary exons values, respectively.
Cut-Off checker calculates the cut-off values of each transcription factor for a given value of allowed matches per 10000bp in secondary exons, which was taken as 3 per 10000bp in this study.

Scanner takes the transcription factor position weight matrix database and list of promoters, to scan for putative cis-acting DNA elements by using the scan algorithm described in the article.
Aligner on the other hand, takes the results file  of the Scanner and list of orthologous promoters, to do pairwise alignment and seek for conserved putative cis-acting DNA elements. 
