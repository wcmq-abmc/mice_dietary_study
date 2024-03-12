Data sets for this analysis came from a mice dietary study where 17 mice models were fed either a controlled diet (5 mice), 50% reduced protein diet (6 mice), or 50% reduced non-essential amino acid (6 mice) for two weeks. For each mouse in the three groups, three types of tissues were collected: i) four pellets of feces to examine changes 
in the metabolome and microbiome, ii) blood samples via cardiac puncture to examine circulating metabolites level, and iii) liver samples. Fecal, serum, and liver samples were sent out for metabolomics profiling. Fecal samples were also sent out for sequencing. All data was processed and analyzed using the “Maplet” package in R.

1.Dietary_study_analysis.R
  This will read metabolites from the three tissues (liver, serum, and fecal) and the microbiome data from the source excel files. Measurements of all data were log-transformed and were normalized by quotient normalization.  Data with >20% missing values were excluded from the data set to prevent artifact analysis. 
  The remaining missing values were imputed to the minimum concentration determined for a given metabolite/microbial species. After processing, the data set comprised of 437 metabolites (161 liver, 158 serum, and 118 fecal) and 64 microbial species. 

2.Boxplots_metabolites_microbiome.R

  This will read the the summarised experiment object containing all the metabolites from Liver/Serum/Fecal matrices and the microbiome data. Generate pdf file containing the Log ion counts for all metabolites (Liver, Serum and Fecal) and Log abundance for microbial species. A Bonferroni level of significance of 0.05 / (437 + 64) was used.
