# Project
Alexandros Kavathas
Hello 
Galaxy
Galaxy was used in order to check the t-cell acute lymphoblastic leukemia or verboom with 4 control and 4 Roels single positive cells.
Verboom: https://haematologica.org/article/view/8715SRA project: SRP132968SRA samples: SRR6727631, SRR6727635, SRR6727639, SRR6727666T-ALL samples, 4 male TAL subtypes
Roels: https://www.nature.com/articles/s41590-020-0747-9SRA project: SRP263016SRA samples:  SRR11832951, SRR11832952, SRR11832962, SRR11832963Normal single-positive CD4 and CD8 cells, two of each


In Galaxy we 
1 add the data in
2 run the data through the fastQ and extract the reads
3 we run the fastQC with the data from fastQ 
4 then we run the fastp process of the fastQ data
5 we run the fast and sensitive alignment program called hisat2 from fastQ data
6 multiQC with the hisat2 data
7 future counts with hisat 2 data 
8 edgeR perform expression of count data with feature counts
9 annotatemyids using the edgeR data
R code
So in R code we read the data that we downloaded from annotatemyids from Galaxy.
We format the data by merging the two files
then we created the metadata file
We create a pca plot 
We make a two sided t-test
Create the log2fc
We visualize the data in three plots 1 kegg 2 data plot 3 scatter plot
