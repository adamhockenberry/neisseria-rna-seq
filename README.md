# neisseria-rna-seq

This repository includes all the code used in a forthcoming manuscript analyzing the transcriptome response of Neisseria gonorrhoeae to oxidative challenge (doi will be provided upon publication).

For the computational side of things, this analysis in the end primarily consisted of a differential expression analysis using Kallisto. Instructions to run Kallisto, and code to create tables of differentially expressed genes from Kallisto/Sleuth output all appear in the `Code` directory. Raw data is not supplied here and users are directed to NCBI GEO GSE114819. A variety of files, however are made available in `Data_release` for those interested in our overall pipeline.

Of note, the dRNA-seq pipeline we ran was primarily manual so code for this portion does not appear here. The relevant `.wig` files (available at GEO) were visualized in IGV and annotation of peaks was performed manually.

-Adam J Hockenberry
