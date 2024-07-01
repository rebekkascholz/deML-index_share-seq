# deML-index_share-seq

Includes R script (deML-index_share-seq.R) to create an index file for SHARE-seq sequencing data (fastq files) demultiplexing with deML (https://github.com/grenaud/deML).
Published in [Scholz et al., 2024](https://doi.org/10.1007/978-1-0716-3437-0_35).

To create the index file, create a .yaml file with the same structure as the example uploaded here (config_SHARE.yaml) for your SHARE-seq runs. 
Include the first and the last Round1 barcode for each sample (= project) as well as up to three P5 (Ad1.X) indices (= up to three different libraries per sample).

Download the R script and start in the directory of the R script via command line:

<code>Rscript deML-index_share-seq.R config.SHARE.yaml</code>

The index.txt file will be saved to the current directory.
