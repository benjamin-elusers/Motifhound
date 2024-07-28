# MotifHound: De Novo Protein Motif Discovery Algorithm

![License](https://img.shields.io/badge/License-MIT-blue.svg)
![Version](https://img.shields.io/badge/version-1.0.0-green.svg)

## Cite us

>Kelil A, Dubreuil B, Levy ED, Michnick SW (2014) 
>Fast and Accurate Discovery of Degenerate Linear Motifs in Protein Sequences. 
>PLOS ONE 9(9): e106081. https://doi.org/10.1371/journal.pone.0106081

## Abstract

Linear motifs mediate a wide variety of cellular functions, which makes their characterization in protein sequences crucial to understanding cellular systems. However, the short length and degenerate nature of linear motifs make their discovery a difficult problem. Here, we introduce MotifHound, an algorithm particularly suited for the discovery of small and degenerate linear motifs. MotifHound performs an exact and exhaustive enumeration of all motifs present in proteins of interest, including all of their degenerate forms, and scores the overrepresentation of each motif based on its occurrence in proteins of interest relative to a background (e.g., proteome) using the hypergeometric distribution. To assess MotifHound, we benchmarked it together with state-of-the-art algorithms. The benchmark consists of 11,880 sets of proteins from S. cerevisiae; in each set, we artificially spiked-in one motif varying in terms of three key parameters, (i) number of occurrences, (ii) length and (iii) the number of degenerate or “wildcard” positions. The benchmark enabled the evaluation of the impact of these three properties on the performance of the different algorithms. The results showed that MotifHound and SLiMFinder were the most accurate in detecting degenerate linear motifs. Interestingly, MotifHound was 15 to 20 times faster at comparable accuracy and performed best in the discovery of highly degenerate motifs. We complemented the benchmark by an analysis of proteins experimentally shown to bind the FUS1 SH3 domain from S. cerevisiae. Using the full-length protein partners as sole information, MotifHound recapitulated most experimentally determined motifs binding to the FUS1 SH3 domain. Moreover, these motifs exhibited properties typical of SH3 binding peptides, e.g., high intrinsic disorder and evolutionary conservation, despite the fact that none of these properties were used as prior information. MotifHound is available (http://michnick.bcm.umontreal.ca or http://tinyurl.com/motifhound) together with the benchmark that can be used as a reference to assess future developments in motif discovery.

## Summary

MotifHound is a state-of-the-art, de novo protein motif discovery algorithm. It is designed to identify and characterize significant motifs in protein sequences. MotifHound is well-suited for large-scale analysis, featuring high efficiency and accuracy.

## Description

MotifHound is an algorithm for the discovery of degenerate linear motifs encoded in protein sequences. MotifHound is based on the fundamental premise that, linear motifs mediating a particular function are enriched in proteins exhibiting that function and rare or absent in other proteins. Accordingly, MotifHound exhaustively enumerates all possible motifs present in proteins of interest, including all degenerated forms, and scores the biological relevance of each motif by evaluating its enrichment relative to the proteome using the cumulative hypergeometric distribution. The package includes also a benchmark of 11,880 protein sets from S. cerevisiae in which we artificially spiked-in motifs varying in terms of three key parameters, (i) occurrences, (ii) length, (iii) and number of degenerate positions. The benchmark enabled to evaluate the impact of these three properties on the performance of MotifHound.

## Features

- Efficient and accurate motif discovery in protein sequences
- Scalable to large datasets
- User-friendly and customizable parameters for advanced users
- Rich and detailed output for downstream analysis
- Open-source and actively maintained

## Installation

### 1.  Download

The original version of MotifHound (as of publication) is available from [Gdrive](https://drive.google.com/drive/folders/1AJKMQAZy9UEbpcGYiWTXcVtRFf5mHcgb?usp=drive_link)

You can download the original tar file [here](https://github.com/benjamin-elusers/Motifhound/raw/main/archive/MotifHound_171013.tar.gz). After downloading, the file needs to be decompressed:

> tar -xvzf MotifHound_171013.tar.gz
> cd MotifHound/

The updated version is available on this current [GitHub repository]([https://github.com/benjamin-elusers/Motifhound]).

### 2. Third party libraries
In order to work with S. cerevisiae and H. sapiens datasets, we also provide precomputed data on disorder, Pfam domains, evolution and function descriptions. These data are not required to run MotifHound but they are recommended to use it with its full potential. These data can be downloaded [here](https://drive.google.com/file/d/1vRcIxc9KJxYCkOmj4xPVEQyG9otnvwf4/view?usp=drive_link) and copied in the "data" directory of MotifHound.

Importantly, MotifHound uses the following programs/libraries that need to be installed:
"Judy arrays" is a C library that can be either installed with the following command on Ubuntu systems:

> sudo apt-get install libjudydebian1

Alternatively, the source of this library can be downloaded [here](http://judy.sourceforge.net/).
blast is required if masking of homologous regions is desired. It can be installed with the following command on Ubuntu systems:

> sudo apt-get install blast2

Alternatively, please refer to the instruction from [NCBI blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi) for additional instructions.

Some Perl modules are also necessary, for this we recommend to use the following commands:
> sudo apt-get install perl cpanminus

followed by:

> sudo cpanm Tk Tk::TableMatrix File::Basename File::Copy List::MoreUtils Cwd Getopt::Long FindBin Benchmark

### 3. Usage

Following these installs, you can then run MotifHound:

Here is a basic example of how to run MotifHound:

> perl ./Scripts/MotifHound.pl --Setfile ./Data/Seq/Set/YEAST_Set_TEST.faa \ 
                               --Proteome ./Data/Seq/Proteome/YEAST_Proteome_TEST.fasta \
                               --Size 3 10 \ 
                               --Scan \ 
                               --WD ./Results \
                               --H --D \
                               --Blastfile ./Data/Blast/Blast_YEAST_Proteome.blast \
                               --Disofile ./Data/Disorder/YEAST_Proteome_DISORDER.dat \ 
                               --Pfam_annot ./Data/Domains/YEAST_Proteome_Pfam_Domains.txt \
                               --Gene_annot ./Data/Genes/YEAST.data \
                               --HTML

To display the help :

> perl ./Scripts/MotifHound.pl --help

### 4. Benchmark data

We benchmarked MotifHound by creating datasets of protein sequences from S. cerevisiae, in which we spiked-in known motifs. The motifs spiked-in vary in length, number of defined positions, and number of repeats. To exhautively cover these three parameters combinations, we created 11,880 datasets. These datasets are available for download, which we hope will help in the development of future algorithms for motif discovery.


## Parallel version

DALEL is a more recent algorithm based on the same principles as MotifHound. 
DALEL can be tested online as a [webserver](http://michnick.bcm.umontreal.ca/dalel/Server/Index.aspx).
The source code is also available for Windows.

It can also be adapted if you have a negative and positive sets of sequence.
DALEL inspects degenerate positions to find the combination of most frequent amino-acids.

To cite DALEL:
>Kelil A, Dubreuil B, Levy ED, Michnick SW. 
>Exhaustive search of linear information encoding protein-peptide recognition. 
>PLoS Comput Biol. 2017 Apr 20;13(4):e1005499. doi: 10.1371/journal.pcbi.1005499. PMID: 28426660; PMCID: PMC5417721.

## License

MotifHound is released under the [MIT License](https://opensource.org/licenses/MIT).

## Acknowledgements

We would like to thank all contributors and users of MotifHound. Your feedback and support make MotifHound a better tool for the scientific community.

## Contact

If you have questions, feature requests, or bug reports, please open an issue on our [GitHub issue tracker](https://github.com/benjamin-elusers/MotifHound/issues). 

For other inquiries, feel free to contact us at: `benjamin.dubreuil@weizmann.ac.il`.

---

Please note that the 'username' in the URLs above should be replaced with your actual GitHub username or the username of the repository owner. The email address should also be replaced with your actual contact email or the contact email for the project.
