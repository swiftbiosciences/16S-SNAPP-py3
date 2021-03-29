16S SNAPP is an analysis workflow to be run on Linux/Mac command line
interface (CLI) for 16S multi-V region amplicon sequences from Swift's 16S SNAP
panel. It generates taxonomic composition tables at and above genus ranks from
multiple demultiplexed fastq files.

16S SNAPP's approach is to associate sequence reads derived from multiple
amplicon regions to their most probable sequences of origin, i.e. the assumed
templates. This is done through database search (blastn) for high identity
matches among reference sequences followd by intersecting aligned reads on the
matching references, read count allocation for multi-mapped reads, and
classification of consensus sequences. It offers higher sensitivity for 16S
multiple-amplicon data compared to tools designed for single 16S amplicon data.

Setup:
   1. Install and setup the following software:
      R (https://www.r-project.org/),
      Java ≥1.8.0_131,
      Python 3.6.8 or later (with Numpy ≥ 1.16.0, Pandas ≥ 1.1.4, Scipy ≥1.4.1),
      DADA2 (https://benjjneb.github.io/dada2/),
      BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/),
      RDPTools (github.com/rdpstaff/RDPTools),
      Cutadapt (https://cutadapt.readthedocs.io/en/stable/),
      VSEARCH (https://github.com/torognes/vsearch),
      MAFFT (https://mafft.cbrc.jp/alignment/software/,
      FASTTREE (http://www.microbesonline.org/fasttree/)
   3. Create a folder , e.g. 'DB', and download to it the reference and primer 
      files from (https://ws.onehub.com/folders/82s4teyk).
   4. Clone this repository (git clone https://github.com/swiftbiosciences/16S-SNAPP-py3)
   5. Edit “config.txt” to enter absolute paths to tools, 'DB' and primer file, 
      and the expected single read length after primer is trimmed


Command to run: snapp.sh config.txt inputdir workdir

   Input: gzipped post-demultiplexing sequence fastq files each (pair)
          representing a single sample and each carrying a Swift 16S primers
          on its 5’ end.

  Output: (“RESDIR” folder under the directory where the command is run):

          lineage-table.tsv,
          feature-table.tsv,
          taxonomy-table.tsv,
          templates_mafft.tree

Limitations:
1. The workflow largely relies on aligning reads from different gene regions
   to their closest reference sequences from the database (currently packed
   with RDP11.5). Nevertheless, sequence pairs that can’t be aligned to any
   reference sequences within the set identity cutoff are treated as individual
   features and classified directly.
2. The resolution of taxonomic assignments depends on the classifier, mainly its
   taxonomic coverage of the samples, and the sequencing depth. The current
   version of 16S SNAPP uses RDP Classifier, but is potentially compatible with 
   other k-mer's based classifier.
3. The phylogenetic tree is built from the sequence alignment of the assumed
   template sequences to approximate the phylogeny of multiple amplicons that
   inherently lack comparable alignment positions.

Additional notes:
   Quality trimming/filtering/denoising may be adjusted in run_dada2.R script
   as needed for specific data sets to improve the performance of the analysis.
