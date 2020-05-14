SNAP APP is an analysis workflow to be run on Linux/Mac command line
interface (CLI) for 16S amplicon sequences from Swift 16S panel. It
generates taxonomic composition tables at and above genus ranks from multiple
demultiplexed fastq files.

SNAP APP's approach is to associate sequence reads derived from multiple
amplicon regions to their most probable sequences of origin, the assumed
templates. This is done through database search (blastn) and intersection of
aligned reads, read count allocation, and classification of consensus sequences.
It offers higher sensitivity compared to tools designed for the common read-level
analysis used for single 16S amplicon data.

Setup:
   1. Install and setup the following software:
      Java 1.8.0_131 or later,
      VSEARCH (https://github.com/torognes/vsearch),
      RDPTools (github.com/rdpstaff/RDPTools),
      BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/),
      Python 2.7,
      Numpy 1.16.2 or later,
      Pandas 0.24.2 or later,
      Scipy 1.2.1 or later
   2. Unzip the package file to obtain this file, config.txt, and two folders:
      “DB” and ”scripts”.
   3. Edit “config.txt” to enter absolute paths to each tool, primer and
      database files, and to ”SNAP-APP” folder. Data processing parameters,
      e.g. expected single read length after primer is trimmed and classification
      confidence cutoff, should be adjusted when needed.

Command to run: snapp.sh config.txt inputdir workdir
   Input: gzipped post-demultiplexing sequence fastq files each (pair)
          representing a single sample and each carrying a Swift 16S primers
          on its 5’ end.
   Output: (“RESDIR” folder under the directory where the command is run):
          lineage-table.tsv
          feature-table.tsv
          taxonomy-table.tsv
          templates_mafft.tree

Limitations:
1. The workflow largely relies on aligning reads from different gene regions
   to their closest reference sequences from the database (currently packed
   with RDP11.5). Nevertheless, sequence pairs that can’t be aligned to any
   reference sequences within the set identity cutoff are treated as individual
   features and classified directly.
2. The resolution of taxonomic assignments depends on the classifier, mainly its
   taxonomic coverage of the samples, and the sequencing depth. The current
   version of SNAP APP uses RDP Classifier, but is potentially compatible with 
   other k-mer's based classifier.
3. The phylogenetic tree is built from the sequence alignment of the assumed
   template sequences to approximate the phylogeny of multiple amplicons that
   inherently lack comparable alignment positions.

Additional notes:
   More trimming/filtering/denoising may be added as needed for specific data
   sets to improve the performance of the analysis.
