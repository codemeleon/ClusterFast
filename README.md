# About

ClusterFast is a scalable, rapid and (CPU and memory) efficient command-line interface for the clustering of orthologous proteins encoded by multiple genomes. ClusterFast was developed by the MiDEP group at the [Malawi Liverpool Wellcome Trust Clinical Research Programme](http://www.mlw.medcol.mw/) (members of [H3ABionet](http://www.h3abionet.org/)), after they became frustrated by the long run time and large memory requirements of currently avaliable clustering programs. 

ClusterFast works by first grouping the most similar protein sequences encoded between a random pair of genomes (files) within the input dataset, then choosing the longest sequence for further comparison to protein sequences selected from another pair of genomes until only one file remains. This is then used to identify less similar or paralogous sequences. The tool is suitable for use with both prokaryotic and eukaryotic genomes. By employing this novel approach, ClusterFast substantially reduces the memory and processing time required for the clustering of orthologous proteins compared to other avaliable programs. 

Using a test dataset of 140 pneumococcal genomes (each 2.5Mbp in size and encoding ~1500 genes), ClusterFast successfully executed in ~5 minutes on a single core. 

ClusterFast is written in Python and uses PBLAT (multicore BLAT), BLAST and MCL programs.

# External tools

_Expected to be in system path_ - PBLAT: <http://icebert.github.io/pblat/> - NCBI BLAST suit - mcl: <http://www.micans.org/mcl/>

# Python and module dependencies

There are a number of dependencies required for ClusterFast, with instructions specific to the type of system you have:

- Python 3+ (Python2+ hasn't been tested)
- NumPy
- Pandas
- Click
- BioPython
- NetworkX

# Installation

python setup.py install

If the installation fails, please contact your system administrator. If you discover any bugs, please let us know by emailing anmol@liv.ac.uk

# Input Files

The input format for ClusterFast is protein sequence files (.faa) of translated amino acid sequences of predicted open reading frames for each genome (sample) in the input dataset. These files can be created using [Prokka](https://github.com/tseemann/prokka). 

# Usage

clusterfast --faaf < protein_seq_folder > --identity < sequence_similarity > --ncor < #_of_cores_to_use > --outfile < outputfile > --blatpath < blat_absolute_path > --identity_close <similarity_for_closely_related_sample> --identity_distant <similarity_for_distant_related_sample> --blastppath <blast_absolute_path> --makeblastdb <makeblastdb_path> --mclpath <mcl_absolute_path>  --sim_algo <blat|anmol|short> --ifl <inflation_rate_for_mcl> --makeblastdb <makeblastdb_path> --identity_distant <similarity_for_distant_related_sample> --identity_close <similarity_for_closely_related_sample> --keepclean

- --help/-h : Help
- --faaf: Folder contain only protein fasta files
- --identity_close : Similarity between sequences for closly related samples. Default 0.8
- --identity_distant : Similarity between sequences for distently related samples. Default 0.25
- --ncor: Number of processor to use
- --outfile: Output file path
- --pblatpath: Path for pblat executable. Default: it will consider it in system path
- --pblatpath: Path for pblat executable. Default: it will consider it in system path (pblat)
- --makeblastdb: Path for makeblastdb executable. Default: it will consider it in system path (makeblastdb)
- --blastp: Path for blastp executable. Default: it will consider it in system path (blastp)
- --mclpath: Path for mcl executable. Default: it will consider it in system path (mcl)
- --distant: Are samples distaly related. Default: False
- --mcl: Run MCL. If not used, result will be provided based on connected sequence list. Default: False
- --ifl: Inflation rate for mcl. Default: 4.0
- --minseq: Sequences in connected group required for running MCL. Default: 1
- --keepclean: Keep deleting interemediate files to minimise disk usage
- --algo: four different Identity calculation method. Use in final interation and mcl. Default: blast

  - **blast**: 2*matches/(sum of length of sequences)
  - **min**: matches/length of shortest sequence in the pair
  - **anmol**: matches/tolal alignment length as following.

    - "*" represents matches in the alignment. Total alignment length includes overhanging sequences, gaps in two sequences, mismatches and matches
    - `ADGTHADT--FGGHJJ---DFGDTJHKJLKSDFHKJLJ`
    - `---*****--******---***-**--******-----`
    - `---THADTFGFGGHJJSDFDFGFTJKHJLKSDF-----`

# License

GPLv3

# Requests

Suggestions for improvements are welcome.

