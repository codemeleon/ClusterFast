# Whats is it?

ClusterFast is a efficiently (memory and CPU) protein sequence clustering pipeline. It is written in python and use PBLAT (multicore BLAT), BLAST and MCL programs.

# External tools

_Expected to be in system path_ - PBLAT: <http://icebert.github.io/pblat/> - NCBI BLAST suit - mcl: <http://www.micans.org/mcl/>

# Python and module dependencies

- Python 3+ (Python2+ hasn't been tested)
- NumPy
- Pandas
- Click
- BioPython
- NetworkX

# Installation

python setup.py install

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

Suggestion and bug reports are welcome.

# Background

Developed out of frustration of slowness of clustering tools such as orthomcl and proteinortho.
