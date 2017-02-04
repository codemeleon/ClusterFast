<H3>Whats is it?</H3>

ClusterFast is a efficiently (memory and CPU) protein sequence clustering pipeline. It is written in python and use PBLAT (multicore BLAT), BLAST and MCL programs.

<H3>External tools</H3>
*Expected to be in system path*

- PBLAT: http://icebert.github.io/pblat/
- NCBI BLAST suit
- mcl: http://www.micans.org/mcl/


<H3>Python and module dependencies</H3>

- Python 3+ (Python2+ hasn't been tested)
- NumPy
- Pandas
- Click
- BioPython
- NetworkX

<H3>Installation</H3>

python setup.py install


<H3>Usage</H3>

clusterfast --faaf < protein_seq_folder > --identity < sequence_similarity > --ncor < #_of_cores_to_use > --outfile < outputfile > --blatpath < blat_absolute_path >
--identity_close <similarity_for_closely_related_sample>
--identity_distant <similarity_for_distant_related_sample>
--blastppath < blast_absolute_path >
--makeblastdb <makeblastdb_path>
--mclpath < mcl_absolute_path >  --sim_algo < blat|anmol|short >
--ifl <Inflation rate for MCL>

* --help/-h : Help
* --faaf: Folder contain only protein fasta files
* --identity_close : Similarity between sequences for closly related samples. Default 0.8
* --identity_distant : Similarity between sequences for distently related samples. Default 0.25
* --ncor: Number of processor to use
* --outfile: Output file path
* --pblatpath: Path for pblat executable. Default: it will consider it in system path
* --pblatpath: Path for pblat executable. Default: it will consider it in system path (pblat)
* --makeblastdb: Path for makeblastdb executable. Default: it will consider it in system path (makeblastdb)
* --blastp: Path for blastp executable. Default: it will consider it in system path (blastp)
* --mclpath: Path for mcl executable. Default: it will consider it in system path (mcl)
* --distant: Are samples distaly related. Default: False
* --mcl: Run MCL. If not used, result will be provided based on connected sequence list. Default: False
* --ifl: Inflation rate for mcl. Default: 4.0
* --minseq: Sequences in connected group required for running MCL. Default: 10
* --algo: four different Identity calculation method. Use in final interation and mcl. Default: blast
  * **blast**: 2*matches/(sum of length of sequences)
  * **anmol**: matches/tolal alignment length as following. **\*** represents matches. Total alignment length includes overhanging sequences, gaps in two sequences, mismatches and matches
  * **min**: matches/length of shortest sequence in the pair
    * ```ADGTHADT--FGGHJJ---DFGDTJHKJLKSDFHKJLJ```
    * ```---*****--******---***-**--******-----```
    * ```---THADTFGFGGHJJSDFDFGFTJKHJLKSDF-----```


<H3>License</H3>

GPLv3

<H3>Requests</H3>
Suggestion and bug reports are welcome.

<H3>Background</H3>

Developed out of frustration of slowness of clustering tools such as orthomcl and proteinortho.
