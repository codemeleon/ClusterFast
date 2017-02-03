#!/usr/bin/env python

# This program is used to genarate cluster of highly similar sequences from
# multiple samples of the same species and very closely replated species where
# sequences differeces is much

# LICENSE
# This is free python script/module released in public domain. Users are
# allowed to copy, modify, publish, use, compile or distribute this script.

# This software is provided without any warranty.

# Please cite ... if you use this script for your research purpose

# For more information please refere to github

import os
from os import path, makedirs, getcwd, system
import subprocess
from subprocess import Popen, PIPE, STDOUT
from glob import glob
from shutil import rmtree
from multiprocessing import cpu_count, Pool
import click
import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO
from functools import partial


np.random.seed(676433767)


def is_tool(name):
    """To check the existance of tool in the path."""
    # Code taken from
    # "http://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script"
    try:
        devnull = open(os.devnull)
        Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as err:
        if err.errno == os.errno.ENOENT:
            return False
    return True


def flat_list(lst):
    """Creating list from a nested list."""
    for x in lst:
        if isinstance(x, list):
            for x in flat_list(x):
                yield x
        else:
            yield x


def blat(blatpath, cor, pair):
    """Pairwise blat comparision to search for the best sequences."""
    # Add option for short sequences
    #pair, cor, blatpath = blat_parameters
    if len(pair) == 1:
        return pair[0], []

    sequences = {}
    for p in pair:
        for rec in SeqIO.parse(p, 'fasta'):
            sequences[rec.id] = rec.seq

    grbg = Popen([blatpath, '-threads=%d' % cor,
                  # "-minIdentity=%f" % identity,
                  '-prot', '-noHead', pair[0], pair[1],
                  pair[0].split('.faa')[0]+'.psl'],
                 stdout=PIPE, stderr=STDOUT).communicate()

    try:
        mapped = pd.read_table(pair[0].split('.faa')[0]+'.psl', header=None,
                               usecols=[0, 9, 10, 13, 14])
        mapped = mapped.sort_values(by=[0], ascending=[False])
        reported = set()
        unreported_indexes = []
        for i, row in mapped.iterrows():
            if row[9] in reported or row[13] in reported:
                continue
            unreported_indexes.append(i)
            reported |= set(row[[9, 13]])
        mapped = mapped.ix[unreported_indexes]
        id_pairs = mapped[[9, 13]].values.tolist()
        long_seq_ids = mapped[mapped[10] > mapped[14]][9].tolist()
        long_seq_ids += mapped[mapped[10] <= mapped[14]][13].tolist()
        ids_to_file = set(long_seq_ids) | (set(sequences.keys()) -
                                           set(flat_list(id_pairs)))

        with open(pair[0], "w") as fout:
            for uid in ids_to_file:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        return pair[0], id_pairs
    except:
        with open(pair[0], "w") as fout:
            for uid in sequences:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        return pair[0], []


def makeblastdbf(makeblastdb, infile):
    """Created blast search database for blastp."""
    # system("makeblastdb -in %s -dbtype prot -out %s.bdb" %
    #        (infile, infile.split('.')[0]))

    grbg = Popen([makeblastdb, "-in", infile, "-dbtype", "prot",
                 "-out", "%s.bdb" % infile.split('.')[0]],
                stdout=PIPE, stderr=STDOUT).communicate()
    return


def blastpf(blastp, algo, identity, evalue, infile):
    """Running blast and seleting best searches."""
    # infile, evalue, identity = parameters
    infile_ = infile.split(".")[0]
    # print(infile_)
    for db in glob("tmp/*.faa"):
        if algo == "blast":
            system("%s -db %s.bdb -query %s -outfmt '6 pident qseqid qlen"
                   " sseqid slen length' -evalue %e |"
                   "awk '{if((2*$1*$6/(100.*($3+$5)) >= %f) && ($2 != $4))"
                   " print $2,$4}' >> %s.bst"
                   % (blastp, db.split('.')[0], infile, evalue,
                      identity, infile_))
        elif algo == "min":
            system("%s -db %s.bdb -query %s -outfmt '6 pident qseqid qlen"
                   " sseqid slen length' -evalue %e |"
                   "awk '{if(($1*$6/(100.*($3<=$5?$3:$5)) >= %f) && ($2 != $4))"
                   " print $2,$4}' >> %s.bst"
                   % (blastp, db.split('.')[0], infile, evalue,
                      identity, infile_))
        # elif algo == "max":
        #     system("%s -db %s.bdb -query %s -outfmt '6 pident qseqid qlen"
        #            " sseqid slen length' -evalue %e |"
        #            "awk '{if(($1*$6/(100.*($3>$5?$3:$5)) >= %f) && ($2 != $4))"
        #            " print $2,$4}' >> %s.bst"
        #            % (blastp, db.split('.')[0], infile, evalue,
        #               identity, infile_))
        else: # add information for anmol
            system("%s -db %s.bdb -query %s -outfmt '6 pident qseqid qlen"
                   " sseqid slen length' -evalue %e |"
                   "awk '{if(($1*$6/(100.*($3<=$5?$3:$5)) >= %f) && ($2 != $4))"
                   " print $2,$4}' >> %s.bst"
                   % (blastp, db.split('.')[0], infile, evalue,
                      identity, infile_))


    try:
        data = pd.read_table("%s.bst" % infile_, sep=" ", header=None)
        return data.values
    except:
        return []


def randompairs(n_count):
    """Generates random pairs for the given list of the sequences."""
    """If one sequence in remaining, touple with single file will"""
    """ be genrated."""
    n_list = list(range(n_count))
    np.random.shuffle(n_list)
    pairs = [n_list[i:i+2] for i in range(0, n_count, 2)]
    return pairs


def randomfilepairs(file_list, pairs):
    """Generating random pair of files."""
    file_pairs = []
    for pair in pairs:
        file_pairs.append([file_list[p] for p in pair])
    return file_pairs


def id_arrange_df(sample_ids, sizes, ids):
    """Arranging clustered sequences ids in dataframe."""
    # ids, sample_ids, sizes = parameters
    to_return = {}
    seq_size = []
    seq_count = len(ids)
    samp_count = 0
    for id_ in ids:
        samp, seq = id_.split('___')
        seq_size.append(sizes[id_])
        if samp in to_return:
            to_return[samp][-1] += ",%s:%d" % (id_.split('___')[1], sizes[id_])
        else:
            samp_count += 1
            to_return[samp] = ["%s:%d" % (id_.split('___')[1], sizes[id_])]
    for id_ in set(sample_ids) - set(to_return.keys()):
        to_return[id_] = ['*']
    to_return["samp_count"] = [samp_count]
    to_return["seq_count"] = [seq_count]
    seq_size = np.array(seq_size)
    to_return['std'] = [np.std(seq_size)]
    to_return['mean'] = [np.mean(seq_size)]
    to_return['median'] = [np.median(seq_size)]
    to_return['min'] = [np.min(seq_size)]
    to_return['max'] = [np.max(seq_size)]
    return pd.DataFrame.from_dict(to_return)



def mclf(prog, algo, mcl, dbcreater, mapper, inflation,
         identity, evalue, infile):
    # infile, prog, inflation = parameters
    # prog = 'blat'
    if prog == 'blast':
        infile_ = infile.split(".")[0]
        # evalue = 1e-10
        # identity = 0.25
        grbg = Popen([dbcreater, "-in", infile, "-dbtype", "prot",
                     "-out", "%sdb" % infile],
                     stdout=PIPE, stderr=STDOUT).communicate()
        # system("makeblastdb -in %s -dbtype prot -out %sdb" % (infile, infile))
        if algo == "blast":
            system("%s -db %sdb -query %s -outfmt '6 pident qseqid qlen"
                   " sseqid slen length' -evalue %e |"
                   "awk '{if((2*$1*$6/(100.*($3+$5)) >= %f) && ($2 != $4))"
                   " print $2\"\t\"$4\"\t\"2*$1*$6/(100.*($3+$5))}' > %s.mclin"
                   % (mapper, infile, infile, evalue, identity, infile))
        elif algo == "min":
            system("%s -db %sdb -query %s -outfmt '6 pident qseqid qlen"
                   " sseqid slen length' -evalue %e |"
                   "awk '{if((2*$1*$6/(100.*($3+$5)) >= %f) && ($2 != $4))"
                   " print $2\"\t\"$4\"\t\"$3<=$5?$3:$5}' > %s.mclin"
                   % (mapper, infile, infile, evalue, identity, infile))
        elif algo == "max":
            system("%s -db %sdb -query %s -outfmt '6 pident qseqid qlen"
                   " sseqid slen length' -evalue %e |"
                   "awk '{if((2*$1*$6/(100.*($3+$5)) >= %f) && ($2 != $4))"
                   " print $2\"\t\"$4\"\t\"$3>$5?$3:$5}' > %s.mclin"
                   % (mapper, infile, infile, evalue, identity, infile))
        else:
            system("%s -db %sdb -query %s -outfmt '6 pident qseqid qstart"
                   " qend qlen sseqid sstart send slen length' -evalue %e"
                   " > %s.blst"
                   % (mapper, infile, infile, evalue, infile))

            data = pd.read_table("%s.blst" % infile, header=None)
            data = data[data[1] != data[5]]
            data['q_r'] = data[4] - data[3]
            data['d_r'] = data[8] - data[7]
            data['iden'] = (data[0] * data[9] / 100.) /(
                data[9] + data[['q_r', 'd_r']].max(axis=1) +
                data[[2, 6]].max(axis=1))
            data = data[data['iden'] >= identity]
            data[[1, 5, 'iden']].to_csv("%s.mclin" % infile, header=False,
                                        index=False, sep="\t")


        # system("mcl %s.mclin --abc -I %s  -o %s.mclout" %
        #        (infile, inflation, infile))
    else:
        grbg = Popen([mapper, "-prot", "-noHead", infile, infile,
                      "%s.psl" % infile],
                     stdout=PIPE, stderr=STDOUT).communicate()
        # system("%s -prot -noHead %s %s %s.psl" % (mapper, infile,
        #                                           infile, infile))
        data = pd.read_table("%s.psl" % infile, header=None)
        data = data[data[9] != data[13]]
        if algo == 'blast':
            data["iden"] = data[0] * 1. / data[[10, 14]].mean(axis=1)
        elif algo == 'min':
            data["iden"] = data[0] * 1. / data[[10, 14]].min(axis=1)
        elif algo == 'max':
            data["iden"] = data[0] * 1. / data[[10, 14]].max(axis=1)
        else:
            data['q_r'] = data[10]-data[12]
            data['d_r'] = data[14]-data[16]
            data['iden'] = data[0]/(data[[0, 1, 5, 7]].sum(axis=1) +
                                    data[[11, 15]].max(axis=1) +
                                    data[['q_r', 'd_r']].max(axis=1)
                                    )
            data = data[data['iden'] >= identity]

        data = data[data['iden'] >= identity]
        data[[9, 13, 'iden']].to_csv("%s.mclin" % infile, header=False,
                                     index=False, sep="\t")
    # system("mcl %s.mclin --abc -I 4  -o %s.mclout" %
    #        (infile, infile))
    grbg = Popen([mcl, "%s.mclin" % infile, "--abc", "-I", str(inflation),
                  "-o", "%s.mclout" % infile],
                 stdout=PIPE, stderr=STDOUT).communicate()
    groups = []
    with open("%s.mclout" % infile) as fin:
        for line in fin:
            groups.append(line[:-1].split())
    return groups


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

#/home/devil/Documents/Tools/Clustering/
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--faaf", help="folder containing protein sequence"
              " file from differnt sample", type=str, default="aa")
@click.option("--identity_close", help="Sequence similarity for closely"
              " related sample. Will be used when --distant option is False",
              type=float, default=0.8, show_default=True)
@click.option("--identity_distant", help="Sequence similarity for closely"
              " related sample. Will be used when --distant option is True",
              type=float, default=0.25, show_default=True)
@click.option("--ncor", help="number of cores", type=int,
              default=1, show_default=True)
@click.option("--outfile", help="output cluster file", type=str,
              default='./clusters.clstr', show_default=True)
@click.option("--pblatpath", help="BLAT path", type=str,
              default="pblat", show_default=True)
@click.option("--makeblastdb", help="makeblastdb path. Functional"
              " When distant option is used", type=str, default="makeblastdb",
              show_default=True)
@click.option("--blastp", help="Blastp path. Function when distant is option"
              " in use", type=str, default="blastp", show_default=True)
@click.option("--mcl", help="Run MCL algorithm", type=bool,
              default=False, show_default=True)
@click.option("--mclpath", help="MCL path", type=str,
              default="mcl", show_default=True)
@click.option("--ifl", help="Inflation factor. Active if mcl algorithm will be"
              " used (between 0 and 4)", type=float, default=4.0,
              show_default=True)
@click.option("--distant", help="Samples from distatly replated organisms",
              type=bool, default=False, show_default=True)
@click.option("--evalue", help="evalue for blast search",
              type=float, default=1e-10, show_default=True)
@click.option("--minseq", help="Run MCL of sequence group more than n"
              " sequences", type=str, default=10, show_default=True)
@click.option("--algo", help="Choose similarity calculation algorithm"
              " It will be used in blast search",
              type=click.Choice(['blast', 'anmol', 'min']),
              default='blast', show_default=True)
def run(faaf, identity_close, identity_distant, ncor, outfile, pblatpath,
        evalue, distant, mcl, algo, blastp, ifl, makeblastdb,
        mclpath, minseq):
    """Cluster Generating function."""
    """The program is to genrate protein sequence cluster from sequences
    distributed in different sample files base on their sequence similarities
    reported by BLAT tool. Depending on different parameters different results
    will be provided.
    Each for contain sequence id in format of orgarmism/sample_sequenceid
    Not Trying to find sequence from different organisms which might have
    same fuction Trying to bring similar sequences together only"""
    if not path.isdir(faaf):
        click.echo("Folder path \"%s\" doesn't exist." % faaf)
        exit(0)

    if not is_tool(pblatpath):
        click.echo("pblat is not in the given path %s" % pblatpath)
        exit(0)

    if distant:
        if not is_tool(makeblastdb):
            click.echo("makeblastdb is not in the given path %s" % makeblastdb)
            exit(0)
        if not is_tool(blastp):
            click.echo("blastp is not in the given path %s" % blastp)
            exit(0)

    if mcl:
        if not is_tool(mclpath):
            click.echo("mcl is not in the given path %s" % mcl)
            exit(0)

    if cpu_count() < ncor:
        click.echo("Number of core on the system is less than given in option.\
                   Therefore %d is being used in the analysis" % cpu_count())
        ncor = cpu_count()

    if path.isdir("tmp"):
        rmtree("tmp")
    makedirs("tmp")
    click.echo("Please have patience. It might take a while to finish ...")
    pool = Pool(ncor)
    base_names = []
    seq_sizes = {}
    file_index = {}
    click.echo("Creating temperory copy of orginal file....")
    for i, source in enumerate(glob("%s/*" % faaf)):
        base_name = source.split("/")[-1].split(".")[0]
        base_names.append(base_name)
        file_index[base_name] = i
        with open("%s/tmp/%d.faa" % (getcwd(), i), "w") as fout:
            for rec in SeqIO.parse(source, 'fasta'):
                fout.write(">%s___%s\n%s\n" % (base_name, rec.id, rec.seq))
                seq_sizes["%s___%s" % (base_name, rec.id)] = len(rec.seq)
    click.echo("%d sequence in %d files....." % (len(seq_sizes),
                                                 len(base_names)))

    aa_files = glob("tmp/*")
    graph = nx.Graph()
    aa_files_count = len(aa_files)
    click.echo("Running BLAT to indetify highly similar sequences .....")
    while aa_files_count > 1:
        aa_files.sort()
        file_pairs = randomfilepairs(aa_files, randompairs(aa_files_count))
        if ncor//aa_files_count == 0:
            cor = 1
        else:
            cor = ncor//aa_files_count

        func = partial(blat, pblatpath, cor)
        file_lists_pair = pool.map(func, file_pairs)
        aa_files = []
        for file_, lists in file_lists_pair:
            aa_files.append(file_)
            graph.add_edges_from(lists)
        aa_files_count = len(aa_files)

    sequences = {}
    files_n_seq = {}
    for rec in SeqIO.parse(aa_files[0], 'fasta'):
        sequences[rec.id] = rec.seq

    if not distant:
        identity = identity_close
        click.echo("Running blast to report distantly related sequences")
        grbg = Popen([pblatpath, '-threads=%d' % ncor,
                      '-prot', '-noHead', aa_files[0], aa_files[0],
                      aa_files[0].split('.faa')[0]+'.psl'],
                     stdout=PIPE, stderr=STDOUT).communicate()

        mapped = pd.read_table(aa_files[0].split('.faa')[0]+'.psl',
                               header=None,
                               usecols=[0, 1, 5, 7, 9, 10, 11, 12,
                                        13, 14, 15, 16])
        mapped = mapped[(mapped[9] != mapped[13])]
        # mapped = mapped.sort_values(by=[0], ascending=[False])
        if algo == 'blast':
            mapped = mapped[(mapped[0] * 1. / mapped[[10, 14]].mean(axis=1)) >=
                            identity]
        elif algo == 'min':
            mapped = mapped[(mapped[0] * 1. / mapped[[10, 14]].min(axis=1)) >=
                            identity]
        # elif algo == 'max':
        #     mapped = mapped[(mapped[0] * 1. / mapped[[10, 14]].max(axis=1)) >=
        #                     identity]
        else:
            mapped['q_r'] = mapped[10]-mapped[12]
            mapped['d_r'] = mapped[14]-mapped[16]
            mapped['iden'] = mapped[0]/(mapped[[0, 1, 5, 7]].sum(axis=1) +
                                        mapped[[11, 15]].max(axis=1) +
                                        mapped[['q_r', 'd_r']].max(axis=1)
                                        )
            mapped = mapped[mapped['iden'] >= identity]

        id_pairs = mapped[[9, 13]].values
        graph.add_edges_from(id_pairs)
        del mapped
    if distant:
        click.echo("Running Blast to indetify distantly related sequences")
        identity = identity_distant
        for uid in sequences:
            if uid.split('___')[0] not in files_n_seq:
                files_n_seq[uid.split('___')[0]] = [uid]
            else:
                files_n_seq[uid.split('___')[0]].append(uid)
        for k in files_n_seq:
            with open("tmp/%s.faa" % file_index[k], "w") as fout:
                for sq in files_n_seq[k]:
                    fout.write(">%s\n%s\n" % (sq, sequences[sq]))
        func = partial(makeblastdbf, makeblastdb)
        grbg = pool.map(func, glob("tmp/*.faa"))
        func = partial(blastpf, blastp, algo, identity, evalue)
        pairs = pool.map(func, glob("tmp/*.faa"))
        for pair in pairs:
            graph.add_edges_from(pair)

    connected_ids = list(nx.algorithms.components.connected.
                         connected_components(graph))

    del graph
    system("rm -rf tmp/*")

    if mcl:
        click.echo("Running MCL over grouped sequences")
        tseq = {}
        for fl in glob("%s/*" % faaf):
            bs = path.split(fl)[1].split('.')[0]
            for rec in SeqIO.parse(fl, "fasta"):
                tseq["%s___%s" % (bs, rec.id)] = rec.seq
        connected_idsx = []
        for i, connected_id in enumerate(connected_ids):
            if len(connected_id) < minseq:
                connected_idsx.append(connected_id)
                continue
            with open("tmp/%d.faa" % i, "w") as fout:
                for id_ in connected_id:
                    fout.write(">%s\n%s\n" % (id_, tseq[id_]))
        del tseq, connected_ids
        connected_ids = connected_idsx.copy()
        del connected_idsx
        if distant:
            func = partial(mclf, "blast", algo, mclpath, makeblastdb, blastp,
                           ifl, identity, evalue)
        else:
            func = partial(mclf, "blat", algo, mclpath, None, pblatpath, ifl,
                           identity, None)
        clusters = pool.map(func, glob("tmp/*.faa"))
        for cluster in clusters:
            connected_ids += cluster
    func = partial(id_arrange_df, base_names, seq_sizes)
    dataframes = pool.map(func, connected_ids)
    del connected_ids
    dataframes = pd.concat(dataframes, ignore_index=True)
    dataframes = dataframes.sort_values(by=["samp_count", "seq_count"],
                                        ascending=[False, True])
    dataframes["cluster"] = ["cluster_%d" % d for d in
                             range(1, dataframes.shape[0]+1)]
    col = ['cluster', 'samp_count', 'seq_count', 'min', 'median', 'mean',
           'std', 'max']
    # Rearanging the columns
    col += list(set(dataframes.columns)-set(col))

    dataframes[col].to_csv(outfile, sep="\t", index=False)
    click.echo("Result file %s generated." % outfile)
    click.echo("Finished")
    # del grbg

    # rmtree("tmp")


if __name__ == '__main__':
    run()
