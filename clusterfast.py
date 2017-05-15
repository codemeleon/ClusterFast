##!/usr/bin/env python

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
from os import path, makedirs, getcwd, system, stat
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
import time
import itertools as its
from scipy.optimize import curve_fit
pd.options.mode.chained_assignment = None


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


def blat(algo, blatpath, cor, keepclean, beginning, identity, minlen,
         mindiff, minmap, mclinfile, pair):
    """Pairwise blat comparision to search for the best sequences."""
    if beginning:
        for i, fl in enumerate(pair):
            flname = path.split(fl)[1]
            base_name = flname.split(".")[0]
            pair[i] = "tmp/%s" % flname
            with open(pair[i], "w") as fout:
                for rec in SeqIO.parse(fl, 'fasta'):
                    if len(rec.seq) < minlen:
                        continue
                    fout.write(">%s___%s___%d\n%s\n" % (base_name, rec.id,
                                                        len(rec.seq), rec.seq))
    if len(pair) == 1:
        return pair[0], []

    sequences = {}
    # TODO: Need to make work only once in case of only one sequence file
    for p in pair:
        for rec in SeqIO.parse(p, 'fasta'):
            sequences[rec.id] = rec.seq
    # print(pair)
    base_name = pair[0].split('.')[0]
    Popen([blatpath, '-threads=%d' % cor,
           '-prot', '-noHead', pair[0], pair[1],
           '%s.psl' % base_name],
          stdout=PIPE, stderr=STDOUT).communicate()
    # exit(1)

    if stat('%s.psl' % base_name).st_size:
        mapped = pd.read_table('%s.psl' % base_name, header=None,
                               usecols=[0, 1, 5, 7, 9, 10, 11, 12,
                                        13, 14, 15, 16])
        # TODO:Fix empty file problem
        # Clsuter and choose
        # add a marker here to avoid researching and ask to generate only bst
        # files for the analysis
        if mindiff:
            mapped = mapped[(mapped[[10, 14]].min(axis=1)) >
                            mindiff * (mapped[[10, 14]].max(axis=1))]
        if minmap:
            # print(mapped.head( ))
            mapped["minmap"] = mapped.apply(lambda x:
                                            min([x[12]-x[11],
                                                 x[16]-x[15]]) /
                                            max([x[10], x[14]]), axis=1)
            mapped = mapped[mapped["minmap"] > minmap]
            del mapped["minmap"]
        if algo == "min":
            mapped.loc[:, "identity"] = mapped[0] / mapped[[10, 14]].min(axis=1)
        elif algo == "max":
            mapped.loc[:, "identity"] = mapped[0] / mapped[[10, 14]].max(axis=1)
        elif algo == "blast":
            mapped.loc[:, "identity"] = mapped[0] / mapped[[10, 14]].mean(axis=1)
        else:
            mapped.loc[:, 'q_r'] = mapped[10] - mapped[12]
            mapped.loc[:, 'd_r'] = mapped[14] - mapped[16]
            mapped.loc[:, "identity"] = mapped[0]/(mapped[[0, 1, 5, 7]].sum(axis=1) +
                                            mapped[[11, 15]].max(axis=1) +
                                            mapped[['q_r', 'd_r']].max(axis=1)
                                            )
        mapped = mapped[mapped["identity"] > identity][[9, 13, "identity"]]
        if mclinfile:
            mapped = mapped[mapped[9] != mapped[13]]
            # TODO: find what kind of output you expect here
            mapped.to_csv("tmp/%s.mclin" % path.split(pair[0])[1].split(".faa")[0], header=False,
                          index=False, sep="\t")
            return
            #f returning some file name here for simplicity

        mapped = mapped[[9, 13]].values.tolist()
        connected_ids = nx.Graph()
        connected_ids.add_edges_from(mapped)
        nxacc = nx.algorithms.components.connected
        connected_ids = list(nxacc.connected_components(connected_ids))
        selected_ref = []
        for connected_id in connected_ids:
            sz = 0
            t_id = ""
            for _id in connected_id:
                sq_sz = len(sequences[_id])
                if sq_sz > sz:
                    sz = sq_sz
                    t_id = _id
            selected_ref.append(t_id)
        # print(list(flat_list(mapped)))
        ids_to_file = set(selected_ref) | (set(sequences.keys()) -
                                           set(flat_list(mapped)))

        with open(pair[0], "w") as fout:
            for uid in ids_to_file:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        if keepclean:
            system("rm %s %s.psl" % (pair[1], base_name))
        return pair[0], mapped  # id_pairs
    else:  # Empty file
        with open(pair[0], "w") as fout:
            for uid in sequences:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        if keepclean:
            system("rm %s %s.psl" % (pair[1], base_name))
        return pair[0], []


def blast_dataframe(mapped2, mindiff, minmap, algo):
    mapped = mapped2.copy()
    mapped = mapped[mapped[0] != mapped[1]]
    mapped.loc[:, "qsize"] = mapped[0].map(lambda x: int(x.split("___")[2]))
    mapped.loc[:, "ssize"] = mapped[1].map(lambda x: int(x.split("___")[2]))
    if minmap:
        mapped["minmap"] = mapped.apply(lambda x:
                                        min([x[7]-x[6],
                                             x[9]-x[8]]) /
                                        max([x['qsize'], x['ssize']]), axis=1)
        mapped = mapped[mapped["minmap"] > minmap]
        del mapped["minmap"]
    if mindiff:
        mapped = mapped[(mapped[['qsize', 'ssize']].min(axis=1)) >
                        (mindiff * mapped[['qsize', 'ssize']].max(axis=1))]

    if algo == "min":
        mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] * 2 / 100.) /
                              mapped[["qsize", "ssize"]].min(axis=1))
    elif algo == "max":
        mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] * 2 / 100.) /
                              mapped[["qsize", "ssize"]].max(axis=1))
    elif algo == "blast":
        mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] / 100.) /
                              mapped[["qsize", "ssize"]].mean(axis=1))
    else:
        mapped.loc[:, 'q_r'] = mapped["qsize"] - mapped[7]
        mapped.loc[:, 'd_r'] = mapped["ssize"] - mapped[9]
        mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] * 2 / 100.)/(
            mapped[3] + mapped[['q_r', 'd_r']].max(axis=1) +
            mapped[[6, 8]].max(axis=1)))
    return mapped
#
# def blast_dataframe(mapped, mindiff, minmap, algo):
#     mapped = mapped[mapped[0] != mapped[1]]
#     mapped.loc[:, "qsize"] = mapped[0].map(lambda x: int(x.split("___")[2]))
#     mapped.loc[:, "ssize"] = mapped[1].map(lambda x: int(x.split("___")[2]))
#     if minmap:
#         mapped.loc[:, "minmap"] = mapped.apply(lambda x:
#                                         max([x[7]-x[6],
#                                              x[9]-x[8]]) /
#                                         max([x['qsize', 'ssize']]))
#         mapped = mapped[mapped["minmap"] > minmap]
#         del mapped["minmap"]
#     if mindiff:
#         mapped = mapped[(mapped[['qsize', 'ssize']].min(axis=1)) >
#                         (mapped[['qsize', 'ssize']].max(axis=1))]
#
#     if algo == "min":
#         mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] * 2 / 100.) /
#                               mapped[["qsize", "ssize"]].min(axis=1))
#     elif algo == "max":
#         mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] * 2 / 100.) /
#                               mapped[["qsize", "ssize"]].max(axis=1))
#     elif algo == "blast":
#         mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] / 100.) /
#                               mapped[["qsize", "ssize"]].mean(axis=1))
#     else:
#         mapped.loc[:, 'q_r'] = mapped["qsize"] - mapped[7]
#         mapped.loc[:, 'd_r'] = mapped["ssize"] - mapped[9]
#         mapped.loc[:, "identity"] = ((mapped[3]*mapped[2] * 2 / 100.)/(
#             mapped[3] + mapped[['q_r', 'd_r']].max(axis=1) +
#             mapped[[6, 8]].max(axis=1)))
#     return mapped


def makeblastdbf(infile):
    """Created blast search database for blastp."""

    Popen(["makeblastdb", "-in", infile, "-dbtype", "prot",
           "-out", "%s.bdb" % infile.split('.')[0]],
          stdout=PIPE, stderr=STDOUT).communicate()
    return


def blastpfxxx(algo, identity, evalue, keepclean, mindiff, minmap,
            mclinfile, infile):
    """Running BLAST and seleting best searches."""
    makeblastdbf(infile)
    infile_ = infile.split(".")[0]
    # TODO: Fix things here
    if mclinfile:
        infile_bs = path.split(infile)[1].split(".fa")[0]
        system("blastp -db %s.bdb -query %s -evalue %e"
               " -outfmt 6 > %s.bst" % (infile_, infile,
                                        evalue, infile_))
        # print("%s.bst" % infile_)
        mapped = pd.read_table("%s.bst" % infile_, header=None,
                               usecols=range(10))
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped[[0, 1, 'identity']].to_csv("tmp/%s.mclin" % infile_bs, header=False,
                      index=False, sep="\t")  # Fix it
        return "tmp/%s.mclin" % infile_bs
    # ncor = 1
    for query in glob("tmp/*.faa"):  # TODO: Need to flip the situation
        query_id = path.split(query)[1].split(".faa")[0]
        outfilebs = "%s_%s" % (infile_, query_id)
        # print(outfilebs)
        system("blastp -db %s.bdb -query %s -evalue %e"
               " -outfmt 6 > %s.bst" % (infile_, query,
                                        evalue, outfilebs))

        mapped = pd.read_table("%s.bst" % outfilebs, header=None,
                               usecols=range(10))
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped[[0, 1, 'identity']].to_csv("%s.bst" % outfilebs, header=None,
                                          index=None, sep="\t")

    return  # mapped.values

def norm_value(df):
    def curve(x, a, b):
        return a*x + b
    #For each multi, find the maximum bit score
    #And normalise according to that
    df["multi"] = df["qsize"]*df["ssize"] # No identical ids
    for_curve = df[["multi", 11]]
    for_curve = for_curve.sort_values(["multi", 11], ascending=[False, False])
    for_curve = for_curve.drop_duplicates(["multi"])
    for_curve["multi"] = np.log10(for_curve["multi"])
    for_curve[11] = np.log10(for_curve[11])
    max_11 = max(for_curve[11])
    max_multi = min(list(for_curve.loc[for_curve[11] == max_11, "multi"]))
    for_curve = for_curve[~((for_curve[11] < max_11)&(for_curve["multi"] > max_multi))]
    multi_min, multi_max = np.ceil(np.min(for_curve["multi"])), np.floor(np.max(for_curve["multi"]))
    diff = multi_max -  multi_min
    if diff < 1:
        diff = 1
    number_of_fragments = 10**diff
    x = []
    y = []

    linespace = np.linspace(multi_min, multi_max, number_of_fragments)
    for i in range(len(linespace)-1): #Make it multicore
        tmp_curve = for_curve[(for_curve["multi"]>=linespace[i]) &
                              (for_curve["multi"]<linespace[i+1])]
        if not len( tmp_curve):
            continue
        max_11 = max(tmp_curve[11])
        max_multi = max(list(tmp_curve.loc[tmp_curve[11] == max_11, "multi"]))
        x.append(max_multi)
        y.append(max_11)


    ## Generate list array
    print(x, y, "Kiran")# If empty, leave the cluster
    pars, _ =  curve_fit(curve, x, y)
    df["identity"] = df[11] /((df["multi"]**pars[0])* (10**pars[1]))
    df.loc[df["identity"]>1.0, "identity"] = 1.0
    return df

def blastpf(algo, identity, evalue, keepclean, mindiff, minmap,
            mclinfile, db_query):
    """Running BLAST and seleting best searches."""
    db, infile = db_query
    infile_ = infile.split(".")[0]
    # TODO: Fix things here
    if mclinfile:
        makeblastdbf(infile)
        # print("anmol")
        infile_bs = path.split(infile)[1].split(".fa")[0]
        system("blastp -db %s.bdb -query %s -evalue %e"
               " -outfmt 6 > %s.bst" % (infile_, infile,
                                        evalue, infile_))
        # print("%s.bst" % infile_)
        if stat("%s.bst" % infile_).st_size:
            mapped = pd.read_table("%s.bst" % infile_, header=None)
            mapped = blast_dataframe(mapped, mindiff, minmap, algo)
            mapped = mapped[mapped["identity"] >= identity]
            # Need some changes here
            mapped = mapped.sort_values([10, 11, "identity"], ascending=[True, False, False])
            mapped = mapped.drop_duplicates([0, 1])
            mapped = norm_value(mapped[[0, 1, "qsize", "ssize", 10, 11]])
            # Remove duplicates, base high percent, low evalue, high bitscore
            # Solve some problem here
            mapped[[0, 1, 'identity']].to_csv("tmp/%s.mclin" % infile_bs, header=False,
                          index=False, sep="\t")  # Fix it
        else:
            open("tmp/%s.mclin" % infile_bs, 'a').close()
        return "tmp/%s.mclin" % infile_bs
    # ncor = 1

    db_ = db.split(".")[0]
    query_id = path.split(infile)[1].split(".faa")[0]
    outfilebs = "%s_%s" % (db_, query_id)
        # print(outfilebs)
    system("blastp -db %s.bdb -query %s -evalue %e"
           " -outfmt 6 > %s.bst" % (db_, infile,
                                    evalue, outfilebs))
    if stat("%s.bst" % outfilebs).st_size:
        mapped = pd.read_table("%s.bst" % outfilebs, header=None)
        print("Anmol")
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        print("Anmol2")
        mapped = mapped[mapped["identity"] >= identity]
        print("Anmol3")
        mapped = mapped.sort_values([10, 11, "identity"], ascending=[True, False, False])
        print("Anmol4")
        mapped = mapped.drop_duplicates([0, 1])
        print("Anmol5")
        mapped[[0, 1, "qsize", "ssize", 10, 11]].to_csv("%s.bst" % outfilebs, header=None,
                                          index=None, sep="\t") # Fix here
    else:
        pass

    return  # mapped.values


# def split_para():



def randompairs(n_count):
    """Generates random pairs for the given list of the sequences."""
    """If one sequence in remaining, touple with single file will"""
    """ be genrated."""
    n_list = list(range(n_count))
    np.random.shuffle(n_list)

    for i in range(0, n_count, 2):
        yield n_list[i:i+2]


def randomfilepairs(file_list, pairs):
    """Generating random pair of files."""
    file_pairs = []
    for pair in pairs:
        yield [file_list[p] for p in pair]


def id_arrange_df(sample_ids, ids):  # sizes,
    """Arranging clustered sequences ids in dataframe."""

    to_return = {}
    seq_size = []
    seq_count = len(ids)
    samp_count = 0
    for id_ in ids:
        print(id_)
        samp, seq, sz = id_.split('___')  # Don't do this. Do it after everything is finished
        # TODO: Put only names separated by comma. nothing else
        sz = int(sz)
        seq_size.append(sz)
        if samp in to_return:
            to_return[samp][-1] += ",%s:%d" % (seq, sz)
        else:
            samp_count += 1
            to_return[samp] = ["%s:%d" % (seq, sz)]
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





def mclf(inflation, ncor, infile, distant):
    if infile:
        Popen(["mcl", infile, "--abc", "-I", str(inflation),
               "-o", "tmp/temp.mclout", "-te", str(ncor)],
              stdout=PIPE, stderr=STDOUT).communicate()
    else:
        system("cat tmp/*.bst > tmp/temp.mclin")
        if distant:
            mclin = pd.read_table("tmp/temp.mclin", header=None)
            mclin = mclin.rename(columns={2:"qsize", 3:"ssize", 4:10, 5:11})
            mclin = norm_value(mclin)
            mclin[[0, 1, "identity"]].to_csv("tmp/temp.mclin", header=False, index=False, sep="\t")

        Popen(["mcl", "tmp/temp.mclin", "--abc", "-I", str(inflation),
               "-o", "tmp/temp.mclout", "-te", str(ncor)],
              stdout=PIPE, stderr=STDOUT).communicate()
    groups = []
    with open("tmp/temp.mclout") as fin:
        for line in fin:
            groups.append(line[:-1].split())
    # system("rm tmp/*.bst tmp/temp.mclin tmp/temp.mclout") # Need to change this but later
    return groups


def paired_list(lst):
    if len(lst) < 2:
        return []
    return [[lst[i], lst[i+1]] for i in range(len(lst)-1)]


def mcl_rerunner(ifl, ncor, df, seq_ids, filename, distant):
    df[df[0].isin(seq_ids) & df[1].isin(seq_ids)].to_csv(filename, index=False, header=False, sep="\t")
    groups = mclf(ifl, ncor, filename, distant)
    if len(groups) == 1:
        return groups
    tgroup = []
    for group in groups:
        # print(group)
        group_seq_sizes = [int(x.split("___")[-1]) for x in group]
        # if np.sqrt(np.mean(group_seq_sizes)) > np.std(group_seq_sizes):
        if np.sqrt(np.min(group_seq_sizes)) > np.std(group_seq_sizes):
                tgroup += [group]
        else:
            tgroup += mcl_rerunner(ifl, ncor, df, group, filename, distant)
    return tgroup


def reanalysis_blast_small(ifl, evalue,
               minseq, minmap, keepclean, minlen, mindiff, algo, identity,
               infile):
    mclfile = blastpf(algo, identity, evalue, keepclean, mindiff, minmap,
                True, (infile, infile))
    if stat(mclfile).st_size:
        seq_sim_df = pd.read_table(mclfile, header=None)
        seq_ids = list(set(seq_sim_df[0]) | set(seq_sim_df[1]))
        return mcl_rerunner(ifl, 1, seq_sim_df, seq_ids, mclfile, True)
    else:
        return []



def reanalysis_blat(ifl, minseq, minmap, keepclean, minlen, mindiff, algo, identity, ncor, infile):
    """Fragments larger cluster in smaller."""
    blat(algo, "pblat", ncor, keepclean, False, identity, minlen,
         mindiff, minmap, True, (infile, infile))
    mclfile = "tmp/%s.mclin" % path.split(infile)[1].split(".faa")[0]
    seq_sim_df = pd.read_table(mclfile, header=None)
    seq_ids = list(set(seq_sim_df[0]) | set(seq_sim_df[1]))
    return mcl_rerunner(ifl, ncor, seq_sim_df, seq_ids, mclfile, False)


def return_sequences(faaf, seqdf, col):
    """Return sequences for each column."""
    seq_dct = SeqIO.to_dict(SeqIO.parse("%s/%s.faa" % (faaf, col), "fasta"))
    toreturn = {}
    for sqids in seqdf[col]:
        if sqids == "*":
            continue
        for seqid in sqids.split(","):
            seqinf = seqid.split(":")
            toreturn["%s___%s___%s" %
                     (col, seqinf[0], seqinf[1])] = seq_dct[seqinf[0]].seq
    return toreturn


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--faaf", help="folder containing protein sequence"
              " file from differnt sample", type=str,
              default="aa2",
              show_default=True)
@click.option("--identity_close", help="Expected minimum sequence similarity"
              " most distant samples for closely related sample. Will be used"
              "when --distant option is False",
              type=float, default=0.8, show_default=True)
@click.option("--identity_distant", help="Expected minimum sequence similarity"
              " most distant samples for distantly related sample. Will be"
              " used when --distant option is True",
              type=float, default=0.25, show_default=True)
@click.option("--ncor", help="number of cores", type=int,
              default=22, show_default=True)
@click.option("--outfile", help="output cluster file", type=str,
              default='clusters.clstr', show_default=True)
@click.option("--pblatpath", help="BLAT path", type=str,
              default="pblat", show_default=True)
@click.option("--makeblastdb", help="makeblastdb path. Functional"
              " When distant option is used", type=str, default="makeblastdb",
              show_default=True)
@click.option("--blastp", help="Blastp path. Function when distant is option"
              " in use", type=str, default="blastp", show_default=True)
@click.option("--ifl", help="Inflation factor. Active if mcl algorithm will be"
              " used (between 0 and 4.0)", type=float, default=4.0,
              show_default=True)
@click.option("--distant", help="Samples from distatly replated organisms",
              type=bool, default=True, show_default=True)
@click.option("--evalue", help="evalue for blast search. Valid for distant"
              " samples only",
              type=float, default=1e-10, show_default=True)
@click.option("--minseq", help="Run MCL of sequence group more than n"
              " sequences", type=int, default=1, show_default=True)
@click.option("--algo", help="Choose similarity calculation algorithm"
              " It will be used in blast search",
              type=click.Choice(['blast', 'anm', 'min']),
              default='blast', show_default=True)
@click.option("--keepclean", help="Keep deleting intermediate files"
              "to suppress disk usage", type=bool, default=True,
              show_default=True)
@click.option("--seed", help="Random seed for pairing of files",
              type=int, default=1234, show_default=True)
@click.option("--minlen", help="Minimum Sequence Length",
              type=int, default=50, show_default=True)
@click.option("--mindiff", help="Sequence length Difference relative to\
                longest in pair",  type=float, default=0.5,
                show_default=True)  # None to be ignored
@click.option("--minmap", help="Minimum mapping relative to longest sequences",
              type=float, default=0.5,
              show_default=True)  # None to be ignored
def run(faaf, identity_close, identity_distant, ncor, outfile, pblatpath,
        evalue, distant, algo, blastp, ifl, makeblastdb, minseq, keepclean,
        seed, minlen, mindiff, minmap):
    """Cluster Generating function."""
    """The program is to genrate protein sequence cluster from sequences
    distributed in different sample files base on their sequence similarities
    reported by BLAT tool. Depending on different parameters different results
    will be provided.
    Each for contain sequence id in format of orgarmism/sample_sequenceid
    Not Trying to find sequence from different organisms which might have
    same fuction Trying to bring similar sequences together only"""
    system("rm clusters.clstr")
    if not path.isdir(faaf):
        click.echo("Folder path \"%s\" doesn't exist." % faaf)
        exit(0)

    if not is_tool(pblatpath):
        click.echo("pblat is not in the given path %s" % pblatpath)
        exit(0)
    np.random.seed(seed)

    if distant:
        if not is_tool(makeblastdb):
            click.echo("makeblastdb is not in the given path %s" % makeblastdb)
            exit(0)
        if not is_tool(blastp):
            click.echo("blastp is not in the given path %s" % blastp)
            exit(0)

    # if mcl:
    #     if not is_tool(mclpath):
    #         click.echo("mcl is not in the given path %s" % mcl)
    #         exit(0)

    if cpu_count() < ncor:
        click.echo("Number of core on the system is less than given in option.\
                   Therefore %d is being used in the analysis" % cpu_count())
        ncor = cpu_count()

    if path.isdir("tmp"):
        rmtree("tmp")
    makedirs("tmp")
    click.echo("Please have patience. It might take a while to finish ...")
    pool = Pool(ncor)
    aa_files = glob("%s/*" % faaf)
    aa_files_count = len(aa_files)
    click.echo("%d files found in input folder" % aa_files_count)  # move this to blat function
    connected_ids = nx.Graph()
    # aa_files_count = len(aa_files)
    if distant:
        identity = identity_distant
    else:
        identity = identity_close
    beginning = True
    click.echo("Running BLAT to indetify highly similar sequences .....")
    cycle = 2
    tmp_idtty = 1.
    for _ in range(cycle):
        tmp_idtty = (tmp_idtty + identity) / 2.
    # tmp_idtty = identity
    n = 0
    while aa_files_count > 1:
        print(n)
        n += 1

        aa_files.sort()
        # TODO: bring halft simililarity concept back

        file_pairs = randomfilepairs(aa_files, randompairs(aa_files_count))
        cor = ncor//aa_files_count + 1

        func = partial(blat, algo, pblatpath, cor, keepclean, beginning,
                       tmp_idtty, minlen, mindiff, minmap, False)
        file_lists_pair = [pool.apply_async(func, ([pair]))
                           for pair in file_pairs]
        beginning = False
        aa_files = []
        for file_n_lists in file_lists_pair:
            file_, lists = file_n_lists.get()
            aa_files.append(file_)
            connected_ids.add_edges_from(lists)
        aa_files_count = len(aa_files)
    sequences = {}
    files_n_seq = {}
    # print("Anmol")
    # exit(1)
    # print(seq_sizes)
    for rec in SeqIO.parse(aa_files[0], 'fasta'):
        sequences[rec.id] = rec.seq
    #
    if not distant:
        base_name = aa_files[0].split('.faa')[0]
        click.echo("Running BLAT to report distantly related sequences...")
        Popen([pblatpath, '-threads=%d' % ncor,
               '-prot', '-noHead', aa_files[0], aa_files[0],
               '%s.psl' % base_name],
              stdout=PIPE, stderr=STDOUT).communicate()
        mapped = pd.read_table(aa_files[0].split('.faa')[0]+'.psl',
                               header=None,
                               usecols=[0, 9, 10, 13, 14])
        if keepclean:
            system("rm %s %s" % (aa_files[0],
                                 '%s.psl' % base_name))
        mapped = mapped[(mapped[9] != mapped[13])]
        if algo == 'blast':  # TODO: move this code in blat function up
            mapped.loc[:, 'identity'] = (mapped[0] * 1. /
                                  mapped[[10, 14]].mean(axis=1))
        elif algo == 'min':
            mapped.loc[:, 'identity'] = (mapped[0] * 1. /
                                  mapped[[10, 14]].min(axis=1))
        elif algo == 'max':
            mapped.loc[:, 'identity'] = (mapped[0] * 1. /
                                  mapped[[10, 14]].max(axis=1))
        else:
            mapped.loc[:, 'q_r'] = mapped[10]-mapped[12]
            mapped.loc[:, 'd_r'] = mapped[14]-mapped[16]
            mapped.loc[:, 'identity'] = mapped[0]/(mapped[[0, 1, 5, 7]].sum(axis=1) +
                                            mapped[[11, 15]].max(axis=1) +
                                            mapped[['q_r', 'd_r']].max(axis=1)
                                            )

        mapped = mapped[mapped['identity'] >= identity]
        mapped[[9, 13, "identity"]].to_csv("%s.bst" % base_name,
                                           sep="\t", header=False,
                                           index=False)
        del mapped
    if distant:
        system("rm tmp/*")
        click.echo("Running BLAST to indetify distantly related sequences")
        seq_count = len(sequences)
        seq_ids = list(sequences.keys())
        sequences_perfile = seq_count//ncor + 1
        # print(identity, algo, evalue)
        for i, j in enumerate(range(0, seq_count, sequences_perfile)):
            with open("tmp/%d.faa" % i, "w") as fout:
                for k in seq_ids[j:j+sequences_perfile]:
                    fout.write(">%s\n%s\n" % (k, sequences[k]))
        # Test Part
        pool.map(makeblastdbf, glob("tmp/*.faa"))
        # Test Part End

        # Test Part
        func = partial(blastpf, algo, identity, evalue, keepclean,
                       mindiff, minmap, False)
        pool.map_async(func, its.product(glob("tmp/*.faa"), repeat=2), chunksize=1)
        pool.close()
        pool.join()
        pool = Pool(ncor)

        # for fl in glob("tmp/*.faa"):
        #     func = partial(blastpf, algo, identity, evalue, keepclean,
        #                    mindiff, minmap, False, fl)
        #     pool.map(func,  glob("tmp/*.faa"))

        # func = partial(blastpf, algo, identity, evalue, keepclean,
        #                mindiff, minmap, False)
        # def blastpf(algo, identity, evalue, keepclean, ncor, mindiff,
        # minmap, infile)
        # pool.map(func,  glob("tmp/*.faa"))
    for connected_list in pool.map(paired_list, mclf(ifl, ncor, None, distant), chunksize=1):
        connected_ids.add_edges_from(connected_list)

    nxacc = nx.algorithms.components.connected
    connected_ids = list(nxacc.connected_components(connected_ids))
    seq_sizes = {}

    base_names = [path.split(fl)[1].split(".")[0]
                  for fl in glob("%s/*" % faaf)]

    func = partial(id_arrange_df, base_names)
    dataframes = pd.concat(pool.map(func, connected_ids, chunksize=1), ignore_index=True)
    del connected_ids
    # Start here *****
    # # Identifying weird groups
    # groups2reanalyse = dataframes[dataframes["std"] >
    #                               np.sqrt(dataframes["mean"])]
    groups2reanalyse = dataframes[(dataframes["seq_count"] >
                                  dataframes["samp_count"]) |
                                  (dataframes["std"] >
                                  np.sqrt(dataframes["min"]))]
    col = ['cluster', 'samp_count', 'seq_count', 'min', 'median', 'mean',
           'std', 'max']
    ifl = ifl * 2
    groups2reanalyse = []
    if len(groups2reanalyse):
        print(len(groups2reanalyse), len(dataframes))
        print("Some sequences")
        time.sleep(5)
        # pass
        dataframes = dataframes.drop(groups2reanalyse.index)
        new_groups = []
        seq_files_col = groups2reanalyse.columns.difference(col)
        sequences = {}
        func = partial(return_sequences, faaf, groups2reanalyse[seq_files_col])
        for seqdict in map(func, seq_files_col):
            sequences.update(seqdict)

        if not distant:
            groups2reanalyse = groups2reanalyse[seq_files_col]
            groups2reanalyse = groups2reanalyse.values.tolist()
            system("rm tmp/*")
            for num, group_lst in enumerate(groups2reanalyse):
                with open("tmp/%d.faa" % num, "w") as fout:
                    for i, grp_gp in enumerate(group_lst):
                        if grp_gp == "*":
                            continue
                        sq_lst = grp_gp.split(",")
                        for sq_ls in sq_lst:
                            sq = sq_ls.split(":")
                            seqid = "%s___%s___%s" % (seq_files_col[i], sq[0], sq[1])
                            fout.write(">%s\n%s\n" % (seqid, sequences[seqid]))
                #

            # Find the resulst and merge
            func = partial(reanalysis_blat, ifl, minseq, minmap, keepclean, minlen, mindiff, algo, identity, 1)
            for lst in pool.map(func,  glob("tmp/*.faa")):
                # print(lst)
                new_groups += lst
            # new_groups += reanalysis_blat(ifl, minseq, minmap, keepclean, minlen, mindiff, algo, identity, "tmp/%d.faa" % num, ncor) # Make in multicore individual



        # Run multi processing system
        else:
            evalue = 10
            multi_groups2reanalyse = groups2reanalyse[groups2reanalyse["seq_count"] >= ncor**2 ]
            multi_groups2reanalyse = multi_groups2reanalyse[seq_files_col]
            single_groups2reanalyse = groups2reanalyse[groups2reanalyse["seq_count"] < ncor**2 ]
            single_groups2reanalyse = single_groups2reanalyse[seq_files_col]
            groups2reanalyse = single_groups2reanalyse[seq_files_col]
            groups2reanalyse = groups2reanalyse.values.tolist()
            system("rm tmp/*")
            for num, group_lst in enumerate(groups2reanalyse):
                with open("tmp/%d.faa" % num, "w") as fout:
                    for i, grp_gp in enumerate(group_lst):
                        if grp_gp == "*":
                            continue
                        sq_lst = grp_gp.split(",")
                        for sq_ls in sq_lst:
                            sq = sq_ls.split(":")
                            seqid = "%s___%s___%s" % (seq_files_col[i], sq[0], sq[1])
                            # print(seqid)
                            fout.write(">%s\n%s\n" % (seqid, sequences[seqid]))
            #(ifl, evalue,
                          # minseq, minmap, keepclean, minlen, mindiff, algo, identity,
                           #infile)
            func = partial(reanalysis_blast_small, ifl, evalue,
                           minseq, minmap, keepclean, minlen, mindiff, algo,
                           identity)
            for lst in pool.map(func,  glob("tmp/*.faa"), chunksize=1):
                print(lst)
                new_groups += lst
            # Some other function
            # Work on removing upper part
            # groups2reanalyse = multi_groups2reanalyse[seq_files_col]
            # groups2reanalyse = groups2reanalyse.values.tolist()
            # for num, group_lst in groups2reanalyse:
            #     seq_ids = []
            #     for i, grp_gp in enumerate(group_lst):
            #         if grp_gp == "*":
            #             continue
            #         sq_lst = grp_gp.split(",")
            #         for sq_ls in sq_lst:
            #             sq = sq_ls.split(":")
            #             seqid = "%s___%s___%s" % (seq_files_col[i], sq[0], sq[1])
            #             # print(seqid)
            #             seq_ids.append(seqid)
            #     chunk_size = len(seq_ids) // ncor + 1
            #     system("tmp/*") # Clear file
            #     for num in range(ncor):
            #         with open("tmp/%d.faa" % num, "w") as fout:
            #             for sq in seq_ids[num * chunk_size: (num + 1) * chunk_size]:
            #                 fout.write(">%s\n%s\n" % (sq, sequences))
            #     func = partial(blastpf, algo, identity, evalue, keepclean,
            #                    1, mindiff, minmap, False)
            #     pool.map(func,  glob("tmp/*.faa"))
            #     # Generate a common file
            #     df = pd.read_table("file", header=None, sep="\t")
            #     new_groups += mcl_rerunner(ifl, ncor, df, seq_ids, filename)  # find a finame for consitancy
                # Send the information to mcl calculator
                # Make it run on multicore sequences are fragmented here
                # Run and
                #


        del sequences
        print(new_groups, base_names)
        func = partial(id_arrange_df, base_names)
        new_groups = pd.concat(pool.map(func, new_groups, chunksize=1), ignore_index=True)
        dataframes = pd.concat([dataframes, new_groups])
        del new_groups
    print(len(dataframes))
    time.sleep(5)
    dataframes = dataframes.sort_values(by=["samp_count", "seq_count"],
                                        ascending=[False, True])
    dataframes.loc[:, "cluster"] = ["cluster_%d" % d for d in
                             range(1, dataframes.shape[0]+1)]

    # Rearanging the columns
    col += list(set(dataframes.columns)-set(col))
    min_seq = 2
    min_samp = 2
    if min_seq:
        dataframes = dataframes[dataframes["seq_count"] >= min_seq]

    if min_samp:
        dataframes = dataframes[dataframes["samp_count"] >= min_samp]

    dataframes[col].to_csv(outfile, sep="\t", index=False)
    click.echo("Result file %s generated." % outfile)
    click.echo("Finished.....")
    # rmtree("tmp")


if __name__ == '__main__':
    run()
