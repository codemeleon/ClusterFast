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



def blast_dataframe(mapped2, mindiff, minmap, algo):
    mapped = mapped2.copy()
    # mapped = mapped[mapped[0] != mapped[1]]
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
        mapped = mapped.drop(['q_r', 'd_r'], axis=1)
        mapped = mapped.rename(columns={0: "db", 1: "qr",
                                        10: "eval", 11: "bits"})
    return mapped

## Write a blast file parser
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
    Popen([blatpath, '-threads=%d' % cor, "-out=blast8",
           '-prot', '-noHead', pair[0], pair[1],
           '%s.bst' % base_name],
          stdout=PIPE, stderr=STDOUT).communicate()
    # exit(1)
    if stat('%s.bst' % base_name).st_size:
        mapped = pd.read_table("%s.bst" % base_name, header=None)
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped = mapped[mapped["identity"] >= identity]

        if mclinfile == "orth":
            mapped = mapped[mapped["db"] != mapped["qr"]]
            return mapped


        mapped = mapped[["db", "qr"]].values.tolist()
        connected_ids = nx.Graph()
        connected_ids.add_edges_from(mapped)
        nxacc = nx.algorithms.components.connected
        connected_ids = list(nxacc.connected_components(connected_ids))
        selected_ref = []
        for connected_id in connected_ids:
            sz = 0
            same_len_seq = []
            for _id in connected_id:
                sq_sz = len(sequences[_id])
                if sq_sz > sz:
                    sz = sq_sz
                    same_len_seq = []
                if sq_sz == sz:
                    same_len_seq.append(_id)
            same_len_seq.sort()
            selected_ref.append(same_len_seq[0])
        ids_to_file = set(selected_ref) | (set(sequences.keys()) -
                                           set(flat_list(mapped)))
        with open(pair[0], "w") as fout:
            for uid in ids_to_file:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        if keepclean:
            system("rm %s %s.bst" % (pair[1], base_name))
        return pair[0], mapped  # id_pairs
    else:  # Empty file
        with open(pair[0], "w") as fout:
            for uid in sequences:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        if keepclean:
            system("rm %s %s.bst" % (pair[1], base_name))
        return pair[0], []

def blast_table_analysis(blastdata, adaptive):
    """This function is based on ProteinOrtho4.0 Algorithm."""
    blast_data = blastdata.copy()#[[0, 1, 10, 11]]
    # blast_data = blast_data.rename(columns={0:'db', 1:'qr', 10:'eval', 11:'bits'})
    blast_data = blast_data[blast_data['db'] != blast_data['qr']]
    blast_data["db_samp"] = [x.split("___")[0] for x in blast_data["db"]]
    blast_data["qr_samp"] = [x.split("___")[0] for x in blast_data["qr"]]
    selected_indexes = []
    for db_samp in set(blast_data["db_samp"]):
        db_blast_data = blast_data[blast_data["db_samp"] == db_samp]
        for db in set(db_blast_data["db"]):
            qr_db_blast_data = db_blast_data[db_blast_data["db"] == db]
            qr_db_blast_data_max = qr_db_blast_data.groupby(['qr_samp']
                                                            )['bits'].max(
                                                            ).reset_index()
            try:
                bits_min = min(qr_db_blast_data_max.loc[
                    qr_db_blast_data_max['qr_samp'] != db_samp, 'bits'])
            except ValueError:
                bits_min = 0

            drop_ix = qr_db_blast_data_max[
                (qr_db_blast_data_max['qr_samp'] == db_samp) &
                (qr_db_blast_data_max['bits'] < bits_min)
                ].index
            if len(drop_ix):
                qr_db_blast_data_max = qr_db_blast_data_max.drop(drop_ix)

            # bits_max = 0.9 * bits_max #- np.sqrt(bits_max) # Choose one of the idea
            # qr_db_blast_data_max = qr_db_blast_data_max.loc[
            #     qr_db_blast_data_max['bits'] >= bits_max,
            # ]
            for i, row in qr_db_blast_data_max.iterrows():
                selected_indexes += list(
                    qr_db_blast_data[
                        (qr_db_blast_data['qr_samp'] == row['qr_samp']) &
                        (qr_db_blast_data['bits'] >= adaptive * row['bits'])
                        ].index)
    blast_data = blast_data.loc[selected_indexes, ]

    blast_data = blast_data.drop(["db_samp", "qr_samp"], axis=1)
    blast_data.loc[blast_data["db"] > blast_data["qr"], ["db", "qr"]
                   ] = blast_data.loc[blast_data["db"] > blast_data["qr"],
                                      ["qr", "db"]].values
    # bits_mean = np.mean(blast_data["bits"])
    # bits_left = bits_mean - np.sqrt(bits_mean)
    # bits_right = bits_mean + np.sqrt(bits_mean)
    # blast_data = blast_data.loc[((blast_data["bits"] > bits_left) &
    #                              (blast_data["bits"] < bits_right)), ]
    blast_data_group = blast_data.groupby(["db", "qr"]).size().reset_index()
    blast_data_group = blast_data_group[blast_data_group[0] > 1]
    blast_data_group = blast_data_group[["db", "qr"]].values.tolist()
    return blast_data_group


def group_seprator(graph, conn_threshold):
    """This function is based on ProteinOrtho4.0 Algorithm."""
    to_return = []
    nxacc = nx.algorithms.components.connected
    connected_ids_lists = list(nxacc.connected_components(graph))
    for connected_ids_list in connected_ids_lists:
        if len(connected_ids_list) == 2:
            to_return.append(connected_ids_list)
    while True:
        node_df = graph.degree(graph.nodes())
        node_df = pd.DataFrame.from_dict({'node': list(node_df.keys()),
                                          'degree': list(node_df.values())})
        node_df = list(node_df.loc[node_df['degree'] < 2, "node"])
        if len(node_df):
            graph.remove_nodes_from(node_df)
        else:
            break
    connected_ids_lists = list(nxacc.connected_components(graph))
    for connected_ids in connected_ids_lists:
        connected_ids_list = list(connected_ids)
        node_count = len(connected_ids_list)
        sub_graph = graph.subgraph(connected_ids_list)
        edge_count = len(sub_graph.edges())
        node_df = sub_graph.degree(connected_ids_list)
        max_degree = np.max(list(node_df.values()))
        max_pair = node_count * (node_count - 1)/2
        if max_pair == edge_count:
            to_return.append(connected_ids_list)
        else:
            x = np.random.random(node_count)
            x_hat = x - np.mean(x)
            last_len = np.sqrt(np.sum(np.square(x_hat)))
            if last_len == 0:
                last_len = 1e-9
            norm = x_hat/last_len
            while True:
                x = np.zeros(node_count)
                for i in range(node_count):
                    for nb in sub_graph.neighbors(connected_ids_list[i]):
                        x[i] += norm[connected_ids_list.index(nb)]
                for k in range(node_count):
                    norm[k] *= ((2 * max_degree) -
                                len(sub_graph.neighbors(
                                    connected_ids_list[i])))
                norm += x
                x_hat = norm - np.mean(norm)
                current_len = np.sqrt(np.sum(np.square(x_hat)))
                if current_len == 0:
                    current_len = 1e-9
                norm = x_hat/current_len
                if abs(current_len - last_len) < 0.001:
                    break
                last_len = current_len

            alg_conn = (-1 * current_len + 2 * max_degree) / node_count
            if alg_conn < conn_threshold:
                new_groups = [np.array(connected_ids_list)[x_hat < 0],
                              np.array(connected_ids_list)[x_hat >= 0]]
                to_return += group_seprator(graph.subgraph(new_groups[0]),
                                            conn_threshold)
                to_return += group_seprator(graph.subgraph(new_groups[1]),
                                            conn_threshold)
            else:
                to_return.append(connected_ids_list)
    return to_return


def makeblastdbf(infile):
    """Created blast search database for blastp."""

    Popen(["makeblastdb", "-in", infile, "-dbtype", "prot",
           "-out", "%s.bdb" % infile.split('.')[0]],
          stdout=PIPE, stderr=STDOUT).communicate()
    return



def blastpf(algo, identity, evalue, keepclean, mindiff, minmap,
            mclinfile, db_query):
    """Running BLAST and seleting best searches."""
    db, infile = db_query
    infile_ = infile.split(".")[0]
    if mclinfile=="orth":
        # makeblastdbf(db)
        infile_bs = path.split(infile)[1].split(".fa")[0]
        system("blastp -db %s.bdb -query %s -evalue %e"
               " -outfmt 6 > %s.bst" % (infile_, infile,
                                        evalue, infile_))
        if stat("%s.bst" % infile_).st_size:
            mapped = pd.read_table("%s.bst" % infile_, header=None)
            mapped = blast_dataframe(mapped, mindiff, minmap, algo)
            mapped = mapped[mapped["identity"] >= identity]
            mapped = mapped.sort_values(['eval', 'bits', "identity"],
                                        ascending=[True, False, False])
            mapped = mapped.drop_duplicates(['db', 'qr'])
            return mapped
        else:
            return pd.DataFrame()

    db_ = db.split(".")[0]
    query_id = path.split(infile)[1].split(".faa")[0]
    outfilebs = "%s_%s" % (db_, query_id)
    system("blastp -db %s.bdb -query %s -evalue %e"
           " -outfmt 6 > %s.bst" % (db_, infile,
                                    evalue, outfilebs))
    if stat("%s.bst" % outfilebs).st_size:
        mapped = pd.read_table("%s.bst" % outfilebs, header=None)
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped = mapped[mapped["identity"] >= identity]
        mapped = mapped.sort_values(["eval", "bits", "identity"],
                                    ascending=[True, False, False])
        mapped = mapped.drop_duplicates(["db", "qr"])
        return mapped[["db","qr"]].values.tolist()

    return []



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


def id_arrange_df(sample_ids, ids):
    """Arranging clustered sequences ids in dataframe."""

    to_return = {}
    seq_size = []
    seq_count = len(ids)
    samp_count = 0
    for id_ in ids:
        samp, seq, sz = id_.split('___')
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



def paired_list(lst):
    if len(lst) < 2:
        return []
    return [[lst[i], lst[i+1]] for i in range(len(lst)-1)]


def blast_big(evalue, minmap, keepclean,
              mindiff, algo, identity, ncor, adaptive, conn_thresh):
    files = glob("tmp/*.faa")
    func = partial(makeblastdbf)
    pool = Pool(ncor)
    pool.map(func, files)
    func = partial(blastpf, algo, identity, evalue, keepclean, mindiff,
                   minmap, "orth")
    df = pool.map(func, its.product(files, repeat=2), chunksize=1)
    df = pd.concat(df)
    if len(df):
        df = blast_table_analysis(df, adaptive)
        connected_ids = nx.Graph()
        connected_ids.add_edges_from(df)
        connected_ids = group_seprator(connected_ids, conn_thresh)
        return connected_ids
    else:
        return []

def blast_small(evalue, minmap, keepclean,
                mindiff, algo, identity, adaptive, conn_thresh, infile):
    makeblastdbf(infile)
    df = blastpf(algo, identity, evalue, keepclean, mindiff, minmap,
                 "orth", (infile, infile))
    splt = "orth"
    if len(df):
        df = blast_table_analysis(df, adaptive)
        connected_ids = nx.Graph()
        connected_ids.add_edges_from(df)
        connected_ids = group_seprator(connected_ids, conn_thresh)
        return connected_ids
    else:
        return []


def reanalysis_blat(minmap, keepclean, mindiff,
                    algo, identity, ncor, adaptive, conn_thresh, infile):
    """Fragments larger cluster in smaller."""
    df = blat(algo, "pblat", ncor, keepclean, False, identity, None,
              mindiff, minmap, "orth", (infile, infile))
    if len(df):
        df = blast_table_analysis(df, adaptive)
        connected_ids = nx.Graph()
        connected_ids.add_edges_from(df)
        connected_ids = group_seprator(connected_ids, conn_thresh)
        return connected_ids
    else:
        return []

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
@click.option("--pblatpath", help="PBLAT path", type=str,
              default="pblat", show_default=True)
@click.option("--makeblastdb", help="makeblastdb path. Functional"
              " When distant option is used", type=str, default="makeblastdb",
              show_default=True)
@click.option("--blastp", help="Blastp path. Function when distant is option"
              " in use", type=str, default="blastp", show_default=True)
@click.option("--distant", help="Samples from distatly replated organisms",
              type=bool, default=True, show_default=True)
@click.option("--evalue", help="evalue for blast search. Valid for distant"
              " samples only",
              type=float, default=1e-10, show_default=True)
@click.option("--algo", help="Choose similarity calculation algorithm"
              " It will be used in blast search",
              type=click.Choice(['blast', 'anm']),
              default='anm', show_default=True)
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
@click.option("--conn_thresh", help="Connection Threshold", type=float,
              default=0.1, show_default=True)
@click.option("--adaptive", help="Adaptive search value", type=float,
              default=0.95, show_default=True)
def run(faaf, identity_close, identity_distant, ncor, outfile, pblatpath,
        evalue, distant, algo, blastp, makeblastdb, keepclean,
        seed, minlen, mindiff, minmap, adaptive, conn_thresh):
    """ClusterFast, to generate quick cluster based on given samples sequences."""
    """The program uses ProteinOrth4 algorithm to genrate protein sequence
    cluster from sequences distributed in different sample files base on their
    sequence similarities reported by BLAT tool. Depending on different
    parameters different results will be provided.
    Each for contain sequence id in format of orgarmism/sample_sequenceid
    Not Trying to find sequence from different organisms which might have
    same fuction Trying to bring similar sequences together only"""
    system("rm clusters.clstr")
    np.random.seed(seed)
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

    if cpu_count() < ncor:
        click.echo("Number of core on the system is less than given in option.\
                   Therefore %d is being used in the analysis" % cpu_count())
        ncor = cpu_count()

    if path.isdir("tmp"):
        rmtree("tmp")
    makedirs("tmp")
    click.echo("Please have patience. It might take a while to finish ...")
    pool = Pool(ncor)
    aa_files = glob("%s/*.faa" % faaf)
    aa_files_count = len(aa_files)
    click.echo("%d files found in input folder" % aa_files_count)
    connected_ids = nx.Graph()
    if distant:
        identity = identity_distant
    else:
        identity = identity_close
    beginning = True
    click.echo("Running BLAT to indetify highly similar sequences .....")
    tmp_idtty = (tmp_idtty + identity) / 2.
    while aa_files_count > 1:
        aa_files.sort()
        cor = ncor//aa_files_count + 1
        func = partial(blat, algo, pblatpath, cor, keepclean, beginning,
                       tmp_idtty, minlen, mindiff, minmap, "X")
        file_pairs = randomfilepairs(aa_files, randompairs(aa_files_count))
        file_lists_pair = pool.map(func, file_pairs)
        beginning = False
        aa_files = []
        for file_n_lists in file_lists_pair:
            file_, lists = file_n_lists
            aa_files.append(file_)
            connected_ids.add_edges_from(lists)
        aa_files_count = len(aa_files)
    sequences = {}
    files_n_seq = {}
    for rec in SeqIO.parse(aa_files[0], 'fasta'):
        sequences[rec.id] = rec.seq
    # ide
    if not distant:
        base_name = aa_files[0].split('.faa')[0]
        click.echo("Running BLAT to report distantly related sequences...")
        grbg = Popen([pblatpath, '-threads=%d' % ncor, "-out=blast8",
                      '-prot', '-noHead', aa_files[0], aa_files[0],
                      '%s.bst' % base_name],
                     stdout=PIPE, stderr=STDOUT).communicate()
        if stat('%s.bst' % base_name).st_size:
            mapped = pd.read_table("%s.bst" % base_name, header=None)
            mapped = blast_dataframe(mapped, mindiff, minmap, algo)
            mapped = mapped[mapped["identity"] >= 0.9 * identity]
            mapped = mapped[mapped["db"] != mapped["qr"]]
            connected_ids.add_edges_from(mapped[["db", "qr"]].values.tolist())
            del mapped
        del grbg
    if distant:
        system("rm tmp/*")
        click.echo("Running BLAST to indetify distantly related sequences")
        seq_count = len(sequences)
        seq_ids = list(sequences.keys())
        sequences_perfile = seq_count//ncor + 1
        for i, j in enumerate(range(0, seq_count, sequences_perfile)):
            with open("tmp/%d.faa" % i, "w") as fout:
                for k in seq_ids[j:j+sequences_perfile]:
                    fout.write(">%s\n%s\n" % (k, sequences[k]))
        files = glob("tmp/*.faa")
        func = partial(makeblastdbf)
        pool.map(func, files)
        func = partial(blastpf, algo, 0.9*identity, evalue, keepclean, mindiff,
                       minmap, "X")
        pairs = pool.map(func, its.product(files, repeat=2), chunksize=1)
        for pair in pairs:
            connected_ids.add_edges_from(pair)

    nxacc = nx.algorithms.components.connected
    connected_ids = list(nxacc.connected_components(connected_ids))
    seq_sizes = {}

    base_names = [path.split(fl)[1].split(".")[0]
                  for fl in glob("%s/*" % faaf)]
    func = partial(id_arrange_df, base_names)
    dataframes = pd.concat(pool.map(func, connected_ids, chunksize=1),
                           ignore_index=True)
    del connected_ids
    groups2reanalyse = dataframes[(
        (dataframes["seq_count"] > dataframes["samp_count"]) |
        (dataframes["std"] > np.sqrt(dataframes["min"])) |
        (dataframes["min"] < mindiff * dataframes["max"])) &
                                  (dataframes["seq_count"] > 2)]
    col = ['cluster', 'samp_count', 'seq_count', 'min', 'median', 'mean',
           'std', 'max']
    # groups2reanalyse = []
    click.echo("Separating Paralogs .....")
    if len(groups2reanalyse):
        print(len(groups2reanalyse), len(dataframes))
        print("Some sequences")
        time.sleep(5)
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
                            seqid = "%s___%s___%s" % (seq_files_col[i], sq[0],
                                                      sq[1])
                            fout.write(">%s\n%s\n" % (seqid, sequences[seqid]))
            func = partial(reanalysis_blat, minmap, keepclean,
                           mindiff, algo, identity, 1, adaptive, conn_thresh)
            for lst in pool.map(func,  glob("tmp/*.faa"), chunksize=1):
                print(len(lst))
                new_groups += lst



        # Run multi processing system
        else:
            evalue = 5
            multi_groups2reanalyse = groups2reanalyse[
                groups2reanalyse["seq_count"] >= ncor**2
                ]
            multi_groups2reanalyse = multi_groups2reanalyse[seq_files_col]

            # Smaller groups for single cores
            single_groups2reanalyse = groups2reanalyse[
                groups2reanalyse["seq_count"] < ncor**2
                ]
            single_groups2reanalyse = single_groups2reanalyse[seq_files_col]
            groups2reanalyse = single_groups2reanalyse.values.tolist()
            system("rm tmp/*")
            for num, group_lst in enumerate(groups2reanalyse):
                with open("tmp/%d.faa" % num, "w") as fout:
                    for i, grp_gp in enumerate(group_lst):
                        if grp_gp == "*":
                            continue
                        sq_lst = grp_gp.split(",")
                        for sq_ls in sq_lst:
                            sq = sq_ls.split(":")
                            seqid = "%s___%s___%s" % (seq_files_col[i], sq[0],
                                                      sq[1])
                            fout.write(">%s\n%s\n" % (seqid, sequences[seqid]))
            func = partial(blast_small, evalue,
                           minmap, keepclean, mindiff, algo,
                           identity, adaptive, conn_thresh)

            #  For smaller group, send it directly to blast and get the result
            for lst in pool.map(func,  glob("tmp/*.faa"), chunksize=1):
                new_groups += lst

            groups2reanalyse = multi_groups2reanalyse.values.tolist()
            for num, group_lst in enumerate(groups2reanalyse):
                temp_seq = {}
                for i, grp_gp in enumerate(group_lst):
                    if grp_gp == "*":
                        continue
                    sq_lst = grp_gp.split(",")
                    for sq_ls in sq_lst:
                        sq = sq_ls.split(":")
                        seqid = "%s___%s___%s" % (seq_files_col[i], sq[0],
                                                  sq[1])
                        temp_seq[seqid] = sequences[seqid]
                seq_count = len(temp_seq)
                equences_perfile = seq_count//ncor + 1
                system("rm tmp/*")
                for i, j in enumerate(range(0, seq_count, sequences_perfile)):
                    with open("tmp/%d.faa" % i, "w") as fout:
                        for k in seq_ids[j:j+sequences_perfile]:
                            fout.write(">%s\n%s\n" % (k, sequences[k]))

                # Simply write a function which could do it for you.
                new_groups += big_blast(evalue, minmap, keepclean, mindiff,
                                        algo, identity, ncor, adaptive,
                                        conn_thresh)




        del sequences
        func = partial(id_arrange_df, base_names)
        new_groups = pd.concat(pool.map(func, new_groups, chunksize=1),
                               ignore_index=True)
        dataframes = pd.concat([dataframes, new_groups])
        del new_groups

    min_seq = 2
    min_samp = 2
    if min_seq:
        dataframes = dataframes[dataframes["seq_count"] >= min_seq]

    # if min_samp:
    #     dataframes = dataframes[dataframes["samp_count"] >= min_samp]

    dataframes = dataframes.sort_values(by=["samp_count", "seq_count"],
                                        ascending=[False, True])
    dataframes.loc[:, "cluster"] = ["cluster_%d" % d for d in
                             range(1, dataframes.shape[0]+1)]

    # Rearanging the columns
    col += list(set(dataframes.columns)-set(col))

    dataframes[col].to_csv(outfile, sep="\t", index=False)
    click.echo("Result file %s generated." % outfile)
    click.echo("Finished.....")
    # rmtree("tmp")


if __name__ == '__main__':
    run()
