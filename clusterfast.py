#!/usr/bin/env python
"""Nothing To Write here."""

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
from os import path, makedirs, system, stat
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
import tempfile
import itertools as its
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
    """BLAST dataframe cleaner."""
    mapped = mapped2.copy()
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


def blatf(algo, pblat, cor, beginning, identity, minlen,
         mindiff, minmap, mclinfile, tmpd, pair):
    """Pairwise blat comparision to search for the best sequences."""
    if beginning:
        for i, fl in enumerate(pair):
            flname = path.split(fl)[1]
            base_name = flname.split(".")[0]
            pair[i] = "%s/%s" % (tmpd, flname)
            with open(pair[i], "w") as fout:
                for rec in SeqIO.parse(fl, 'fasta'):
                    if len(rec.seq) < minlen or len(rec.seq) > 20000:
                        continue
                    fout.write(">%s___%s___%d\n%s\n" % (base_name, rec.id,
                                                        len(rec.seq), rec.seq))
    if len(pair) == 1:
        return pair[0], []

    sequences = {}
    for p in pair:
        for rec in SeqIO.parse(p, 'fasta'):
            sequences[rec.id] = rec.seq
    base_name = pair[0].split('.')[0]
    if mclinfile != "orth":
        outfile = '%s.bst' % base_name
    else:
        f1 = path.split(pair[1])[1].split(".faa")[0]
        outfile = '%s_%s.bst' % (base_name, f1)
    Popen([pblat, '-threads=%d' % cor, "-out=blast8",
           '-prot', '-noHead', pair[0], pair[1],
           outfile],
          stdout=PIPE, stderr=STDOUT).communicate()
    if stat(outfile).st_size:
        mapped = pd.read_table(outfile, header=None)
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped = mapped[mapped["identity"] >= identity]

        if mclinfile == "orth":
            mapped = mapped[mapped["db"] != mapped["qr"]]
            return mapped[["db", "qr", "bits"]]

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
        # if keepclean:
        #     system("rm %s %s.bst" % (pair[1], base_name))
        return pair[0], mapped  # id_pairs
    else:  # Empty file
        with open(pair[0], "w") as fout:
            for uid in sequences:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        return pair[0], []


def blast_table_analysis(blastdata, adaptive):
    """This function is based on ProteinOrtho4.0 Algorithm."""
    blast_data = blastdata.copy()
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

            for i, row in qr_db_blast_data_max.iterrows():
                selected_indexes += list(
                    qr_db_blast_data[
                        (qr_db_blast_data['qr_samp'] == row['qr_samp']) &
                        (qr_db_blast_data['bits'] >= adaptive * row['bits'])
                        ].index)
    blast_data = blast_data.loc[list(set(selected_indexes)), ]

    blast_data = blast_data.drop(["db_samp", "qr_samp"], axis=1)
    blast_data.loc[blast_data["db"] > blast_data["qr"], ["db", "qr"]
                   ] = blast_data.loc[blast_data["db"] > blast_data["qr"],
                                      ["qr", "db"]].values
    blast_data_group = blast_data.groupby(["db", "qr"]).size().reset_index()
    blast_data_group = blast_data_group[blast_data_group[0] > 1]
    blast_data_group = blast_data_group[["db", "qr"]].values.tolist()
    return blast_data_group

def norm_mult(max_degree, graph, connected_ids_list, node_count):
    to_return = np.ones(node_count)
    for i in range(node_count):
        to_return[i] = ((2 * max_degree) -
                        len(graph.neighbors(
                            connected_ids_list[i])))

    return to_return


def neigh_index(df, graph, connected_ids_list, node_count):
    to_return = {}
    for i in range(node_count):
        to_return[i] = df['ids'].isin(graph.neighbors(connected_ids_list[i]))
    return to_return



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
        connected_ids_list.sort()
        connected_ids_df = pd.DataFrame.from_dict({'ids': connected_ids_list})
        node_count = len(connected_ids_list)
        sub_graph = graph.subgraph(connected_ids_list)
        edge_count = len(sub_graph.edges())
        node_df = sub_graph.degree(connected_ids_list)
        max_degree = np.max(list(node_df.values()))
        max_pair = node_count * (node_count - 1)/2
        if max_pair == edge_count:
            to_return.append(connected_ids_list)
        else:
            normmult = norm_mult(max_degree, sub_graph, connected_ids_list,
                                 node_count)
            neighindex = neigh_index(connected_ids_df, sub_graph,
                                     connected_ids_list, node_count)
            x = np.random.random(node_count)
            x_hat = x - np.mean(x)
            last_len = np.linalg.norm(x_hat)
            if last_len == 0:
                last_len = 1e-9
            norm = x_hat/last_len
            while True:
                x = np.zeros(node_count)
                for i in range(node_count):
                    x[i] += np.sum(norm[neighindex[i]])
                norm *= normmult
                norm += x
                x_hat = norm - np.mean(norm)
                current_len = np.linalg.norm(x_hat)
                if current_len == 0:
                    current_len = 1e-9
                norm = x_hat/current_len
                if np.abs(current_len - last_len) < 0.001:
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


def makeblastdbf(makeblastdb, infile):
    """Created blast search database for blastp."""

    Popen([makeblastdb, "-in", infile, "-dbtype", "prot",
           "-out", "%s.bdb" % infile.split('.')[0]],
          stdout=PIPE, stderr=STDOUT).communicate()
    return


def blastpf(blastp, algo, identity, evalue, mindiff, minmap,
            mclinfile, db_query):
    """Running BLAST and seleting best searches."""
    db, infile = db_query

    db_ = db.split(".faa")[0]
    query_id = path.split(infile)[1].split(".faa")[0]
    outfilebs = "%s_%s" % (db_, query_id)
    if mclinfile == "orth":
        system("%s -db %s.bdb -query %s -evalue %e"
               " -outfmt 6 > %s.bst" % (blastp, db_, infile,
                                        evalue, outfilebs))
        if stat("%s.bst" % outfilebs).st_size:
            mapped = pd.read_table("%s.bst" % outfilebs, header=None)
            mapped = blast_dataframe(mapped, mindiff, minmap, algo)
            mapped = mapped[mapped["identity"] >= identity]
            mapped = mapped.sort_values(['eval', 'bits', "identity"],
                                        ascending=[True, False, False])
            mapped = mapped.drop_duplicates(['db', 'qr'])
            return mapped[["db", "qr", "bits"]]
        else:
            return pd.DataFrame()


    system("%s -db %s.bdb -query %s -evalue %e"
           " -outfmt 6 > %s.bst" % (blastp, db_, infile,
                                    evalue, outfilebs))
    if stat("%s.bst" % outfilebs).st_size:
        mapped = pd.read_table("%s.bst" % outfilebs, header=None)
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped = mapped[mapped["identity"] >= identity]
        mapped = mapped.sort_values(["eval", "bits", "identity"],
                                    ascending=[True, False, False])
        mapped = mapped.drop_duplicates(["db", "qr"])
        return mapped[["db", "qr"]].values.tolist()
    return []


def randompairs(n_count):
    """Generate random pairs for the given list of the files."""
    n_list = list(range(n_count))
    np.random.shuffle(n_list)
    for i in range(0, n_count, 2):
        yield n_list[i:i+2]


def randomfilepairs(file_list, pairs):
    """Generating random pair of files."""
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
    """Generate list of paires."""
    if len(lst) < 2:
        return []
    return [[lst[i], lst[i+1]] for i in range(len(lst)-1)]


def multi_sep(adaptive, conn_thresh, df):
    if len(df):
        df = blast_table_analysis(df, adaptive)
        connected_ids = nx.Graph()
        connected_ids.add_edges_from(df)
        connected_ids = group_seprator(connected_ids, conn_thresh)
        return connected_ids
    else:
        return []

def blast_big(blastp, makeblastdb, tmpd, evalue, minmap,
              mindiff, algo, identity, ncor, adaptive, conn_thresh):
    """Fragment large cluster in smaller."""
    files = glob("%s/*.faa" % tmpd)
    func = partial(makeblastdbf, makeblastdb)
    pool = Pool(ncor)
    pool.map(func, files)
    func = partial(blastpf, blastp, algo, identity, evalue, mindiff,
                   minmap, "orth")
    df = pool.map(func, its.product(files, repeat=2), chunksize=1)
    df = pd.concat(df, ignore_index=True)
    return df


def blast_small(blastp, makeblastdb, evalue, minmap,
                mindiff, algo, identity, adaptive, conn_thresh, infile):
    """Fragment large cluster in smaller."""
    makeblastdbf(makeblastdb, infile)
    df = blastpf(blastp, algo, identity, evalue, mindiff, minmap,
                 "orth", (infile, infile))
    return df


def reanalysis_blat(pblat, minmap, mindiff,
                    algo, identity, ncor, adaptive, conn_thresh, tmpd, infile):
    """Fragment large cluster in smaller."""
    df = blatf(algo, pblat, ncor, False, identity, None,
               mindiff, minmap, "orth", tmpd, (infile, infile))
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
@click.option("-faaf", help="folder containing .faa protein sequence files"
              " file from differnt sample", type=str,
              default="aa2",
              show_default=True)
@click.option("-identity", help="Expected minimum sequence similarity"
              " For distant 0.25 will be used"
              " For close 0.8 will be used as default",
              type=float, default=None, show_default=True)
@click.option("-ncor", help="number of cores", type=int,
              default=22, show_default=True)
@click.option("-outfile", help="output cluster file", type=str,
              default='clusters.clstr', show_default=True)
@click.option("-pblat", help="PBLAT path", type=str,
              default="pblat", show_default=True)
@click.option("-makeblastdb", help="makeblastdb path. Functional"
              " When distant option is used", type=str, default="makeblastdb",
              show_default=True)
@click.option("-blastp", help="Blastp path. Function when distant is option"
              " in use", type=str, default="blastp", show_default=True)
@click.option("-distant", help="Samples from distatly replated organisms"
              " (for close pblat and for distant)",
              type=bool, default=True, show_default=True)
@click.option("-evalue", help="evalue for blast search. Valid for distant"
              " samples only",
              type=float, default=1e-10, show_default=True)
@click.option("-algo", help="Choose similarity calculation algorithm"
              " It will be used in blast search",
              type=click.Choice(['blast', 'anm']),
              default='anm', show_default=True)
@click.option("-seed", help="Random seed for pairing of files",
              type=int, default=1234, show_default=True)
@click.option("-minlen", help="Minimum Sequence Length",
              type=int, default=50, show_default=True)
@click.option("-mindiff", help="Sequence length Difference relative to"
              " longest in pair", type=float, default=0.5,
              show_default=True)  # None to be ignored
@click.option("-minmap", help="Minimum mapping relative to longest sequences",
              type=float, default=0.5,
              show_default=True)  # None to be ignored
@click.option("-conn_thresh", help="Connection Threshold", type=float,
              default=0.1, show_default=True)
@click.option("-adaptive", help="Adaptive search value", type=float,
              default=0.95, show_default=True)
def run(faaf, identity, ncor, outfile, pblat,
        evalue, distant, algo, blastp, makeblastdb,
        seed, minlen, mindiff, minmap, adaptive, conn_thresh):
    """Generate quick cluster based on given sequence samples."""
    """The program uses ProteinOrth4 algorithm to genrate protein sequence
    cluster from sequences distributed in different sample files base on their
    sequence similarities reported by BLAT tool. Depending on different
    parameters different results will be provided.
    Each for contain sequence id in format of orgarmism/sample_sequenceid
    Not Trying to find sequence from different organisms which might have
    same fuction Trying to bring similar sequences together only"""
    np.random.seed(seed)
    # Checking ........
    if path.isfile(outfile):
        click.echo("Given output file already exists. Exiting ....")
        exit(0)

    if not path.isdir(faaf):
        click.echo("Folder path \"%s\" doesn't exist." % faaf)
        exit(0)

    if not is_tool(pblat):
        click.echo("pblat is not in the given path %s" % pblat)
        exit(0)

    if distant:
        if not is_tool(makeblastdb):
            click.echo("makeblastdb is not in the given path %s" % makeblastdb)
            exit(0)
        if not is_tool(blastp):
            click.echo("blastp is not in the given path %s" % blastp)
            exit(0)
    if not identity:
        identity = 0.25 if distant else 0.8

    if cpu_count() < ncor:
        click.echo("Number of core on the system is less than given in option.\
                   Therefore %d is being used in the analysis" % cpu_count())
        ncor = cpu_count()

    # Check done ...
    tmpd = tempfile.mkdtemp(dir=os.getcwd())

    pool = Pool(ncor)
    aa_files = glob("%s/*.faa" % faaf)
    aa_files_count = len(aa_files)
    click.echo("%d files found in input folder" % aa_files_count)
    connected_ids = nx.Graph()
    beginning = True
    click.echo("Please have patience. It might take a while to finish ...")
    click.echo("Running BLAT to indetify highly similar sequences .....")
    tmp_idtty = (1. + identity) / 2.
    # Looking for highly similar sequences .....
    while aa_files_count > 1:
        aa_files.sort()
        cor = ncor//aa_files_count + 1
        func = partial(blatf, algo, pblat, cor, beginning,
                       tmp_idtty, minlen, mindiff, minmap, "X", tmpd)
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
    for rec in SeqIO.parse(aa_files[0], 'fasta'):
        sequences[rec.id] = rec.seq

    # Looking for remotely related sequences ....
    if not distant:
        system("rm %s/*" % tmpd)
        click.echo("Running BLAST to indetify distantly related sequences")
        seq_count = len(sequences)
        seq_ids = list(sequences.keys())
        sequences_perfile = seq_count//ncor + 1
        for i, j in enumerate(range(0, seq_count, sequences_perfile)):
            with open("%s/%d.faa" % (tmpd, i), "w") as fout:
                for k in seq_ids[j:j+sequences_perfile]:
                    fout.write(">%s\n%s\n" % (k, sequences[k]))
        files = glob("%s/*.faa" % tmpd)
        func = partial(blatf, algo, pblat, 1, False, identity, None,
                   mindiff, minmap, "orth", tmpd)  # , (infile, infile)
        pairs = pool.map(func, its.product(files, repeat=2), chunksize=1)
        for pair in pairs:
            connected_ids.add_edges_from(pair[["db", "qr"]].values.tolist())

    if distant:
        system("rm %s/*" % tmpd)
        click.echo("Running BLAST to indetify distantly related sequences")
        seq_count = len(sequences)
        seq_ids = list(sequences.keys())
        sequences_perfile = seq_count//ncor + 1
        for i, j in enumerate(range(0, seq_count, sequences_perfile)):
            with open("%s/%d.faa" % (tmpd, i), "w") as fout:
                for k in seq_ids[j:j+sequences_perfile]:
                    fout.write(">%s\n%s\n" % (k, sequences[k]))
        files = glob("%s/*.faa" % tmpd)
        func = partial(makeblastdbf, makeblastdb)
        pool.map(func, files)
        # TODO: make it bidirectional to avoid below replications
        func = partial(blastpf, blastp, algo, identity, evalue,
                       mindiff, minmap, "X")
        pairs = pool.map(func, its.product(files, repeat=2), chunksize=1)
        # TODO: Use dataframes here.
        for pair in pairs:
            connected_ids.add_edges_from(pair)

    nxacc = nx.algorithms.components.connected
    connected_ids = list(nxacc.connected_components(connected_ids))

    base_names = [path.split(fl)[1].split(".")[0]
                  for fl in glob("%s/*" % faaf)]
    func = partial(id_arrange_df, base_names)
    dataframes = pd.concat(pool.map(func, connected_ids, chunksize=1),
                           ignore_index=True)
    del connected_ids
    full = True
    if full:
        groups2reanalyse = dataframes[(
            (dataframes["seq_count"] > dataframes["samp_count"]) |
            (dataframes["std"] > np.sqrt(dataframes["min"])) |
            (dataframes["min"] < mindiff * dataframes["max"])) &
                                      (dataframes["seq_count"] > 2)]
    else:
        # This part may be used for some other purpose in future
        groups2reanalyse = dataframes[dataframes["seq_count"] > 2]

    col = ['cluster', 'samp_count', 'seq_count', 'min', 'median', 'mean',
           'std', 'max']
        # Separating bigger and diverse clusters .....
    evalue **= 2 # Highly similar sequences expected in the
    if len(groups2reanalyse):
        click.echo("Separating Paralogs .....")
        click.echo(str(len(groups2reanalyse)) + " Clusters are being tried to"
                   " be fragmented again...... ")
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
            system("rm %s/*" % tmpd)
            for num, group_lst in enumerate(groups2reanalyse):
                with open("%s/%d.faa" % (tmpd, num), "w") as fout:
                    for i, grp_gp in enumerate(group_lst):
                        if grp_gp == "*":
                            continue
                        sq_lst = grp_gp.split(",")
                        for sq_ls in sq_lst:
                            sq = sq_ls.split(":")
                            seqid = "%s___%s___%s" % (seq_files_col[i], sq[0],
                                                      sq[1])
                            fout.write(">%s\n%s\n" % (seqid, sequences[seqid]))

            func = partial(reanalysis_blat, pblat, minmap, mindiff,
                           algo, identity, 1, adaptive, conn_thresh, tmpd)
            for lst in pool.map(func,  glob("%s/*.faa" % tmpd), chunksize=1):
                # print(len(lst))
                new_groups += lst

        # Run multi processing system
        else:   # chosen dues to small number sequences
            multi_groups2reanalyse = groups2reanalyse[
                groups2reanalyse["seq_count"] >= ncor**2
                ]
            multi_groups2reanalyse = multi_groups2reanalyse[seq_files_col]

            # Smaller groups for single cores
            single_groups2reanalyse = groups2reanalyse[
                groups2reanalyse["seq_count"] < ncor**2
                ]
            single_groups2reanalyse = single_groups2reanalyse[seq_files_col]
            dfs = []

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
                sequences_perfile = seq_count//ncor + 1
                system("rm %s/*" % tmpd)
                seq_ids = list(temp_seq)
                for i, j in enumerate(range(0, seq_count, sequences_perfile)):
                    with open("%s/%d.faa" % (tmpd, i), "w") as fout:
                        for k in seq_ids[j:j+sequences_perfile]:
                            fout.write(">%s\n%s\n" % (k, sequences[k]))
                # print("Kiran")

                # Simply write a function which could do it for you.
                # def blast_big(blastp, makeblastdb, tmpd, evalue, minmap,
                #               mindiff, algo, identity, ncor, adaptive, conn_thresh):
                # ncor = 1
                dfs.append(blast_big(blastp, makeblastdb, tmpd, evalue,
                                        minmap, mindiff, algo, identity, ncor,
                                        adaptive, conn_thresh))
            groups2reanalyse = single_groups2reanalyse.values.tolist()
            system("rm %s/*" % tmpd)
            for num, group_lst in enumerate(groups2reanalyse):
                with open("%s/%d.faa" % (tmpd, num), "w") as fout:
                    for i, grp_gp in enumerate(group_lst):
                        if grp_gp == "*":
                            continue
                        sq_lst = grp_gp.split(",")
                        for sq_ls in sq_lst:
                            sq = sq_ls.split(":")
                            seqid = "%s___%s___%s" % (seq_files_col[i], sq[0],
                                                      sq[1])
                            fout.write(">%s\n%s\n" % (seqid, sequences[seqid]))
            # def blast_small(blastp, makeblastdb, evalue, minmap,
            #                 mindiff, algo, identity, adaptive, conn_thresh, infile):
            func = partial(blast_small, blastp, makeblastdb, evalue, minmap,
                           mindiff, algo, identity, adaptive, conn_thresh)

            #  For smaller group, send it directly to blast and get the result
            for lst in pool.map(func,  glob("%s/*.faa" % tmpd), chunksize=1):
                dfs.append(lst)
            func = partial(multi_sep, adaptive, conn_thresh)
            for lst in pool.map(func, dfs, chunksize=1):
                new_groups += lst


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
    #     pass
        # dataframes = dataframes[dataframes["samp_count"] >= min_samp]

    dataframes = dataframes.sort_values(by=["samp_count", "seq_count"],
                                        ascending=[False, True])
    dataframes.loc[:, "cluster"] = ["cluster_%d" % d for d in
                                    range(1, dataframes.shape[0]+1)]

    # Rearanging the columns
    col += list(set(dataframes.columns)-set(col))

    dataframes[col].to_csv(outfile, sep="\t", index=False)
    click.echo("Result file %s generated." % outfile)
    click.echo("Finished.....")
    rmtree(tmpd)


if __name__ == '__main__':
    run()
