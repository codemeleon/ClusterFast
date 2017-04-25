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
    base_name = pair[0].split('.')[0]
    Popen([blatpath, '-threads=%d' % cor,
           '-prot', '-noHead', pair[0], pair[1],
           '%s.psl' % base_name],
          stdout=PIPE, stderr=STDOUT).communicate()
    # exit(1)

    try:
        mapped = pd.read_table('%s.psl' % base_name, header=None,
                               usecols=[0, 9, 10, 11, 12,
                                        13, 14, 15, 16])
        # TODO:Fix empty file problem
        # Clsuter and choose
        # add a marker here to avoid researching and ask to generate only bst
        # files for the analysis
        if mindiff:
            mapped = mapped[(mapped[[10, 14]].min(axis=1)) >
                            mindiff * (mapped[[10, 14]].max(axis=1))]
        if minmap:
            mapped.loc[:, "minmap"] = mapped.apply(lambda x:
                                            max([x[12]-x[11],
                                                 x[16]-x[15]]) /
                                            max([x[10, 14]]))
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
            mapped.to_csv("tmp/mcl.mclin", header=False,
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
    except EOFError:  # Empty file
        with open(pair[0], "w") as fout:
            for uid in sequences:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        if keepclean:
            system("rm %s %s.psl" % (pair[1], base_name))
        return pair[0], []


def blast_dataframe(mapped, mindiff, minmap, algo):
    mapped = mapped[mapped[0] != mapped[1]]
    mapped.loc[:, "qsize"] = mapped[0].map(lambda x: int(x.split("___")[2]))
    mapped.loc[:, "ssize"] = mapped[1].map(lambda x: int(x.split("___")[2]))
    if minmap:
        mapped.loc[:, "minmap"] = mapped.apply(lambda x:
                                        max([x[7]-x[6],
                                             x[9]-x[8]]) /
                                        max([x['qsize', 'ssize']]))
        mapped = mapped[mapped["minmap"] > minmap]
        del mapped["minmap"]
    if mindiff:
        mapped = mapped[(mapped[['qsize', 'ssize']].min(axis=1)) >
                        (mapped[['qsize', 'ssize']].max(axis=1))]

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


def makeblastdbf(infile):
    """Created blast search database for blastp."""

    Popen(["makeblastdb", "-in", infile, "-dbtype", "prot",
           "-out", "%s.bdb" % infile.split('.')[0]],
          stdout=PIPE, stderr=STDOUT).communicate()
    return


def blastpf(algo, identity, evalue, keepclean, ncor, mindiff, minmap,
            mclinfile, infile):
    """Running BLAST and seleting best searches."""
    makeblastdbf(infile)
    infile_ = infile.split(".")[0]
    # TODO: Fix things here
    if mclinfile:
        system("blastp -num_threads %d -db %s.bdb -query %s -evalue %e"
               " -outfmt 6 > %s.bst" % (ncor, infile_, infile,
                                        evalue, infile_))
        mapped = pd.read_table("%s.bst" % infile_, header=None,
                               usecols=range(10))
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped.to_csv("tmp/mcl.mclin", header=False,
                      index=False, sep="\t")
        return
    # ncor = 1
    for query in glob("tmp/*.faa"):  # TODO: Need to flip the situation
        query_id = path.split(query)[1].split(".faa")[0]
        outfilebs = "%s_%s" % (infile_, query_id)
        print(outfilebs)
        system("blastp -num_threads %d -db %s.bdb -query %s -evalue %e"
               " -outfmt 6 > %s.bst" % (ncor, infile_, query,
                                        evalue, outfilebs))

        mapped = pd.read_table("%s.bst" % outfilebs, header=None,
                               usecols=range(10))
        mapped = blast_dataframe(mapped, mindiff, minmap, algo)
        mapped[[0, 1, 'identity']].to_csv("%s.bst" % outfilebs, header=None,
                                          index=None, sep="\t")

    return  # mapped.values



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


def mclf(inflation, ncor, infile):
    if infile:
        Popen(["mcl", infile, "--abc", "-I", str(inflation),
               "-o", "tmp/temp.mclout", "-te", str(ncor)],
              stdout=PIPE, stderr=STDOUT).communicate()
    else:
        system("cat tmp/*.bst > tmp/temp.mclin")
        Popen(["mcl", "tmp/temp.mclin", "--abc", "-I", str(inflation),
               "-o", "tmp/temp.mclout", "-te", str(ncor)],
              stdout=PIPE, stderr=STDOUT).communicate()
    groups = []
    with open("tmp/temp.mclout") as fin:
        for line in fin:
            groups.append(line[:-1].split())
    system("rm tmp/*.bst tmp/temp.mclin tmp/temp.mclout")
    return groups


def paired_list(lst):
    if len(lst) < 2:
        return []
    return [[lst[i], lst[i+1]] for i in range(len(lst)-1)]



def reanalysis(clusterframe, sequences, sample_ids, distant, ncor, ifl, evalue,
               minseq, minmap, keepclean, minlen, mindiff, algo, identity):
    """Fragments larger cluster in smaller."""
    final_list = []
    cls_lsts = []
    for group_lst in clusterframe:
        if not distant:
            infile = "tmp/frame_seq.fa"
            with open(infile, "w") as fout:
                for i, grp_gp in enumerate(group_lst):
                    if grp_gp == "*":
                        continue
                    sq_lst = grp_gp.split(",")
                    for sq_ls in sq_lst:
                        sq = sq_ls.split(":")
                        seqid = "%s___%s___%s" % (sample_ids[i], sq[0], sq[1])
                        # print(seqid)
                        fout.write(">%s\n%s\n" % (seqid, sequences[seqid]))

            blat(algo, "pblat", ncor, keepclean, False, identity, minlen,
                 mindiff, minmap, True, (infile, infile))
            cls_lsts = mclf(ifl, ncor, "tmp/mcl.mclin")
        else:
            ## ***** # Count sequences, split and combine
            system("rm tmp/*")
            ## *****
            blastpf(algo, identity, evalue, keepclean, ncor, mindiff, minmap,
                    True, infile)

        if len(cls_lsts) == 1:
            final_list += cls_lsts
        else:
            for lst in cls_lsts:
                # lst_seq_sizes = [int(x.split("___")[:-1]) for x in lst]
                col = ['cluster', 'samp_count', 'seq_count', 'min', 'median', 'mean',
                       'std', 'max']
                tlst = id_arrange_df(sample_ids, lst)
                if len(tlst[np.sqrt(tlst['mean']) > tlst['std']]):
                        final_list.append(lst)
                else:
                    tlst = tlst[tlst.columns.difference(col)].values.tolist()
                    newlist = reanalysis(tlst, sequences, sample_ids, distant,
                                         ncor, ifl, evalue, minseq, minmap,
                                         keepclean, minlen, mindiff, algo,
                                         identity)
                    final_list += newlist
    # print(final_list)
    return final_list


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
              default="aa",
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
              default=20, show_default=True)
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
                longest in pair",  type=float, default=None,
                show_default=True)  # None to be ignored
@click.option("--minmap", help="Minimum mapping relative to longest sequences",
              type=float, default=None,
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
    print("Anmol")
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
        func = partial(blastpf, algo, identity, evalue, keepclean,
                       1, mindiff, minmap, False)
        # def blastpf(algo, identity, evalue, keepclean, ncor, mindiff,
        # minmap, infile)
        pool.map(func,  glob("tmp/*.faa"))
    for connected_list in pool.map(paired_list, mclf(ifl, ncor, None)):
        connected_ids.add_edges_from(connected_list)

    nxacc = nx.algorithms.components.connected
    connected_ids = list(nxacc.connected_components(connected_ids))
    seq_sizes = {}

    base_names = [path.split(fl)[1].split(".")[0]
                  for fl in glob("%s/*" % faaf)]

    func = partial(id_arrange_df, base_names)
    dataframes = pd.concat(pool.map(func, connected_ids), ignore_index=True)
    del connected_ids
    # Start here *****
    # # Identifying weird groups
    groups2reanalyse = dataframes[dataframes["std"] >
                                  np.sqrt(dataframes["mean"])]
    col = ['cluster', 'samp_count', 'seq_count', 'min', 'median', 'mean',
           'std', 'max']
    if len(groups2reanalyse):
        # Remove sequences from original list
        dataframes = dataframes.drop(groups2reanalyse.index)
        seq_files_col = groups2reanalyse.columns.difference(col)
        func = partial(return_sequences, faaf,
                       groups2reanalyse[seq_files_col])
        seq_dicts = map(func, seq_files_col)
        fileseq = {}
        for dct in seq_dicts:
            fileseq.update(dct)
        del seq_dicts
        groups2reanalyse = groups2reanalyse[seq_files_col]
        groups2reanalyse = groups2reanalyse.values.tolist()
        newlist = reanalysis(groups2reanalyse, fileseq, seq_files_col,
                             distant, ncor, ifl, evalue, minseq, minmap, keepclean,
                             minlen, mindiff, algo, identity)
        func = partial(id_arrange_df, base_names)
        newlist = pd.concat(pool.map(func, newlist), ignore_index=True)
        dataframes = pd.concat(dataframes, newlist)
        del newlist
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
