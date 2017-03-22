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


def blat(blatpath, cor, keepclean, beginning, identity, pair):
    """Pairwise blat comparision to search for the best sequences."""
    if beginning:
        for i, fl in enumerate(pair):
            flname = path.split(fl)[1]
            base_name = flname.split(".")[0]
            pair[i] = "tmp/%s" % flname
            with open(pair[i], "w") as fout:
                for rec in SeqIO.parse(fl, 'fasta'):
                    fout.write(">%s___%s___%d\n%s\n" % (base_name, rec.id,
                                                        len(rec.seq), rec.seq))

    if len(pair) == 1:
        return pair[0], []

    sequences = {}
    for p in pair:
        for rec in SeqIO.parse(p, 'fasta'):
            sequences[rec.id] = rec.seq
    base_name = pair[0].split('.')[0]
    grbg = Popen([blatpath, '-threads=%d' % cor,
                  '-prot', '-noHead', pair[0], pair[1],
                  '%s.psl' % base_name],
                 stdout=PIPE, stderr=STDOUT).communicate()
    # exit(1)

    try:
        mapped = pd.read_table('%s.psl' % base_name, header=None,
                               usecols=[0, 9, 10, 13, 14])
        # TODO:Fix empty file problem
        # Clsuter and choose
        mapped = mapped[mapped[0] >
                        (mapped[[10, 14]].min(axis=1) *
                         identity)][[9, 13]].values.tolist()
        connected_ids = nx.Graph()
        connected_ids.add_edges_from(mapped)
        nxacc = nx.algorithms.components.connected
        connected_ids = list(nxacc.connected_components(connected_ids))
        selected_ref = []
        for connected_id in  connected_ids:
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

        # mapped = mapped.sort_values(by=[0], ascending=[False])
        # reported = set()
        # unreported_indexes = []
        # for i, row in mapped.iterrows():
        #     if row[9] in reported or row[13] in reported:
        #         continue
        #     unreported_indexes.append(i)
        #     reported |= set(row[[9, 13]])
        # mapped = mapped.ix[unreported_indexes]
        # id_pairs = mapped[[9, 13]].values.tolist()
        # long_seq_ids = mapped[mapped[10] > mapped[14]][9].tolist()
        # long_seq_ids += mapped[mapped[10] <= mapped[14]][13].tolist()
        # ids_to_file = set(long_seq_ids) | (set(sequences.keys()) -
        #                                    set(flat_list(id_pairs)))

        with open(pair[0], "w") as fout:
            for uid in ids_to_file:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        if keepclean:
            system("rm %s %s.psl" % (pair[1], base_name))
        return pair[0], mapped #id_pairs
    except EOFError:  # Empty file
        with open(pair[0], "w") as fout:
            for uid in sequences:
                fout.write(">%s\n%s\n" % (uid, sequences[uid]))
        if keepclean:
            system("rm %s %s.psl" % (pair[1], base_name))
        return pair[0], []


def makeblastdbf(infile):
    """Created blast search database for blastp."""

    Popen(["makeblastdb", "-in", infile, "-dbtype", "prot",
           "-out", "%s.bdb" % infile.split('.')[0]],
          stdout=PIPE, stderr=STDOUT).communicate()
    return


def blastpf(algo, identity, evalue, keepclean, infile):
    """Running BLAST and seleting best searches."""
    makeblastdbf(infile)
    # TODO: Fix things here
    infile_ = infile.split(".")[0]
    # ncor = 1
    for query in glob("tmp/*.faa"): # TODO: Need to flip the situation
        system("blastp -db %s.bdb -query %s -evalue %e -outfmt 6 >> %s.bst"
               % (infile_, query, evalue,  infile_))
    data = pd.read_table("%s.bst" % infile_, header=None,
                         usecols=range(10))
    data = data[data[0] != data[1]]
    data["qsize"] = data[0].map(lambda x: int(x.split("___")[2]))
    data["ssize"] = data[1].map(lambda x: int(x.split("___")[2]))

    if algo == "blast":
        data = data[((data[3]*data[2] / 100.) /
                     data[["qsize", "ssize"]].mean(axis=1)) > identity][[0, 1]]
    elif algo == "min":
        data = data[((data[3]*data[2] * 2 / 100.) /
                     data[["qsize", "ssize"]].min(axis=1)) > identity][[0, 1]]

    elif algo == "max":
        data = data[((data[3]*data[2] * 2 / 100.) /
                     data[["qsize", "ssize"]].max(axis=1)) > identity][[0, 1]]
    else:
        data['q_r'] = data["qsize"] - data[7]
        data['d_r'] = data["ssize"] - data[9]
        data = data[((data[3]*data[2] * 2 / 100.)/(
            data[3] + data[['q_r', 'd_r']].max(axis=1) +
            data[[6, 8]].max(axis=1))) > identity][[0, 1]]

    return data.values


def blastpf2(algo, identity, evalue, keepclean, infile):
    """Running BLAST and seleting best searches."""
    makeblastdbf(infile)
    # TODO: Fix things here
    infile_ = infile.split(".")[0]
    # ncor = 1
    for query in glob("tmp/*.faa"): # TODO: Need to flip the situation
        system("blastp -db %s.bdb -query %s -evalue %e -outfmt 6 >> %s.bst"
               % (infile_, query, evalue,  infile_))
    data = pd.read_table("%s.bst" % infile_, header=None,
                         usecols=range(10))
    data = data[data[0] != data[1]]
    data["qsize"] = data[0].map(lambda x: int(x.split("___")[2]))
    data["ssize"] = data[1].map(lambda x: int(x.split("___")[2]))

    if algo == "blast":
        data["iden"] = ((data[3]*data[2] / 100.) /
                             data[["qsize", "ssize"]].mean(axis=1))
    elif algo == "min":
        data["iden"] = ((data[3]*data[2] * 2 / 100.) /
                        data[["qsize", "ssize"]].min(axis=1))

    elif algo == "max":
        data["iden"] = ((data[3]*data[2] * 2 / 100.) /
                        data[["qsize", "ssize"]].max(axis=1))
    else:
        data['q_r'] = data["qsize"] - data[7]
        data['d_r'] = data["ssize"] - data[9]
        data["iden"] = ((data[3]*data[2] * 2 / 100.)/(
            data[3] + data[['q_r', 'd_r']].max(axis=1) +
            data[[6, 8]].max(axis=1)))
    data[[0, 1, "iden"]].to_csv("%s.mclin" % infile_,
                                header=False, index=False)

    return












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


def id_arrange_df(sample_ids, ids): # sizes,
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


def mclf2(inflation):
    system("cat tmp/*.mclin > tmp/temp.mclinall")
    Popen(["mcl", "tmp/temp.mclinall", "--abc", "-I", str(inflation),
           "-o", "tmp/temp.mcloutall"],
          stdout=PIPE, stderr=STDOUT).communicate()
    groups = []
    with open("tmp/temp.mcloutall") as fin:
        for line in fin:
            groups.append(line[:-1].split())
    return groups


def mclf(prog, algo, mcl, dbcreater, mapper, inflation,
         identity, evalue, tseq, minseq, keepclean, connected_id):
    connected_id = list(connected_id)
    if len(connected_id) < minseq:
        for id_ in connected_id:
            del tseq[id_]
        return [connected_id]
    infile = "tmp/%s.faa" % connected_id[0]
    with open(infile, "w")  as fout:
        for id_ in connected_id:
            fout.write(">%s\n%s\n" % (id_, tseq[id_]))
            del tseq[id_]

    if prog == 'blast':
        # TODO: Add file deletion option
        # infile_ = infile.split(".")[0]
        grbg = Popen([dbcreater, "-in", infile, "-dbtype", "prot",
                     "-out", "%sdb" % infile],
                     stdout=PIPE, stderr=STDOUT).communicate()

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
        # elif algo == "max":
        #     system("%s -db %sdb -query %s -outfmt '6 pident qseqid qlen"
        #            " sseqid slen length' -evalue %e |"
        #            "awk '{if((2*$1*$6/(100.*($3+$5)) >= %f) && ($2 != $4))"
        #            " print $2\"\t\"$4\"\t\"$3>$5?$3:$5}' > %s.mclin"
        #            % (mapper, infile, infile, evalue, identity, infile))
        else:
            system("%s -db %sdb -query %s -outfmt '6 pident qseqid qstart"
                   " qend qlen sseqid sstart send slen length' -evalue %e"
                   " > %s.blst"
                   % (mapper, infile, infile, evalue, infile))

            data = pd.read_table("%s.blst" % infile, header=None)
            data = data[data[1] != data[5]]
            data['q_r'] = data[4] - data[3]
            data['d_r'] = data[8] - data[7]
            data['iden'] = (data[0] * data[9] / 100.) / (
                data[9] + data[['q_r', 'd_r']].max(axis=1) +
                data[[2, 6]].max(axis=1))
            data = data[data['iden'] >= identity]
            if not len(data):
                return []
            data[[1, 5, 'iden']].to_csv("%s.mclin" % infile, header=False,
                                        index=False, sep="\t")

    else:
        grbg = Popen([mapper, "-prot", "-noHead", infile, infile,
                      "%s.psl" % infile],
                     stdout=PIPE, stderr=STDOUT).communicate()
        # system("%s -prot -noHead %s %s %s.psl" % (mapper, infile,
        #                                           infile, infile))
        data = pd.read_table("%s.psl" % infile, header=None)
        if keepclean:
            system("rm %s.psl" % infile)
        data = data[data[9] != data[13]]
        # if '101008_00560' in data[9] or '101008_00560' in data[13]:
        #     print(data)
        if algo == 'blast':
            data["iden"] = data[0] * 1. / data[[10, 14]].mean(axis=1)
        elif algo == 'min':
            data["iden"] = data[0] * 1. / data[[10, 14]].min(axis=1)
        # elif algo == 'max':
        #     data["iden"] = data[0] * 1. / data[[10, 14]].max(axis=1)
        else:
            data['q_r'] = data[10]-data[12]
            data['d_r'] = data[14]-data[16]
            data['iden'] = data[0]/(data[[0, 1, 5, 7]].sum(axis=1) +
                                    data[[11, 15]].max(axis=1) +
                                    data[['q_r', 'd_r']].max(axis=1)
                                    )

        data = data[data['iden'] >= identity]
        if not len(data):
            if keepclean:
                system("rm %s*" % infile)
            return []
        data[[9, 13, 'iden']].to_csv("%s.mclin" % infile, header=False,
                                     index=False, sep="\t")

    Popen([mcl, "%s.mclin" % infile, "--abc", "-I", str(inflation),
           "-o", "%s.mclout" % infile],
          stdout=PIPE, stderr=STDOUT).communicate()
    groups = []
    with open("%s.mclout" % infile) as fin:
        for line in fin:
            groups.append(line[:-1].split())
    if keepclean:
        system("rm %s*" % infile)
    return groups


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--faaf", help="folder containing protein sequence"
              " file from differnt sample", type=str,
              default="/home/devil/Documents/Tools/Clustering/aa2",
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
              default=1, show_default=True)
@click.option("--outfile", help="output cluster file", type=str,
              default='clusters.clstr', show_default=True)
@click.option("--pblatpath", help="BLAT path", type=str,
              default="pblat", show_default=True)
@click.option("--makeblastdb", help="makeblastdb path. Functional"
              " When distant option is used", type=str, default="makeblastdb",
              show_default=True)
@click.option("--blastp", help="Blastp path. Function when distant is option"
              " in use", type=str, default="blastp", show_default=True)
@click.option("--mcl", help="Run MCL algorithm", type=bool,
              default=True, show_default=True)
@click.option("--mclpath", help="MCL path", type=str,
              default="mcl", show_default=True)
@click.option("--ifl", help="Inflation factor. Active if mcl algorithm will be"
              " used (between 0 and 4.0)", type=float, default=4.0,
              show_default=True)
@click.option("--distant", help="Samples from distatly replated organisms",
              type=bool, default=False, show_default=True)
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
def run(faaf, identity_close, identity_distant, ncor, outfile, pblatpath,
        evalue, distant, mcl, algo, blastp, ifl, makeblastdb,
        mclpath, minseq, keepclean, seed):
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
    aa_files = glob("%s/*" % faaf)
    aa_files_count = len(aa_files)
    click.echo("%d files found in input folder" % aa_files_count)

    # for i, source in enumerate(glob("%s/*" % faaf)):
    #     # TODO: Merging this code with above one to reduce the use of
    #     # intermediate files
    #     base_name = source.split("/")[-1].split(".")[0]
    #     base_names.append(base_name)
    #     file_index[base_name] = i
    #     with open("%s/tmp/%d.faa" % (getcwd(), i), "w") as fout:
    #         for rec in SeqIO.parse(source, 'fasta'):
    #             fout.write(">%s___%s\n%s\n" % (base_name, rec.id, rec.seq))
    #             seq_sizes["%s___%s" % (base_name, rec.id)] = len(rec.seq)
    # click.echo("%d sequence in %d files....." % (len(seq_sizes),
    #                                              len(base_names)))

    # aa_files = glob("tmp/*")
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

    while aa_files_count > 1:

        aa_files.sort()
        # TODO: bring halft simililarity concept back
        # TODO: Choose similarities based on smallest sequences

        file_pairs = randomfilepairs(aa_files, randompairs(aa_files_count))
        cor = ncor//aa_files_count + 1

        func = partial(blat, pblatpath, cor, keepclean, beginning, tmp_idtty)
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
    # exit(1)
    # print(seq_sizes)
    for rec in SeqIO.parse(aa_files[0], 'fasta'):
        sequences[rec.id] = rec.seq
    #
    if not distant:
        base_name = aa_files[0].split('.faa')[0]
        click.echo("Running BLAT to report distantly related sequences...")
        grbg = Popen([pblatpath, '-threads=%d' % ncor,
                      '-prot', '-noHead', aa_files[0], aa_files[0],
                      '%s.psl' % base_name],
                     stdout=PIPE, stderr=STDOUT).communicate()

        # mapped = pd.read_table(aa_files[0].split('.faa')[0]+'.psl',
        #                        header=None,
        #                        usecols=[0, 1, 5, 7, 9, 10, 11, 12,
        #                                 13, 14, 15, 16])
        mapped = pd.read_table(aa_files[0].split('.faa')[0]+'.psl',
                               header=None,
                               usecols=[0, 9, 10, 13, 14])
        # TODO:Fix empty file problem
        # Clsuter and choose
        if keepclean:
            system("rm %s %s" % (aa_files[0],
                                 '%s.psl' % base_name))
        mapped = mapped[(mapped[9] != mapped[13])]
        mapped = mapped[mapped[0] >
                        (mapped[[10, 14]].min(axis=1) *
                         identity)][[9, 13]].values.tolist()
        # if algo == 'blast':
        #     mapped = mapped[(mapped[0] * 1. / mapped[[10, 14]].mean(axis=1)) >=
        #                     identity]
        # elif algo == 'min':
        #     mapped = mapped[(mapped[0] * 1. / mapped[[10, 14]].min(axis=1)) >=
        #                     identity]
        # elif algo == 'max':
        #     mapped = mapped[(mapped[0] * 1. / mapped[[10, 14]].max(axis=1)) >=
        #                     identity]
        # else:
        #     mapped['q_r'] = mapped[10]-mapped[12]
        #     mapped['d_r'] = mapped[14]-mapped[16]
        #     mapped['iden'] = mapped[0]/(mapped[[0, 1, 5, 7]].sum(axis=1) +
        #                                 mapped[[11, 15]].max(axis=1) +
        #                                 mapped[['q_r', 'd_r']].max(axis=1)
        #                                 )
        #     mapped = mapped[mapped['iden'] >= identity]

        # id_pairs = mapped[[9, 13]].values
        connected_ids.add_edges_from(mapped)
        # connected_ids.add_edges_from(id_pairs)
        del mapped, grbg
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

        # for uid in sequences:
        #     file_id = uid.split('___')[0]
        #     if file_id not in files_n_seq:
        #         files_n_seq[file_id] = [uid]
        #     else:
        #         files_n_seq[file_id].append(uid)
        # for k in files_n_seq:
        #     with open("tmp/%s.faa" % k, "w") as fout:
        #         for sq in files_n_seq[k]:
        #             fout.write(">%s\n%s\n" % (sq, sequences[sq]))
        # func = partial(makeblastdbf, makeblastdb)
        # # TODO: Merge above function with blastf to delete the dbs
        # pool.map(func, glob("tmp/*.faa"))
        # del grbg
        # print("Anmol")
        # need this area to be rearranged
        func = partial(blastpf,"min", identity, evalue, keepclean)
        # Ask for seq sizes while copying the files
        # use for loop with core
        # for fl in glob("tmp/*.faa"):
        #     connected_ids.add_edges_from(func(fl))
        # for fl in glob("tmp/*.faa"):
        #     connected_ids.add_edges_from(func(fl))
        pairs = pool.map(func,  glob("tmp/*.faa"))
        for pair in pairs:
            connected_ids.add_edges_from(pair)
    #
    nxacc = nx.algorithms.components.connected
    connected_ids = list(nxacc.connected_components(connected_ids))
    # # Replace graph withn gereator not with the list to reduce memory usage
    #
    # del graph
    #
    # exit(1)
    # base_names = [] # Fix this issue
    seq_sizes = {}
    if mcl:
        try:
            system("rm tmp/*")
        except:
            pass
        # Consider that 20 core will be used
        click.echo("Running MCL over grouped sequences")
        tseq = {}
        for fl in glob("%s/*" % faaf):
            base_name = path.split(fl)[1].split('.')[0]
            for rec in SeqIO.parse(fl, "fasta"):
                tseq["%s___%s___%d" %
                     (base_name, rec.id, len(rec.seq))] = str(rec.seq)
                seq_sizes["%s___%s___%d" %
                          (base_name, rec.id, len(rec.seq))] = len(rec.seq)

        if distant:
            func = partial(mclf, "blast", algo, mclpath, makeblastdb, blastp,
                           ifl, identity, evalue, tseq, minseq, keepclean)
            low_seq = []
            high_seq = []
            for connected_id in connected_ids:
                if len(connected_id) <= 400:
                    low_seq.append(connected_id)
                else:
                    high_seq.append(connected_id)
            connected_ids = []
            clusters = pool.map(func, low_seq)
            for cluster in clusters:
                connected_ids += cluster
            for connected_id in high_seq:
                try:
                    system("rm tmp/*")
                except:
                    pass
                seq_count = len(connected_id)
                seq_ids = list(connected_id)
                sequences_perfile = seq_count//ncor + 1
                # print(identity, algo, evalue)
                for i, j in enumerate(range(0, seq_count, sequences_perfile)):
                    with open("tmp/%d.faa" % i, "w") as fout:
                        for k in seq_ids[j:j+sequences_perfile]:
                            fout.write(">%s\n%s\n" % (k, tseq[k]))
                func = partial(blastpf2, algo, identity, evalue, keepclean)
                pool.map(func,  glob("tmp/*.faa"))
                connected_ids += mclf2(ifl)


        else:
            func = partial(mclf, "blat", algo, mclpath, None, pblatpath, ifl,
                           identity, None, tseq, minseq, keepclean)

            clusters = pool.map(func, connected_ids)
            del connected_ids, tseq
            connected_ids = []
            for cluster in clusters:
                connected_ids += cluster
    # else:
    #     for fl in glob("%s/*" % faaf):
    #         base_name = path.split(fl)[1].split(".")[0]
    #         base_names.append(base_name)
    #         for rec in SeqIO.parse(fl, "fasta"):
    #             seq_sizes["%s___%s" % (base_name, rec.id)] = len(rec.seq)
    # print(seq_sizes, base_names)

    base_names = [path.split(fl)[1].split(".")[0]
                  for fl in glob("%s/*" % faaf)]

    func = partial(id_arrange_df, base_names)  # , seq_sizes
    # dataframest = [pool.apply_async(func, ([fl])) for fl in connected_ids]
    dataframes = pd.concat(pool.map(func, connected_ids), ignore_index=True)
    del connected_ids

    # dataframes = dataframest[0].get()
    # for dt in dataframest[1:]:
    #     dataframes = pd.concat([dataframes, dt.get()], ignore_index=True)
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
    click.echo("Finished.....")
    rmtree("tmp")


if __name__ == '__main__':
    run()
