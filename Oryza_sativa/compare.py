#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time       : 2022/9/27
# @Author     : zzhen
# @File       : compare.py.py
# @Software   : PyCharm
# @Description: A script for calculating the identity and coverage between the recovered sequences and the gold standard sequences.
# @Copyright  : Copyright (c) 2022 by sculab, All Rights Reserved.

import os
import argparse
import re
import csv
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# get the abs_path of file or folder
def absolute_path(path: str) -> str:
    if not os.path.isabs(path):
        return os.path.abspath(path)
    return path


# get seq from files -> [{id1:seq1},{id2:seq2}]
def get_seq(fasta_file: str, max_seq_number: int = 100, seq_count_limit: int = False):
    seq_list = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_list.append({record.description: str(record.seq)})
    return seq_list


def get_alignment_score(query_id: str, query_seq: str, gold_seq: str):
    # score matrix
	# Global alignment with free end gaps   5/-4/-12/-3 (65%)   1/-1/-3/-2 (75%)
    match_score = 5
    mismatch_score = -4.0
    gap_open_score = -12
    gap_extend_score = -3
	# Global alignment, the same residue will be given 1 point, and different and gap will not deduct points
    alignments = pairwise2.align.globalms(query_seq, gold_seq, match_score, mismatch_score, gap_open_score,
                                          gap_extend_score)  
    # return score
    score = alignments[0][2]
    return {score: (query_id, query_seq)}


# alignment
def do_alignment(combined_file, aligned_file):
    cmd = "muscle3 -in {} -out {} -quiet".format(combined_file, aligned_file)
    subprocess.call(cmd, shell=True)
    return combined_file


def get_best_query_seq(query_path: str, gold_path: str, paired_dir:str):
    query_seq_list = get_seq(query_path, max_seq_number=100, seq_count_limit=False)
    gold_seq_list = get_seq(gold_path, max_seq_number=100, seq_count_limit=False)
    file_name = os.path.splitext(os.path.basename(query_path))[0]
    if len(query_seq_list) > 1 and len(gold_seq_list) > 1:
        print("Both query and ref contain multiple sequences, which are not supported in the current version")
        sys.exit()
    elif len(query_seq_list) < len(gold_seq_list):
        print("The number of sequences in the query is smaller than the number of sequences in the gold file.")
        print("Maybe the input directory is wrong")
        sys.exit()
    else:
        pass
    gold_seq = list(gold_seq_list[0].values())[0]
    # 1 vs 1
    if len(query_seq_list) == 1 and len(gold_seq_list) == 1:
        path = os.path.join(paired_dir, file_name + ".fasta")
        with open(query_path, "r") as f1, open(path, "w") as f:
            f.writelines(f1.readlines())
    elif len(query_seq_list) > 1 and len(gold_seq_list) == 1:
        score_info = dict()
        for i in range(len(query_seq_list)):
            score_info.update(
                get_alignment_score(list(query_seq_list[i].keys())[0], list(query_seq_list[i].values())[0], gold_seq))
        max_score = max(score_info.keys())
        path = os.path.join(paired_dir, file_name + ".fasta")
        with open(path, "w") as f:
            f.writelines([">" + score_info[max_score][0] + "\n", str(score_info[max_score][1]) + "\n"])
    else:
        pass


def get_best_query_seq_parallel(query_dir: str, gold_dir: str, pair_dir: str, thread_number: int):
    extension = (".fasta", ".fas", ".fa", ".fna", ".ffn", ".frn", ".faa", ".fna")
    query_file_list = [i for i in os.listdir(query_dir) if
                       os.path.isfile(os.path.join(query_dir, i)) and i.endswith(extension)]
    gold_file_list = [i for i in os.listdir(gold_dir) if
                      os.path.isfile(os.path.join(gold_dir, i)) and i.endswith(extension)]

    query_dict, gold_dict = dict(), dict()
    for i in query_file_list:
        query_dict.update({os.path.splitext(i)[0]: os.path.join(query_dir, i)})
    for i in gold_file_list:
        gold_dict.update({os.path.splitext(i)[0]: os.path.join(gold_dir, i)})
    shared_gene = [i for i in query_dict.keys() if i in gold_dict.keys()]
    combine_files = []  # [(query_path, gold_path, gene_name)]
    for i in shared_gene:
        combine_files.append((query_dict[i], gold_dict[i], i))
    results = []

    with ThreadPoolExecutor(max_workers=thread_number) as executor:
        futures = [executor.submit(get_best_query_seq, query_path=file[0], gold_path=file[1], paired_dir=pair_dir)
                   for file in combine_files]
        for future in as_completed(futures):
            results.append(future.result())


def do_combine(query_file: str, gold_file: str, combined_dir: str):
    extension = (".fasta", ".fas", ".fa", ".fna", ".ffn", ".frn", ".faa", ".fna")
    file_list_query = [i for i in os.listdir(query_file) if
                       os.path.isfile(os.path.join(query_file, i)) and i.endswith(extension)]

    files_list_gold = [i for i in os.listdir(gold_file) if
                       os.path.isfile(os.path.join(gold_file, i)) and i.endswith(extension)]
    query_dict, gold_dict = dict(), dict()
    for i in file_list_query:
        query_dict.update({os.path.splitext(i)[0]: os.path.join(query_file, i)})
    for i in files_list_gold:
        gold_dict.update({os.path.splitext(i)[0]: os.path.join(gold_file, i)})

    shared_gene = [i for i in query_dict.keys() if i in gold_dict.keys()]

    combine_files = []  # [(query_path, gold_path, gene_name)]
    for i in shared_gene:
        combine_files.append((query_dict[i], gold_dict[i], i))
    if combine_files:
        for i in combine_files:
            path1 = i[0]
            path2 = i[1]
            gene_name = i[2]
            path = os.path.join(combined_dir, gene_name + ".fasta")
            with open(path1, "r") as f1, open(path2, "r") as f2, open(path, "a") as f:
                f.writelines(f1.readlines())
                f.writelines(f2.readlines())


def do_alignment_parallel(combined_dir: str, aligned_dir: str, thread_number: int):
    gene_names = [i for i in os.listdir(combined_dir) if os.path.isfile(os.path.join(combined_dir, i))]
    results = []
    with ThreadPoolExecutor(max_workers=thread_number) as executor:
        futures = [executor.submit(do_alignment, combined_file=os.path.join(combined_dir, i),
                                   aligned_file=os.path.join(aligned_dir, i))
                   for i in gene_names]
        for future in as_completed(futures):
            results.append(future.result())


def my_log(log_file: str, info: list):
    with open(log_file, "a", newline='') as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(info)


# [identity,coverage,mutation,insert,deletion,compressed_gap]
def do_parse_alignment(aligned_file: str) -> list:
    #[{id1:seq1},{id2:seq2}]
    seq_list = get_seq(aligned_file, max_seq_number=2, seq_count_limit=True)
    if len(seq_list) < 2:
        print(seq_list)
        print(aligned_file)

    query_seq = list(seq_list[0].values())[0]
    gold_seq = list(seq_list[1].values())[0]
    query_name = list(seq_list[0].keys())[0].replace(" ", "_")
    gold_name = list(seq_list[1].keys())[0].replace(" ", "_")
    gene_name = os.path.splitext(os.path.basename(aligned_file))[0]
    gold_seq_length = len(gold_seq.replace("-", ""))
    query_start_pos = len(query_seq) - len(query_seq.lstrip("-"))
    query_end_pos = len(query_seq.rstrip("-"))
    gold_start_pos = len(gold_seq) - len(gold_seq.lstrip("-"))
    gold_end_pos = len(gold_seq.rstrip("-"))
    start_pos = max(query_start_pos, gold_start_pos)
    end_pos = min(query_end_pos, gold_end_pos)

    # match  mismatch insert deletion mutation compressed_gap
    match_num, mismatch_num, insert, deletion, mutation, compressed_gap = 0, 0, 0, 0, 0, 0

    gap_flag = False
    for i in zip(query_seq[start_pos:end_pos], gold_seq[start_pos:end_pos]):
        if i[0] == "-" and not gap_flag:
            compressed_gap += 1
            gap_flag = True
        elif i[0] != "-" and gap_flag:
            gap_flag = False
        else:
            pass
        if i[0] != i[1]:
            mismatch_num += 1
            if i[0] != "-" and i[1] != "-":
                mutation += 1
            elif i[0] == "-" and i[1] != "-":
                deletion += 1
            else:
                insert += 1
        else:
            match_num += 1
    # coverage = (match_num + mutation) / gold_seq_length
    coverage = format((match_num + mutation) / gold_seq_length * 100, ".2f")
    # gap_compress
    identity = format(match_num / (match_num + mutation + deletion + compressed_gap) * 100, ".2f")

    info = [gene_name, query_name, gold_name, identity, coverage, mutation, deletion, insert, compressed_gap,
            gold_seq_length]
    return info


def do_parse_alignment_parallel(aligned_dir: str, log_file: str, quality_file: str, thread_number: int):
    header = "gene_name,query_name,gold_name,identity,coverage,mutation,deletion,insert,compressed_gap,seq_length"
    with open(log_file, "w") as f:
        f.write(header+"\n")
    aligned_file_list = [os.path.join(aligned_dir, i) for i in os.listdir(aligned_dir) if i.endswith(".fasta") and
                         os.path.isfile(os.path.join(combined_dir, i))]
    results = []
    with ThreadPoolExecutor(max_workers=thread_number) as executor:
        futures = [executor.submit(do_parse_alignment, aligned_file=i)
                   for i in aligned_file_list]
        for future in as_completed(futures):
            results.append(future.result())

    level1, level2, level3, level4, total_mutation, total_deletion, total_compressed_gap, total_seq_len = 0, 0, 0, 0, 0, 0, 0, 0
    for info in results:
        my_log(log_file, info)
        # level1 level2 level3 level4
        # identity coverage mutation deletion insert
        if float(info[3]) == 100 and float(info[4]) == 100:
            level1 += 1
        elif float(info[3]) >= 99 and 100 > float(info[4]) >= 50 and not info[5]:
            level2 += 1
        elif float(info[3]) >= 99 and 50 > float(info[4]) > 0 and not info[5]:
            level3 += 1
        else:
            level4 += 1
        # mutation
        total_mutation += info[5]
        # deletion
        total_deletion += info[6]
        # compressed_gap
        total_compressed_gap += info[8]
        # total_seq_len
        total_seq_len += info[9]
    level1_info = "level1: 100% identity and 100% coverage ,no insert:{}".format(level1)
    level2_info = "level2: 99-100% identity and 50-100% coverage ,no insert:{}".format(level2)
    level3_info = "level3: 99-100% identity and 0~50% coverage ,no insert:{}".format(level3)
    level4_info = "level4: other:{}".format(level4)

    number_info = "mutations: {}, insert: {} ,deletion: {}, seq_length: {}".format(total_mutation, total_compressed_gap,
                                                                                   total_deletion, total_seq_len)
    stat_info = "mutations: {}%, indel: {}%".format(round(total_mutation / total_seq_len * 100, 2),
                                                    round((total_deletion + total_compressed_gap) / total_seq_len * 100,
                                                          2))

    print(
        level1_info + "\n" + level2_info + "\n" + level3_info + "\n" + level4_info + "\n" + number_info + "\n" + stat_info + "\n")
    with open(quality_file, "w") as f:
        f.write(
            level1_info + "\n" + level2_info + "\n" + level3_info + "\n" + level4_info + "\n" + number_info + "\n" + stat_info + "\n")


if __name__ == '__main__':
    t1 = time.time()
    # -i -r filename should be like xx.fasta and be the same. 
	# muscle3 should have been added to PATH	
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="reports",
                                   usage="%(prog)s <-q> <-r> <-o> [options]")

    pars.add_argument("-o", "--out", dest="out_dir", help="Specify the result folder <dir>",
                      metavar="", required=True)
    pars.add_argument("-t", "--thread", metavar="", dest="thread_number", help="Thread", type=int, default=1)
    pars.add_argument("-q", "--query", metavar="", dest="query", type=str, help="query <dir>", required=True)
    pars.add_argument("-r", "--reference", metavar="", dest="gold", type=str, help="gold references <dir>",
                      required=True)
    args = pars.parse_args()
    query_file = absolute_path(args.query)
    gold_file = absolute_path(args.gold)
    out_dir = absolute_path(args.out_dir)
    thread_number = args.thread_number

    paired_dir = os.path.join(out_dir, "paired")
    combined_dir = os.path.join(out_dir, "combined")
    aligned_dir = os.path.join(out_dir, "aligned")
    log_file = os.path.join(out_dir, "align_report.csv")
    quality_file = os.path.join(out_dir, "quality.txt")

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.isdir(paired_dir):
        os.makedirs(paired_dir)
    if not os.path.isdir(combined_dir):
        os.makedirs(combined_dir)
    if not os.path.isdir(aligned_dir):
        os.makedirs(aligned_dir)
    get_best_query_seq_parallel(query_file, gold_file, paired_dir, thread_number)
    do_combine(paired_dir, gold_file, combined_dir)
    do_alignment_parallel(combined_dir, aligned_dir, thread_number)
    do_parse_alignment_parallel(aligned_dir, log_file, quality_file, thread_number)
    t2 = time.time()
    t = format(t2 - t1, ".2f")
    print("Done:{}s".format(t))
