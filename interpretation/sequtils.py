#!/usr/bin/env python

#Some functions implemented in this helper script may not have ended up being used in the final interpretation notebooks

import numpy as np
import gzip
import re
from collections import OrderedDict
import random

def load_motif_matches(motif_match_file, doprint=False):
        motif_matches = OrderedDict()
        fp = open(motif_match_file, "r")
        if doprint:
                print("#Loading " + motif_match_file + " ...")
        numlines = 0
        for line in fp:
                match = re.match("((\w|\-)+)\s+((\w|\:|\-)+)\s+(\d+)\s+(\d+)\s+(\+|\-)\s+.+\s+(\w+)$", line)
                if match:
                        numlines = numlines + 1
                        motif = match.group(1)
                        sequence = match.group(3)
                        begin = int(match.group(5))
                        end = int(match.group(6))
                        strand = match.group(7)
                        seqval = match.group(8)
                        entry = dict()
                        entry['motif'] = motif
                        entry['sequence'] = sequence
                        entry['begin'] = begin-1 # Homer motif match file is 1 indexed, convert to 0
                        entry['end'] = end # Homer motif match file is 1 indexed AND inclusive, convert to 0 and exclusive
                        entry['strand'] = strand
                        entry['seqval'] = seqval
                        if sequence not in motif_matches:
                                motif_matches[sequence] = list()
                        motif_matches[sequence].append(entry)
        fp.close()
        if doprint:
                print("#Loaded " + str(numlines) + " motif matches in " + str(len(motif_matches.keys())) + " sequences")
        return motif_matches

def load_fimo_motif_matches(motif_match_file, doprint=False):
        motif_matches = OrderedDict()
        fp = open(motif_match_file, "r")
        if doprint:
                print("#Loading " + motif_match_file + " ...")
        numlines = 0
        fp.readline()
        for line in fp:
            line = line.split()
            numlines = numlines + 1
            motif = line[0]
            sequence = line[1]
            begin = int(line[2])
            end = int(line[3])
            strand = line[4]
            seqval = line[8]
            entry = dict()
            entry['motif'] = motif
            entry['sequence'] = sequence
            entry['begin'] = begin-1 # Fimo motif match file is 1 indexed, convert to 0
            entry['end'] = end # Fimo motif match file is 1 indexed AND inclusive, convert to 0 and exclusive
            entry['strand'] = strand
            entry['seqval'] = seqval
            if sequence not in motif_matches:
                motif_matches[sequence] = list()
            motif_matches[sequence].append(entry)
        fp.close()
        if doprint:
                print("#Loaded " + str(numlines) + " motif matches in " + str(len(motif_matches.keys())) + " sequences")
        return motif_matches

def rename(label):
    match=re.match('.*_(chr.*)$',label)
    if match:
        return match.group(1)
    else:
        return ""

def renameAll(x):
    new_x=[]
    for element in x:
        new_x.append(rename(element))
    return new_x

def get_relevant_labels_in_order_of_scores(labels, motif_matches):
    relevant_labels_list=[]
    relevant_indices_list=[]
    sequence_index=0
    positive_labels=[]
    for label in motif_matches.keys():
        positive_labels.append(label)
    positive_labels_set = set(positive_labels)
    print("Motif matches sequences are " + str(len(positive_labels_set)))
    print("Supplied labels are " + str(len(labels)))
    for sequence_label in labels:
        if sequence_label in positive_labels_set:
            relevant_indices_list.append(sequence_index)
            relevant_labels_list.append(sequence_label)
        else:
            print("Did not find this label in motif matches: " + sequence_label)
        sequence_index=sequence_index+1
    #print (len(relevant_indices_list))
    #print (len(relevant_labels_list))
    return relevant_indices_list, relevant_labels_list

def get_relevant_scores(relevant_indices_list, scores, seq_len=400):
    relevant_scores=np.zeros((len(relevant_indices_list),seq_len))
    index=0
    for scores_index in relevant_indices_list:
        relevant_scores[index]=scores[scores_index]
        index=index+1
    return relevant_scores

def removeUnsupportedChars(sequences, labels, labeled_sequences):
    removed=[]
    chars=['R','Y','S','W','K','M','B','D','H','V','N']
    print(len(sequences))
    for seq in sequences:
        if any((c in chars) for c in seq):
            print(seq)
            removed.append(seq)
            sequences.remove(seq)
    print(len(sequences))
    for i in removed:
        key=labeled_sequences.keys()[labeled_sequences.values().index(i)]
        print(key)
        del labeled_sequences[key]
        labels.remove(key)
    print (len(labels))
    print(len(labeled_sequences))

def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return

def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: " + str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1

def get_random_set(seqdict, num, sort=True):
    newlist=seqdict.items()
    if sort:
        newlist = [newlist[i] for i in sorted(random.sample(range(len(newlist)), num))]
    else:
        newlist = [newlist[i] for i in random.sample(range(len(newlist)), num)]
    return dict(newlist)

def get_value(label):
    value = -1
    match = re.match("((dinuc_shuff_|dinuc_shuffled_).+)$", label)
    if match:            
        chrom = match.group(1)
        if match.group(2) == 'dinuc_shuff_':
            value = 0
        else:
            value = 1
    return value 

def load_sequences_from_bedfile(seqfile):
    seqs = OrderedDict()
    fp = gzip.open(seqfile, "rb")
    print("#Loading " + seqfile + " ...")
    for line in fp:
        parsed = (line.decode()).split()
        seqs[parsed[0]]=parsed[1]
    fp.close()
    print("#Loaded " + str(len(seqs.keys())) + " sequences from " + seqfile)
    return seqs

def count_seqs(seqfile):
    return len(load_sequences_from_bedfile(seqfile))

def load_sequences(seqfile):
        seqs = OrderedDict()
        fp = open(seqfile, "rb")
        print("#Loading " + seqfile + " ...")
        expecting = "label"
        label=''
        for line in fp:
                if expecting == "label":
                        match = re.match(">(.*)$", line)
                        if match:
                                label = match.group(1)
                                expecting = "sequence"
                        else:
                                print("Expecting LABEL but found (!!): " + line)
                                continue
                else:
                        match = re.match("(\w+)$", line)
                        if match:
                                sequence = match.group(1)
                                seqs[label]=sequence
                        else:
                                print("Expecting SEQUENCE but found (!!): " + line)
                        expecting = "label"
                        label=''
        fp.close()
        print("#Loaded " + str(len(seqs.keys())) + " sequences from " + seqfile)
        return seqs

def load_labels(labels_file):
    f = open(labels_file, "r")
    l = f.read().splitlines()
    print ("Read " + str(len(l)) + " labels from " + str(labels_file))
    f.close()
    return l

def save_labels(labeled_sequences, filename):
    print(type(labeled_sequences))
    with open(filename, "w") as ff:
        for ss in labeled_sequences.keys():
            ff.write(str(ss) +"\n")
    ff.close()
