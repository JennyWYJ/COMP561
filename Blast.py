# Second attempt at alignment using BLAST (more time and space efficient)
import numpy as np
from numpy.random import choice
import random

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline

import subprocess


# Testing Bio.BLAST. "test.fasta" can be changed to whatever fasta filename we set for query. 
# Database here is solely NCBI database.
def test_BLAST():
    my_query = SeqIO.read("test.fasta", format="fasta")
    result_handle = NCBIWWW.qblast("blastn", "nt", my_query.seq)
    blast_result = open("my_blast.xml", "w")
    blast_result.write(result_handle.read())
    blast_result.close()
    result_handle.close()

# Initialize (returns) database nucleotide probability list. 
def sequence_probs(probsfile, seqfile):
    pf = open(probsfile)
    sf = open(seqfile)

    n = 0

    lines = pf.readlines()
    probs = []
    for line in lines:
        l = line.split()
        for prob in l:
            if prob != '\n':
                probs.append(float(prob))

    seqs = []
    while True:
        c = sf.read(1)
        nprobs = [0,0,0,0]
        if not c:
            print("End of file")
            break
        if n > len(probs)-1:
            break
        if c == 'A':
            nprobs[0] = probs[n]
            other_probs = (1.0-probs[n])/3
            nprobs[1] = other_probs
            nprobs[2] = other_probs
            nprobs[3] = other_probs
            seqs.append(nprobs)

        elif c == 'T':
            nprobs[1] = probs[n]
            other_probs = (1.0-probs[n])/3
            nprobs[0] = other_probs
            nprobs[2] = other_probs
            nprobs[3] = other_probs
            seqs.append(nprobs)

        elif c == 'G':
            nprobs[2] = probs[n]
            other_probs = (1.0-probs[n])/3
            nprobs[0] = other_probs
            nprobs[1] = other_probs
            nprobs[3] = other_probs
            seqs.append(nprobs)

        elif c == 'C':
            nprobs[3] = probs[n]
            other_probs = (1.0-probs[n])/3
            nprobs[0] = other_probs
            nprobs[1] = other_probs
            nprobs[2] = other_probs
            seqs.append(nprobs)

        n+=1

    pf.close()
    sf.close()

    return seqs

# Semi-randomly selects nucleotides based probabilities list (seq_probs) to form database sequence for alignment.
def random_seq(seq_probs):
    list_of_candidates = ['A','T', 'G', 'C']
    seq_copy = []
    for i in range(len(seq_probs)):
        probability_distribution = seq_probs[i]
        draw = choice(list_of_candidates, 1, p = probability_distribution)
        seq_copy.append(draw[0])
    return(seq_copy)

# Get the list of probabilities depending on a specified seqence (seq) and the probability matrix (based on specific nucleotides).
def datab_probs(seq, prob):

    probs_list = []
    for i, n in enumerate(seq):
        
        if n == 'A':
            probs_list.append(prob[i][0])
        elif n == 'T':
            probs_list.append(prob[i][1])
        elif n == 'G':
            probs_list.append(prob[i][2])
        elif n == 'C':
            probs_list.append(prob[i][3])

    return probs_list

# Alignment using BLAST

# Running BLAST on specified database with specific query
# Params: .fasta query filename (q - string), database name (d - string), .xml outfile name (o - string)
def run_BLASTn(q, d, o):
    blastn_cline = NcbiblastnCommandline(query=q, db=d, \
    evalue=1e-20, outfmt=5, out=o)
    stdout, stderr = blastn_cline()

# Create query .fasta file using randomized database (no mutations, no indels). Control testing BLAST (query with identical nucleotides).
# Params: .fasta filename (string), query number (int), dbseq database sequence from random_seq (list)
def create_query(filename, qnum, dbseq):
    fname = "ControlQueries/"+filename+".fasta" #path might need to change depending on user
    qname = "Control query "+str(qnum)
    qfile = open(fname, "w")

    random_startpos = random.randint(0, len(dbseq)-3)
    random_endpos = 0
    if (random_startpos == (len(dbseq)-3)):
        random_endpos = len(dbseq)-1 # query needs to be at least 2 nucleotides long
    else:
        random_endpos = random.randint(random_startpos+2, len(dbseq)-1)

    if random_endpos == len(dbseq)-1:
        qfile.write(">" + qname + "\n")
        for i in range(random_startpos,len(dbseq)):
            qfile.write(dbseq[i])
        qfile.write("\n")
    else:
        qfile.write(">" + qname + "\n")
        for i in range(random_startpos,random_endpos+1):
            qfile.write(dbseq[i])
        qfile.write("\n")

    qfile.close()

# Create query .fasta file containing mutations using randomized database.
# Params: .fasta filename (string), query number (int), number of mutations (int - at least 1), dbseq database sequence from random_seq (list)
def mut_query(filename, qnum, mutNum, dbseq):
    fname = "MutationQueries/"+filename+".fasta" #path might need to change depending on user
    qname = "Mutation query "+str(qnum)
    qfile = open(fname, "w")
    
    numMuts = mutNum
    list_of_candidates = ['A','T', 'G', 'C']

    random_startpos = random.randint(0, len(dbseq)-3)
    random_endpos = 0
    if (random_startpos == (len(dbseq)-3)):
        random_endpos = len(dbseq)-1 # query needs to be at least 2 nucleotides long
    else:
        random_endpos = random.randint(random_startpos+2, len(dbseq)-1)

    numNucs = random_endpos - random_startpos
    if numMuts >= numNucs:
        numMuts = numNucs - 1 #suppose at least one nucleotide is not corrupted
    
    tmp_query = []
    if random_endpos == len(dbseq)-1:
        qfile.write(">" + qname + "\n")
        for i in range(random_startpos,len(dbseq)):
            tmp_query.append(dbseq[i])
    else:
        qfile.write(">" + qname + "\n")
        for i in range(random_startpos,random_endpos+1):
            tmp_query.append(dbseq[i])
    
    for m in range(numMuts):
        mutPos = random.randint(0, len(tmp_query))
        nAtPos = tmp_query[mutPos]
        prob_distribution = []
        for n in list_of_candidates:
            if nAtPos == n:
                prob_distribution.append(0)
            else:
                prob_distribution.append(1/3)
        mutNuc = choice(list_of_candidates, 1, p = prob_distribution)
        tmp_query[mutPos] = mutNuc[0]

    for nu in tmp_query:
        qfile.write(nu)
    qfile.write("\n")

    qfile.close()

# Creating a BLAST database by calling bash script
def create_datab(script_path):
    subprocess.call(script_path)

#create_datab('https://github.com/JennyWYJ/COMP561/createDB.sh')
seq = sequence_probs('probsa.txt', 'seqa.txt')
y = random_seq(seq)
print(y)
#create_query('test', 1, y)
#mut_query('test', 1, 2, y)