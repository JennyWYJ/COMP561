# Second attempt at alignment using BLAST (more time and space efficient)
import numpy as np
from numpy.random import choice
import random

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast import NCBIXML
E_VALUE_THRESH = 1000#1e-20

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

# ALIGNMENT SECTION BELOW USING BLAST

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

# Create query .fasta file containing mutations using randomized database.
# Params: .fasta filename (string), query number (int), number of indels (int - at least 1), dbseq database sequence from random_seq (list)
def indel_query(filename, qnum, indelNum, dbseq):
    fname = "IndelQueries/"+filename+".fasta" #path might need to change depending on user
    qname = "Indel query "+str(qnum)
    qfile = open(fname, "w")
    
    numIndels = indelNum
    list_of_candidates = ['A','T', 'G', 'C']
    prob_distribution = [0.25, 0.25, 0.25, 0.25]
    indels = ['I', 'D'] #insert or delete
    indelProbs = [0.5, 0.5]

    random_startpos = random.randint(0, len(dbseq)-3)
    random_endpos = 0
    if (random_startpos == (len(dbseq)-3)):
        random_endpos = len(dbseq)-1 # query needs to be at least 2 nucleotides long
    else:
        random_endpos = random.randint(random_startpos+2, len(dbseq)-1)

    numNucs = random_endpos - random_startpos
    if numIndels >= numNucs//2: #suppose we have single indels, so can't be more than half of the og query
        if numNucs//2 > 0 and numNucs != 3:
            numIndels = numNucs//2
        else:
            numIndels = 1
    
    tmp_query = []
    if random_endpos == len(dbseq)-1:
        qfile.write(">" + qname + "\n")
        for i in range(random_startpos,len(dbseq)):
            tmp_query.append(dbseq[i])
    else:
        qfile.write(">" + qname + "\n")
        for i in range(random_startpos,random_endpos+1):
            tmp_query.append(dbseq[i])

    for m in range(numIndels):
        indelPos = random.randint(0, len(tmp_query))
        inOrDel = choice(indels, 1, p = indelProbs)
        inNuc = choice(list_of_candidates, 1, p = prob_distribution)

        if inOrDel[0] == 'I':
            tmp_query = tmp_query[:indelPos] + [inNuc[0]] + tmp_query[indelPos:]
        elif inOrDel[0] == 'D':
            if indelPos+1 < len(tmp_query):
                tmp_query = tmp_query[:indelPos] + tmp_query[indelPos+1:]
            elif indelPos+1 == len(tmp_query):
                tmp_query = tmp_query[:indelPos]

    for nu in tmp_query:
        qfile.write(nu)
    qfile.write("\n")

    qfile.close()

# Create fasta file for generated database
# Params: generated database nucleotides (list), .fasta filename (string), database number (int)
# Returns: .fasta filename
def create_fasta_datab(generated_db, fname, dbNum):
    dbname = "Databases/"+fname+".fasta" #path might need to change depending on user
    dbfile = open(dbname, "w")

    dbfile.write(">Database " + str(dbNum) + "\n")
    for n in generated_db:
        dbfile.write(n)
    dbfile.write("\n")

    dbfile.close()
    return dbname

# Creating a BLAST database using database .fasta filename (string)
# Param: .fasta filename
# Returns: custom database
def create_datab(filename):
    cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=filename)
    return cline()

# Running BLAST on specified database with specific query
# Params: .fasta query filename (q - string), database name (d - string), .xml outfile name (o - string)
def run_BLASTn(q, d, o):
    blastn_cline = NcbiblastnCommandline(query=q, db=d, \
    evalue=1e-20, outfmt=5, out=o)
    return blastn_cline()

# Parsing XML file
# Param: XML filename (string)
def parse_XML(filename):
    for record in NCBIXML.parse(open(filename)):
        print(record.alignments)
        if record.alignments: #skip queries with no matches
            print("Hi")
            print("QUERY: %s" % record.query[:60])
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print("MATCH: %s " % align.title[:60])
                        print(hsp.expect)

# TESTS
seq = sequence_probs('probsa.txt', 'seqa.txt')
y = random_seq(seq)
print(y)
create_query('test', 1, y)
#mut_query('test', 1, 2, y)
#indel_query('test', 1, 3, y)
dbname = create_fasta_datab(y, 'test', 1)
#d = create_datab(dbname)
run_BLASTn("ControlQueries/test.fasta", dbname, 'out.xml')
parse_XML('out.xml')