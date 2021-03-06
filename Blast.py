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

def create_query_file(filename, fileType):
    fname = ''
    if fileType == "Control":
        fname = "ControlQueries/"+filename+".fasta"
    elif fileType == "Mutation":
        fname = "MutationQueries/"+filename+".fasta"
    elif fileType == "Indels":
        fname = "IndelQueries/"+filename+".fasta" 
    qfile = open(fname, "w")
    qfile.close()

# Create query .fasta file using randomized database (no mutations, no indels). Control testing BLAST (query with identical nucleotides).
# Params: .fasta filename (string), query number (int), dbseq database sequence from random_seq (list), cap query length at 20 nucleotides (bool)
def create_query(filename, qnum, dbseq, capLength):
    fname = "ControlQueries/"+filename+".fasta" #path might need to change depending on user
    qname = "Control query "+str(qnum)
    qfile = open(fname, "a")

    random_startpos = random.randint(0, len(dbseq)-3)
    random_endpos = 0
    if (random_startpos == (len(dbseq)-3)):
        random_endpos = len(dbseq)-1 # query needs to be at least 2 nucleotides long
    else:
        random_endpos = random.randint(random_startpos+2, len(dbseq)-1)

    endPos = random_startpos + 20
    if capLength == True and endPos<len(dbseq):
        random_endpos = endPos
        
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
# Params: .fasta filename (string), query number (int), number of mutations (int - at least 1), dbseq database sequence from random_seq (list), cap query length at 20 nucleotides (bool)
def mut_query(filename, qnum, mutNum, dbseq, capLength):
    fname = "MutationQueries/"+filename+".fasta" #path might need to change depending on user
    qname = "Mutation query "+str(qnum)
    qfile = open(fname, "a")
    
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
    
    endPos = random_startpos + 20
    if capLength == True and endPos<len(dbseq):
        random_endpos = endPos

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
        mutPos = random.randint(0, len(tmp_query)-1)
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
# Params: .fasta filename (string), query number (int), number of indels (int - at least 1), dbseq database sequence from random_seq (list), cap query length at 20 nucleotides (bool)
def indel_query(filename, qnum, indelNum, dbseq, capLength):
    fname = "IndelQueries/"+filename+".fasta" #path might need to change depending on user
    qname = "Indel query "+str(qnum)
    qfile = open(fname, "a")
    
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
    
    endPos = random_startpos + 20
    if capLength == True and endPos<len(dbseq):
        random_endpos = endPos

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
    cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=filename, )
    cline()

# Running BLAST on specified database with specific query
# Params: .fasta query filename (q - string), database name (d - string), .xml outfile name (o - string)
def run_BLASTn(q, d, o):
    blastn_cline = NcbiblastnCommandline(query=q, db=d, \
    evalue=0.05, outfmt=5, out=o, word_size=7, task='blastn')
    print(blastn_cline)
    stdout, stderr = blastn_cline()

# Calculates probability of matched sequence database
# Params: database sequence (list), database probability matrix (list(list)), start match pos (int), end match pos (int)
def calculate_probs(db, dbprobs, start, end):
    final_prob = 1
    for i in range(start-1, end):
        if db[i] == 'A':
            final_prob = final_prob * dbprobs[i][0]
        elif db[i] == 'T':
            final_prob = final_prob * dbprobs[i][1]
        elif db[i] == 'G':
            final_prob = final_prob * dbprobs[i][2]
        elif db[i] == 'C':
            final_prob = final_prob * dbprobs[i][3]
    return final_prob

# Parsing XML file
# Param: XML filename (string), results filename (string), database probabilities (list(list)), database sequence (list)
def parse_XML(filename, resultsfile, dbprobs, db):
    fname = "Results/"+resultsfile+".txt" #path might need to change depending on user
    result = open(fname, "a")
    for record in NCBIXML.parse(open(filename)):
        #print(record.alignments)
        if record.alignments: #skip queries with no matches
            #print("Hi")
            q = "QUERY: %s\n" % record.query[:60]
            print("QUERY: %s" % record.query[:60])
            result.write(q)
            for align in record.alignments:
                for hsp in align.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        m = "MATCH: %s\n" % align.title[:60]
                        print("MATCH: %s " % align.title[:60])
                        result.write(m)
                        
                        start = int(hsp.sbjct_start)
                        end = int(hsp.sbjct_end)
                        datab_prob = calculate_probs(db, dbprobs, start, end)
                        print("Match database sequence probability: " + str(datab_prob) + '\n')
                        result.write("Match database sequence probability: " + str(datab_prob) + '\n')
                        
                        h = str(hsp)+'\n'
                        print(hsp)
                        result.write(h)
    result.close()
                    

# Permutates for numPerms amount of times control queries (not permutation for BLAST)
def permutate_control(numPerms, datab):
    create_query_file('controls', "Control")
    for i in range(numPerms):
        create_query('controls', i, datab, True)

# Permutates for numPerms amount of times mutation queries (not permutation for BLAST)
def permutate_mut(numPerms, datab):
    create_query_file('mutations', "Mutation")
    for i in range(numPerms):
        numMuts = random.randint(1,5)
        mut_query('mutations', i, numMuts, datab, True)

# Permutates for numPerms amount of times mutation queries (not permutation for BLAST)
def permutate_indel(numPerms, datab):
    create_query_file('indels', "Indels")
    for i in range(numPerms):
        numIndels = random.randint(1,8)
        indel_query('indels2', i, numIndels, datab, True)

# Permutates datab for a same control query        
def permutate_datab(numPerms, seqProbs):
    datab = random_seq(seqProbs)
    dbname = create_fasta_datab(datab, 'permutedb50', 1)
    create_datab(dbname)

    create_query_file("control_for_db_perm50", "Control")
    create_query("control_for_db_perm50", 1, datab, True)

    result = open("Results/control_for_db_perm50_results.txt", "w")

    for i in range(numPerms):
        result.write("DATABASE"+str(i)+'\n')
        run_BLASTn("ControlQueries/control_for_db_perm50.fasta", dbname, 'controldb50_out.xml')
        parse_XML('controldb50_out.xml','control_for_db_perm50_results', seqProbs, datab)

        datab = random_seq(seqProbs)
        dbname = create_fasta_datab(datab, 'permutedb50', 1)
        create_datab(dbname)

    result.close()

# Permutates datab for a same mutation query        
def permutate_datab_mut(numPerms, seqProbs):
    datab = random_seq(seqProbs)
    dbname = create_fasta_datab(datab, 'permutedb10-mut', 1)
    create_datab(dbname)

    create_query_file("mut_for_db_perm10", "Mutation")
    mut_query("mut_for_db_perm10", 1, 4, datab, True)

    result = open("Results/mut_for_db_perm10_results.txt", "w")

    for i in range(numPerms):
        result.write("DATABASE"+str(i)+'\n')
        run_BLASTn("MutationQueries/mut_for_db_perm10.fasta", dbname, 'mutdb10_out.xml')
        parse_XML('mutdb10_out.xml','mut_for_db_perm10_results', seqProbs, datab)

        datab = random_seq(seqProbs)
        dbname = create_fasta_datab(datab, 'permutedb10-mut', 1)
        create_datab(dbname)

    result.close()

# Permutates datab for a same indel query        
def permutate_datab_indel(numPerms, seqProbs):
    datab = random_seq(seqProbs)
    dbname = create_fasta_datab(datab, 'permutedb10-indel', 1)
    create_datab(dbname)

    create_query_file("indel_for_db_perm10", "Indels")
    indel_query("indel_for_db_perm10", 1, 2, datab, True)

    result = open("Results/indel_for_db_perm10_results.txt", "w")

    for i in range(numPerms):
        result.write("DATABASE"+str(i)+'\n')
        run_BLASTn("IndelQueries/indel_for_db_perm10.fasta", dbname, 'indeldb10_out.xml')
        parse_XML('indeldb10_out.xml','indel_for_db_perm10_results', seqProbs, datab)

        datab = random_seq(seqProbs)
        dbname = create_fasta_datab(datab, 'permutedb10-indel', 1)
        create_datab(dbname)

    result.close()

# TESTS
def test_controlq():
    seq = sequence_probs('full_probs.txt', 'full_seq.txt')
    y = random_seq(seq)
    permutate_control(3, y)
    #print(y)
    #create_query('control3', 1, y, True)
    #mut_query('test', 1, 2, y)
    #indel_query('test', 1, 3, y)
    dbname = create_fasta_datab(y, 'controldb', 1)
    create_datab(dbname)
    run_BLASTn("ControlQueries/controls.fasta", dbname, 'control_out.xml')
    parse_XML('control_out.xml','control_results', seq, y)

def test_mutq():
    seq = sequence_probs('full_probs.txt', 'full_seq.txt')
    y = random_seq(seq)
    permutate_mut(5, y)
    dbname = create_fasta_datab(y, 'mutdb', 1)
    create_datab(dbname)
    run_BLASTn("MutationQueries/mutations.fasta", dbname, 'mut_out.xml')
    parse_XML('mut_out.xml','mut_results', seq, y)

def test_indelq():
    seq = sequence_probs('full_probs.txt', 'full_seq.txt')
    y = random_seq(seq)
    permutate_indel(5, y)
    dbname = create_fasta_datab(y, 'indeldb2', 1)
    create_datab(dbname)
    run_BLASTn("IndelQueries/indels2.fasta", dbname, 'indel_out2.xml')
    parse_XML('indel_out2.xml','indel_results2', seq, y)

def test_db_perms():
    seq = sequence_probs('full_probs.txt', 'full_seq.txt')
    permutate_datab(50, seq)

def test_db_perms_mut():
    seq = sequence_probs('full_probs.txt', 'full_seq.txt')
    permutate_datab_mut(10, seq)

def test_db_perms_indel():
    seq = sequence_probs('full_probs.txt', 'full_seq.txt')
    permutate_datab_indel(10, seq)

#test_controlq()
#test_mutq()
#test_indelq()
test_db_perms()
#test_db_perms_mut()
#test_db_perms_indel()