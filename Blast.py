# Second attempt at alignment using BLAST (more time and space efficient)
import numpy as np
from numpy.random import choice

# Initialize (returns) query nucleotide probability list. 
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