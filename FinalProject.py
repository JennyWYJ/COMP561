#Code for COMP 561 Final Project initial attempt (using S-W algorithm)
import numpy as np
from numpy.random import choice

MATCH = 1
MISMATCH = -1
GAP = -1

STOP = 0
LEFT = 1 
UP = 2
DIAGONAL = 3

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

# Semi-randomly selects nucleotides based probabilities to form database sequence for alignment.
def random_seq(seq_probs):
    list_of_candidates = ['A','T', 'G', 'C']
    seq_copy = []
    for i in range(len(seq_probs)):
        probability_distribution = seq_probs[i]
        draw = choice(list_of_candidates, 1, p = probability_distribution)
        seq_copy.append(draw[0])
    return(seq_copy)  

def best_prob(probs):
    seq_copy = []
    list_of_candidates = ['A','T', 'G', 'C']
    for set in probs:
        best = 0
        index = ''
        for i in range(len(set)):
            if set[i] > best:
                best = set[i]
                index = i
        seq_copy.append(list_of_candidates[index])
    #print(seq_copy)
    return seq_copy

def best_score(probs):
    seq_copy = []
    list_of_candidates = ['A','T', 'G', 'C']
    for set in probs:
        if 1.0 in set:
            i = set.index(1.0)
            seq_copy.append(list_of_candidates[i])
        else:
            seq_copy.append('M')
    print(seq_copy)
    return seq_copy




def datab_probs(seq, prob):

    probs_list = []
    #print(prob)
    #print(seq)
    for i, n in enumerate(seq):
        #print(prob[i][0])
        
        if n == 'A':
            probs_list.append(prob[i][0])
        elif n == 'T':
            probs_list.append(prob[i][1])
        elif n == 'G':
            probs_list.append(prob[i][2])
        elif n == 'C':
            probs_list.append(prob[i][3])

    return probs_list


# Alignment
def smith_waterman(query, datab, probs_list):

    row = len(query) + 1
    col = len(datab) + 1
    scores = np.zeros((row, col))
    trace = np.zeros((row, col))  

    for pos, n in enumerate(datab):
        if isinstance(n, list):
            for r in range(0, len(query)+1):
                scores[r][pos] = [0, 0, 0, 0]
    
    max_score = -1
    max_index = (-1, -1)
    
    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            if datab[j-1] == 'M':
                match_value = MATCH
            elif query[i - 1] == datab[j - 1]:
                match_value = MATCH 
            else:
                match_value = MISMATCH
            diagonal_score = scores[i - 1, j - 1] + match_value
            
            # Calculating the vertical gap score
            vertical_score = scores[i - 1, j] + GAP
            
            # Calculating the horizontal gap score
            horizontal_score = scores[i, j - 1] + GAP
            
            # Taking the highest score 
            scores[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)
            
            # Tracking where the cell's value is coming from    
            if scores[i, j] == 0: 
                trace[i, j] = STOP
                
            elif scores[i, j] == horizontal_score: 
                trace[i, j] = LEFT
                
            elif scores[i, j] == vertical_score: 
                trace[i, j] = UP
                
            elif scores[i, j] == diagonal_score: 
                trace[i, j] = DIAGONAL 
                
            # Tracking the cell with the maximum score
            if scores[i, j] >= max_score:
                max_index = (i,j)
                max_score = scores[i, j]
    
    # Initialising the variables for tracing
    aligned_query = ""
    aligned_datab = ""   
    current_aligned_query = ""   
    current_aligned_datab = ""  
    (max_i, max_j) = max_index
    total_prob = 1
    # Tracing and computing the pathway with the local alignment
    nucleotide = ['A','T', 'G', 'C']
    while trace[max_i, max_j] != STOP:
        if trace[max_i, max_j] == DIAGONAL:
            current_aligned_query = query[max_i - 1]
            if datab[max_j - 1] == 'M':
                current_aligned_datab = query[max_i - 1]
            else:
                current_aligned_datab = datab[max_j - 1]
            if isinstance(probs_list[max_j - 1], list):
                index = nucleotide.index(query[max_i - 1])
                total_prob = total_prob*probs_list[max_j - 1][index]
            else:
                total_prob = total_prob*probs_list[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1
            
        elif trace[max_i, max_j] == UP:
            current_aligned_query = query[max_i - 1]
            current_aligned_datab = '-'
            max_i = max_i - 1    
            
        elif trace[max_i, max_j] == LEFT:
            current_aligned_query = '-'
            current_aligned_datab = datab[max_j - 1]
            if isinstance(probs_list[max_j - 1], list):
                index = nucleotide.index(datab[max_j - 1])
                total_prob = total_prob*probs_list[max_j-1][index]
            else:
                total_prob = total_prob*probs_list[max_j - 1]
            max_j = max_j - 1
            
        aligned_query = aligned_query + current_aligned_query
        aligned_datab = aligned_datab + current_aligned_datab
    
    # Reversing the order of the sequences
    aligned_query = aligned_query[::-1]
    aligned_datab = aligned_datab[::-1]
    
    return aligned_query, aligned_datab, max_score, total_prob



seq = sequence_probs('probsa.txt', 'seqa.txt')
best = best_score(seq)
print(best)
#y = random_seq(seq)
#print(seq)
#print(y)
#z = datab_probs(y,seq)
#z = datab_probs(best,seq)
#print(z)
query = ['A', 'C', 'A', 'T']
print(smith_waterman(query, best, seq))