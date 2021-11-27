from numpy.random import choice
list_of_candidates = ['A','T', 'G', 'C']
probability_distribution = [0.4, 0.0, 0.6, 0.0]

for i in range(10):
    draw = choice(list_of_candidates, 1, p = probability_distribution)
    print(draw)