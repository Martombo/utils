def random_seq(length):
    "returns a random nucleotide sequence of the desired length"
    import random as rn

    seq = ''
    for k in range(length):
        rand = rn.random()
        if rand <= 0.25:
            next_letter = 'A'
        elif rand <= 0.5:
            next_letter = 'C'
        elif rand <= 0.75:
            next_letter = 'T'
        else:
            next_letter = 'G'

        seq += next_letter

    return seq