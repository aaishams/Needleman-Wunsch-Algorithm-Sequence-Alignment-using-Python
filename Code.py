'''
THIS IS MY FIRST BIOINFORMATICS PROJECT USING PYTHON.

This code can be used for the sequence alignment of 2 sequences using Needleman-Wunsch Algorithm.

HOPE IT MAKE YOUR BIOINFORMATICS BASED ACTIVITIES EFFICIENT!! 

'''
# define the function to alert if the input sequence is not a nucleotide
def is_nucleotide(seq):
    nucleotides = set('ATGC')
    return all(nuc in nucleotides for nuc in seq)

# define the function to calculate the alignment score
def alignment_score(alignment1, alignment2, match_score, mismatch_score, gap_penalty):
    score = 0
    for char1, char2 in zip(alignment1, alignment2):
        if char1 == '-' or char2 == '-':
            score += gap_penalty
        elif char1 == char2:
            score += match_score
        else:
            score += mismatch_score
    return score

# define the function to display the matrix
def display_matrix(matrix, seq1, seq2):
    print("- - " + " ".join(seq2))
    for i in range(len(matrix)):
        if i == 0:
            print("- ", end=" ")
        else:
            print(seq1[i - 1], end=" ")
        for j in range(len(matrix[i])):
            print(matrix[i][j], end=" ")
        print()

# define the function
def needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty):
    # initialise the matrix
    matrix = [[0 for i in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]

    # fill the penalties
    for i in range(1, len(seq1) + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(1, len(seq2) + 1):
        matrix[0][j] = j * gap_penalty

    # fill the scoring matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            a = matrix[i - 1][j] + gap_penalty
            b = matrix[i][j - 1] + gap_penalty
            c = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            matrix[i][j] = max(a, b, c)

    # traceback to find the optimal alignment by recursion
    align1 = ''
    align2 = ''
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        elif j > 0 and matrix[i][j] == matrix[i][j - 1] + gap_penalty:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1
        else:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2 
            i -= 1
            j -= 1

    return align1, align2, matrix

# get the user input
seq1 = input("Enter the first (Query) sequence:")
seq1 = seq1.upper()
if not is_nucleotide(seq1):
    print("Warning! The given sequence is not a nucleotide..")
seq2 = input("Enter the second (Subject) sequence:")
seq2 = seq2.upper()
if not is_nucleotide(seq2):
    print("Warning! The given sequence is not a nucleotide..")
match_score = int(input("Enter the match score:"))
mismatch_score = int(input("Enter the mismatch score:"))
gap_penalty = int(input("Enter the gap penalty:"))

# get and display the output
alignment1, alignment2, final_matrix = needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty)
display_matrix(final_matrix, seq1, seq2)
print("Alignment 1:", alignment1)
print("Alignment 2:", alignment2)
score = alignment_score(alignment1, alignment2, match_score, mismatch_score, gap_penalty)
print("Alignment score:", score)