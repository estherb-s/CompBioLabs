import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Defines the order of the amino acids used for the score matrix
alphabet = "ARNDCQEGHILKMFPSTWYV"

# BLOSUM50 matrix stored as a 2D array
blosum50 = np.array([
    [ 5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0],
    [-2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3],
    [-1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3],
    [-2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4],
    [-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],
    [-1, 1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3],
    [-1, 0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3],
    [ 0,-3, 0,-1,-3,-2,-3, 8,-2,-4,-4,-2,-3,-4,-2, 0,-2,-3,-3,-4],
    [-2, 0, 1,-1,-3, 1, 0,-2,10,-4,-3, 0,-1,-1,-2,-1,-2,-3, 2,-4],
    [-1,-4,-3,-4,-2,-3,-4,-4,-4, 5, 2,-3, 2, 0,-3,-3,-1,-3,-1, 4],
    [-2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,-3, 3, 1,-4,-3,-1,-2,-1, 1],
    [-1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,-2,-4,-1, 0,-1,-3,-2,-3],
    [-1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7, 0,-3,-2,-1,-1, 0, 1],
    [-3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,-4,-3,-2, 1, 4,-1],
    [-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3],
    [ 1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5, 2,-4,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,-3,-2, 0],
    [-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15, 2,-3],
    [-2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,-1],
    [ 0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5],
])

#Linear gap penalty of -8
gap_penalty = -4

# The two sequences to align
X = "HEAGAWGHEE"
Y = "PAWHEAE"

def needlemanWunsch(X,Y):

    # 2D Score matrix length of the sequences +1
    score_matrix = np.zeros((len(X)+1, len(Y)+1), dtype=int)      

    for i, Xa in enumerate(X):
        score_matrix[i+1, 0] = gap_penalty*(i+1)
    for j, Ya in enumerate(Y):
        score_matrix[0, j+1] = gap_penalty*(j+1)
        
    for j, Ya in enumerate(Y):
        for i, Xa in enumerate(X):
            # index (position) of X and Y amino acids in the alphabet
            Xai = alphabet.index(Xa)
            Yai = alphabet.index(Ya)
            
            # indices to X and Y pos
            fi = i + 1
            fj = j + 1

            # Match score for aligning X and Y
            match_score = blosum50[Xai, Yai]

            # The two sequences are aligned         
            match = score_matrix[fi-1, fj-1] + match_score  # (fi-1, fj-1) = diagonally up & left cell
            # X is aligned with a gap         
            delete = score_matrix[fi-1, fj] + gap_penalty    # (fi-1, fj) = cell to the left of the current cell
            # Y is aligned with a gap         
            insert = score_matrix[fi, fj-1] + gap_penalty  # (fi, fj - 1) = cell above the current cell
            # This gives the max score among all potential alignments
            max_score = max(match, delete, insert)
            score_matrix[fi, fj] = max_score
    
    # Compute the alignment
    alignX, alignY = '', ''
    # Start from the bottom right cell
    i,j = len(X), len(Y)
    while i>0 or j>0:
        Xai = alphabet.index(X[i - 1])
        Yai = alphabet.index(Y[j - 1])
        match_Score = blosum50[Xai, Yai]
        print(i, j, "Match score:", match_Score)
        # Align 2 residues         
        if i > 0 and j >0 and score_matrix[i, j] == score_matrix[i-1, j-1] + match_Score:
            alignX += X[i-1]
            alignY += Y[j-1]
            i -= 1
            j -= 1
        # Align X with a gap             
        elif i > 0 and score_matrix[i, j] == score_matrix[i-1, j] + gap_penalty:
            alignX += X[i-1]
            alignY += '-'
            i -= 1
        # Align Y with a gap             
        elif j > 0:
            alignX += '-'
            alignY += Y[j-1]
            j -= 1
            
    while i >0:
        alignX += X[i-1]
        alignY += '-'
        i -= 1
    while j > 0:
        alignX += '-'
        alignY += Y[j-1]
        j-=1
        
    alignX = alignX[::-1]
    alignY = alignY[::-1]
    return(alignX, alignY)

alignX, alignY = needlemanWunsch(X, Y)
print("\n" + alignX + "\n" + alignY)

A = "SALPQPTTPVSSFTSGSMLGRTDTALTNTYSAL" 
B = "PSPTMEAVTSVEASTASHPHSTSSYFATTYYHLY"
alignA, alignB = needlemanWunsch(A, B)
print("\n" + alignA + "\n" + alignB)
