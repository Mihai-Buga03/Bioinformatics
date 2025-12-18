import numpy as np
import matplotlib.pyplot as plt

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=0):
    n = len(seq1)
    m = len(seq2)

    score_matrix = np.zeros((n + 1, m + 1))

    for i in range(n + 1):
        score_matrix[i][0] = i * gap
    for j in range(m + 1):
        score_matrix[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch
            
            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap
            
            score_matrix[i][j] = max(diagonal_score, up_score, left_score)

    align1 = ""
    align2 = ""
    visual_match = ""
    
    i, j = n, m

    path_x = [j]
    path_y = [i]

    matches = 0
    
    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        
        if seq1[i - 1] == seq2[j - 1]:
            diagonal_score = score_matrix[i - 1][j - 1] + match
        else:
            diagonal_score = score_matrix[i - 1][j - 1] + mismatch
            
        up_score = score_matrix[i - 1][j] + gap
        left_score = score_matrix[i][j - 1] + gap

        if current_score == diagonal_score:
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            if seq1[i-1] == seq2[j-1]:
                visual_match += "|"
                matches += 1
            else:
                visual_match += " "
            i -= 1
            j -= 1
        elif current_score == up_score:
            align1 += seq1[i - 1]
            align2 += "-"
            visual_match += " "
            i -= 1
        else:
            align1 += "-"
            align2 += seq2[j - 1]
            visual_match += " "
            j -= 1
            
        path_x.append(j)
        path_y.append(i)

    while i > 0:
        align1 += seq1[i - 1]
        align2 += "-"
        visual_match += " "
        path_x.append(j)
        path_y.append(i)
        i -= 1
    while j > 0:
        align1 += "-"
        align2 += seq2[j - 1]
        visual_match += " "
        path_x.append(j)
        path_y.append(i)
        j -= 1
        
    path_x.append(0)
    path_y.append(0)

    align1 = align1[::-1]
    align2 = align2[::-1]
    visual_match = visual_match[::-1]

    return align1, align2, visual_match, score_matrix, path_x, path_y, matches

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"
GAP = 0
MATCH = 1
MISMATCH = -1

a1, a2, bars, matrix, px, py, match_count = needleman_wunsch(S1, S2, MATCH, MISMATCH, GAP)

print(f"Sequence 1: {S1}")
print(f"Sequence 2: {S2}")
print("-" * 30)
print("Alignment Result:")
print(a1)
print(bars)
print(a2)
print("-" * 30)
print(f"Matches: {match_count}")
print(f"Length: {len(a1)}")
print(f"Similarity: {(match_count/len(a1))*100:.2f}%")

fig, ax = plt.subplots(1, 2, figsize=(12, 6))

cax = ax[0].imshow(matrix, cmap='magma', interpolation='nearest')
ax[0].set_title('Alignment Matrix Heatmap')
ax[0].set_xlabel('Sequence 2')
ax[0].set_ylabel('Sequence 1')
fig.colorbar(cax, ax=ax[0])

ax[1].imshow(np.zeros_like(matrix), cmap='Pastel1') 
ax[1].set_xticks(np.arange(-.5, len(S2)+1, 1), minor=True)
ax[1].set_yticks(np.arange(-.5, len(S1)+1, 1), minor=True)
ax[1].grid(which='minor', color='black', linestyle='-', linewidth=1)
ax[1].set_xticks([])
ax[1].set_yticks([])

ax[1].plot(px, py, color='red', linewidth=3, marker='s', markersize=5)

ax[1].invert_yaxis() 
ax[1].set_title('Traceback Path')

plt.tight_layout()
plt.show()