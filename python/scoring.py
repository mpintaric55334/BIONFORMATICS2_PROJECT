import os
import numpy as np
import matplotlib.pyplot as plt


EDNA = {"A": {"A": 5, "C": -4, "G": -4, "T": -4, "-": -4},
        "C": {"A": -4, "C": 5, "G": -4, "T": -4, "-": -4},
        "T": {"A": -4, "C": -4, "G": -4, "T": 5, "-": -4},
        "G": {"A": -4, "C": -4, "G": 5, "T": -4, "-": -4},
        "-": {"A": -4, "C": -4, "G": -4, "T": -4, "-": -4}}


def score(pair: tuple[str, str], match: int, mismatch: int, gap: int):
    """
    Function that calculates normal score.

    Arguments:
        - pair : tuple[str, str] => pair of aligned sequences
        - match: int => score update on match
        - mismatch: int => score update on mismatch
        - gap: int => score update on gap

    Returns:
        - int: total score of alignment
    """
    n = len(pair[0])
    total_score = 0
    for i in range(n):
        symbol1 = pair[0][i]
        symbol2 = pair[1][i]
        if symbol1 != "-" and symbol2 != "-":
            if symbol1 == symbol2:
                total_score += match
            else:
                total_score += mismatch
        else:
            total_score += gap
    return total_score


def EDNAscore(pair: tuple[str, str]):
    """
    Function that calculates normal score.

    Arguments:
        - pair : tuple[str, str] => pair of aligned sequences

    Returns:
        - int: total score of alignment
    """
    n = len(pair[0])
    total_score = 0
    for i in range(n):
        symbol1 = pair[0][i]
        symbol2 = pair[1][i]
        total_score += EDNA[symbol1][symbol2]
    return total_score


def scores_from_folder(folder_path: str, match: int, mismatch: int, gap: int):
    """
    Function that return score arrays for every alignment folder.

    Arguments:
        - pair : tuple[str, str] => pair of aligned sequences
        - match: int => score update on match
        - mismatch: int => score update on mismatch
        - gap: int => score update on gap

    Returns:
        - scores: list[int]  => normal scores
        - edna_scores: list[int] => edna scores
    """
    scores = []
    edna_scores = []
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            seq1 = lines[0].strip()
            seq2 = lines[1].strip()
            pair = (seq1, seq2)
            scores.append(score(pair, match, mismatch, gap))
            edna_scores.append(EDNAscore(pair))
    scores = np.array(scores)
    edna_scores = np.array(edna_scores)
    return scores, edna_scores


list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20]
rand_scores, rand_edna_scores = scores_from_folder(
    "../data/alignments_hmm_random", 1, -1, -1)
estim_scores, estim_edna_scores = scores_from_folder(
    "../data/alignments_hmm_estim", 1, -1, -1)
rand_train_scores, rand_train_edna_scores = scores_from_folder(
    "../data/alignments_hmm_rand_train", 1, -1, -1)
estim_train_scores, estim_train_edna_scores = scores_from_folder(
    "../data/alignments_hmm_estim_train", 1, -1, -1)
needle_scores, needle_edna_scores = scores_from_folder(
    "../data/alignments_needle", 1, -1, -1)
print("Average normal score for random: ", np.mean(rand_scores))
print("Average EDNA score for random: ", np.mean(rand_edna_scores))

print("Average normal score for estim: ", np.mean(estim_scores))
print("Average EDNA score for estim: ", np.mean(estim_edna_scores))

print("Average normal score for random train: ", np.mean(rand_train_scores))
print("Average EDNA score for random train: ", np.mean(rand_train_edna_scores))

print("Average normal score for estim train: ", np.mean(estim_train_scores))
print("Average EDNA score for random: ", np.mean(estim_train_edna_scores))

print("Average normal score for needleman: ", np.mean(needle_scores))
print("Average EDNA score for needleman: ", np.mean(needle_edna_scores))

# normal scores
plt.scatter(list, rand_scores, label="Random", color='blue')
plt.scatter(list, estim_scores, label="Estimation", color='red')
plt.scatter(list, estim_train_scores, label="Trained", color='green')
plt.scatter(list, needle_scores, label="Needleman", color='black')
plt.legend()
plt.show()

#EDNA scores
plt.scatter(list, rand_edna_scores, label="Random", color='blue')
plt.scatter(list, estim_edna_scores, label="Estimation", color='red')
plt.scatter(list, estim_train_edna_scores, label="Trained", color='green')
plt.scatter(list, needle_edna_scores, label="Needleman", color='black')
plt.legend()
plt.show()