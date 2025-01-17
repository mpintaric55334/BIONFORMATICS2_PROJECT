import os


def score(pair: tuple[str, str], match: int, mismatch: int, gap: int):
    n = len(pair[0])
    total_score = 0
    counter = 0
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
            if symbol1 == "-" and pair[1][i-1] == "-":
                counter +=1
    return total_score


def scores_from_folder(folder_path: str, match: int, mismatch: int, gap: int):
    scores = []
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            seq1 = lines[0].strip()
            seq2 = lines[1].strip()
            pair = (seq1, seq2)
            scores.append(score(pair, match, mismatch, gap))
    return scores

print(scores_from_folder("../data/alignments_hmm", 3, 2, -1))