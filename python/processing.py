from Bio import SeqIO
import random
import os


def generate_data_indices(n, n_training_data: int, n_testing_data: int,
                          n_estimating_data: int):
    """
    Helper function that generates random indices in desired range.

    Arguments:
        - n: int => total number of datapoints
        - n_training_data => number of training datapoints
        - n_testing_data => number of testing datapoints
        - n_estimating_data => number of estimating_data

    Returns:
        - training_indices: list[int] => list of datapoint indices
        for training
        - testing_indices: list[int] => list of datapoint indices
        for testing
        - estimating_indices: list[int] => list of datapoint indices
        for estimating
    """

    all_values = list(range(n))
    random.shuffle(all_values)

    training_indices = all_values[:n_training_data]
    testing_indices = all_values[n_training_data:n_training_data
                                 + n_testing_data]
    estimating_indices = all_values[n_training_data +
                                    n_testing_data:n_training_data +
                                    n_testing_data + n_estimating_data]

    return training_indices, testing_indices, estimating_indices


def load_data(fasta_file: str, n_training_data: int, n_testing_data: int,
              n_estimating_data: int):
    """
    Function that loads a predetermined number of  random
    sequences for training, testing, and estimation.

    Arguments:
        - fasta_file: str => fasta file from which the data is loaded
        - n_training_data: int => number of training sequences
        - n_testing_data: int => number of testing sequences
        - n_estimating_data: int => number of sequences for estimation

    Returns:
        - training_seq: list[str] => list of training sequences
        - testing_seq: list[str] => list of testing sequences
        - estimating_seq: list[str] => list of estimating sequences

    """
    # discard sequences which have other nucleotides(not acgt)
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        set_seq = set(seq)
        allowed = True
        for s in set_seq:
            if s not in "ACGT-":
                allowed = False
        if not allowed:
            continue
        sequences.append(seq)
    i = 0
    training_indices, testing_indices, estimating_indices = generate_data_indices(
        len(sequences), n_training_data, n_testing_data, n_estimating_data)
    training_seq = []
    testing_seq = []
    estimating_seq = []
    for seq in sequences:
        if i in training_indices:
            training_seq.append(seq)
        elif i in testing_indices:
            testing_seq.append(seq)
        elif i in estimating_indices:
            estimating_seq.append(seq)
        i += 1
    return training_seq, testing_seq, estimating_seq


def make_pairs(sequences: list[str]):
    """
    Helper function that makes pairs of sequences
    from list of sequences.

    Arguments:
        - sequences: list[str] => list of sequences

    Returns:
        - seq_pairs: list[tuple(str, str)] => list
        of pair of sequences
    """
    seq_pairs = []
    N = len(sequences)
    for i in range(N):
        for j in range(i + 1, N):
            seq_pairs.append((sequences[i], sequences[j]))
    return seq_pairs


def MSA_TO_PSA(seq_pairs: list[str]):
    """
    Function that prepares sequences from multiple sequence alignment
    to pairwise sequence alignment. It removes unnecessary dashes.
    Ex:
        C-AGG
        A-AGG
    Becomes:
        CAGG
        AAGG

    Arguments:
        - seq_pairs: list[tuple(str, str)] => list of sequence pairs

    Returns:
        - new_pairs: list[tuple(str, str)] => list of sequence pairs
    """
    new_pairs = []
    for pair in seq_pairs:
        first, second = pair
        first_new, second_new = "", ""
        N = len(first)
        for i in range(N):
            if first[i] == "-" and second[i] == "-":
                continue
            first_new += first[i]
            second_new += second[i]
        assert len(first_new) == len(second_new), (
            "Sequences must be of same length after dash removal"
        )
        new_pairs.append((first_new, second_new))
    return new_pairs


def remove_all_dashes(seq_pairs: list[str]):
    """
    Function that removes all dashes from list
    of sequence pairs. This is mainly for testing data.

    Arguments:
        - seq_pairs: list[str] => list of sequence pairs

    Returns:
        - new_pairs: list[str] => list of sequence pairs
    """

    new_pairs = []
    for pair in seq_pairs:
        first, second = pair
        first = first.replace("-", "")
        second = second.replace("-", "")
        new_pairs.append((first, second))
    return new_pairs


def write_data(folder_path, seq_pairs):
    """
    Function that writes each pair of strings in the specified folder.

    Args:
        folder_path: str => The directory of pair files.
        seq_pairs: list[tuple(str, str)]: List of seq pairs.

    """
    os.makedirs(folder_path, exist_ok=True)
    N = len(seq_pairs)
    for i in range(N):
        file_path = os.path.join(folder_path, "pair" + str(i))
        with open(file_path, "w") as file:
            file.write(f"{seq_pairs[i][0]}\n{seq_pairs[i][1]}")


# load data
train_data, test_data, estimate_data = load_data("C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\data\\HIV1_ALL_2022_genome_DNA.fasta", 3, 3, 3)
# make pairs
train_pairs = make_pairs(train_data)
test_pairs = make_pairs(test_data)
estimate_pairs = make_pairs(test_data)
# MSA to PSA
train_pairs = MSA_TO_PSA(train_pairs)
test_pairs = MSA_TO_PSA(test_pairs)
estimate_pairs = MSA_TO_PSA(estimate_pairs)
# remove dashes for test
test_pairs = remove_all_dashes(test_pairs)
# write
write_data("C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\data\\train_data\\baum_welch_train", train_pairs)
write_data("C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\data\\train_data\\estimate", estimate_pairs)
write_data("C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\data\\test_data", test_pairs)
