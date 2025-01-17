from minineedle import needle, core
import os


def generate_needleman_alignments(test_data_path: str, alignment_path: str):
    i = 0
    for file_name in os.listdir(test_data_path):
        file_path = os.path.join(test_data_path, file_name)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            seq1 = lines[0].strip()
            seq2 = lines[1].strip()
            alignment = needle.NeedlemanWunsch(seq1, seq2)
            alignment.align()
            al1, al2 = alignment.get_aligned_sequences(core.AlignmentFormat.str)
            align_file_path = os.path.join(alignment_path, "pair" + str(i))
            with open(align_file_path, "w") as align_file:
                align_file.write(f"{al1}\n{al2}")
        i += 1

generate_needleman_alignments("../data/data_test",
                              "../data/alignments_needle")
