import os


def check_state(symbol1, symbol2):
    if symbol1 != "-" and symbol2 != "-":
        return "MM"
    if symbol1 != "-" and symbol2 == "-":
        return "Ix"
    if symbol1 == "-" and symbol2 != "-":
        return "Iy"


def initialize_A():
    """
    Function that initializes the transmission
    matrix A (state to state transmission)
    """
    A = {
        "MM": {
            "MM": 0,
            "Ix": 0,
            "Iy": 0
        },
        "Ix": {
            "MM": 0,
            "Ix": 0,
            "Iy": 0
        },
        "Iy": {
            "MM": 0,
            "Ix": 0,
            "Iy": 0
        }
    }
    return A


def initialize_pi():
    """
    Function that initializes initial probability
    array.
    """
    pi = {
        "MM": 0,
        "Ix": 0,
        "Iy": 0
    }
    return pi


def initialize_E():
    """
    Function that initializes initial emission probabilities.
    """
    E = {
        "MM": {
            "AA": 0,
            "AC": 0,
            "AT": 0,
            "AG": 0,
            "CA": 0,
            "CC": 0,
            "CT": 0,
            "CG": 0,
            "TA": 0,
            "TC": 0,
            "TT": 0,
            "TG": 0,
            "GA": 0,
            "GC": 0,
            "GT": 0,
            "GG": 0,
        },
        "Ix": {
            "A-": 0,
            "C-": 0,
            "T-": 0,
            "G-": 0
        },
        "Iy": {
            "-A": 0,
            "-C": 0,
            "-T": 0,
            "-G": 0
        }

    }
    return E


def update_pi(pi, pair):
    first, second = pair
    state = check_state(first[0], second[0])
    pi[state] += 1
    return pi


def update_A(A, pair):
    first, second = pair
    N = len(first)
    state = check_state(first[0], second[0])

    for i in range(1, N):
        new_state = check_state(first[i], second[i])
        A[state][new_state] += 1
        state = new_state
    return A


def update_E(E, pair):
    first, second = pair
    N = len(first)
    for i in range(N):
        state = check_state(first[i], second[i])
        emission = first[i] + second[i]
        E[state][emission] += 1
    return E


def calculate_pi(pi):
    n = sum(pi.values())
    for state in pi:
        pi[state] /= n
    return pi


def calculate_A(A):
    for state in A:
        n = sum(A[state].values())
        for transmission_state in A[state]:
            A[state][transmission_state] /= n
    return A


def calculate_E(E):
    for state in E:
        n = sum(E[state].values())
        for emission in E[state]:
            E[state][emission] /= n
    return E


def calculate_values(A, E, pi, estimate_path):
    for file_name in os.listdir(estimate_path):
        file_path = os.path.join(estimate_path, file_name)

        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                lines = file.readlines()
                seq1 = lines[0].strip()
                seq2 = lines[1].strip()
                pair = (seq1, seq2)
                A = update_A(A, pair)
                E = update_E(E, pair)
                pi = update_pi(pi, pair)
    A = calculate_A(A)
    E = calculate_E(E)
    pi = calculate_pi(pi)
    return A, E, pi


def write_probs(A, E, pi, transmission_path):
    os.makedirs(transmission_path, exist_ok=True)
    # A
    A_path = os.path.join(transmission_path, "A.txt")
    with open(A_path, "w") as file:
        to_write = ""
        for state in A:
            for new_state in A[state]:
                to_write += str(A[state][new_state]) + "\n"
        to_write = to_write.strip()
        file.write(f"{to_write}")
    # E
    E_path = os.path.join(transmission_path, "E.txt")
    with open(E_path, "w") as file:
        to_write = ""
        for state in E:
            for new_state in E[state]:

                to_write += str(E[state][new_state]) + "\n"
        to_write = to_write.strip()
        file.write(f"{to_write}")
    # pi
    pi_path = os.path.join(transmission_path, "pi.txt")
    with open(pi_path, "w") as file:
        to_write = ""
        for state in pi:
            to_write += str(pi[state]) + "\n"
        to_write = to_write.strip()
        file.write(f"{to_write}")


A = initialize_A()
E = initialize_E()
pi = initialize_pi()
A, E, pi = calculate_values(A, E, pi, "C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\data\\train_data\\estimate")
write_probs(A, E, pi, "C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\transmission_values\\estimate")