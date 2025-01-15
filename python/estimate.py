import os


def check_state(symbol1: str, symbol2: str):
    """
    Function that checks the current state
    and returns it.

    Arguments:
        - symbol1: str => first element of state
        - symbol2: str => second element of state

    Returns:
        - [empty]: str => current state
    """
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


def update_pi(pi: dict, pair: tuple[str, str]):
    """
    Function that updates pi values for every
    pair in estimation.

    Arguments:
        - pi: dict => dictionary of transmission probabilites
        - pair: tuple(str, str) => pair of sequences

    Returns:
        - pi: dict => updated dictionary of probabilities
    """
    first, second = pair
    state = check_state(first[0], second[0])
    pi[state] += 1
    return pi


def update_A(A: dict, pair: tuple[str, str]):
    """
    Function that updates A values for every
    pair in estimation.

    Arguments:
        - A: dict => dictionary of transmission probabilites
        - pair: tuple(str, str) => pair of sequences

    Returns:
        - A: dict => updated dictionary of probabilities
    """
    first, second = pair
    N = len(first)
    state = check_state(first[0], second[0])

    for i in range(1, N):
        new_state = check_state(first[i], second[i])
        A[state][new_state] += 1
        state = new_state
    return A


def update_E(E: dict, pair: tuple[str, str]):
    """
    Function that updates E values for every
    pair in estimation.

    Arguments:
        - E: dict => dictionary of transmission probabilites
        - pair: tuple(str, str) => pair of sequences

    Returns:
        - E: dict => updated dictionary of probabilities
    """
    first, second = pair
    N = len(first)
    for i in range(N):
        state = check_state(first[i], second[i])
        emission = first[i] + second[i]
        E[state][emission] += 1
    return E


def calculate_pi(pi: dict):
    """
    Function that calculates final pi probabilities.

    Arguments:
        - pi: dict => dict of pi probabilities

    Returns:
        - pi: dict => final dict of pi probabilities
    """
    n = sum(pi.values())
    for state in pi:
        pi[state] /= n
    return pi


def calculate_A(A: dict):
    """
    Function that calculates final A probabilities.

    Arguments:
        - A: dict => dict of A probabilities

    Returns:
        - A: dict => final dict of A probabilities
    """
    for state in A:
        n = sum(A[state].values())
        for transmission_state in A[state]:
            A[state][transmission_state] /= n
    return A


def calculate_E(E: dict):
    """
    Function that calculates final E probabilities.

    Arguments:
        - E: dict => dict of E probabilities

    Returns:
        - E: dict => final dict of E probabilities
    """
    for state in E:
        n = sum(E[state].values())
        for emission in E[state]:
            E[state][emission] /= n
    return E


def calculate_values(A: dict, E: dict, pi: dict, estimate_path: str):
    """
    Function that calculates estimate probabilities fro all estimate
    pairs.

    Arguments:
        - A: dict => dict of A probabilities
        - E: dict => dict of E probabilities
        - pi: dict => dict of pi probabilities
        - estimate_path: str => path for estimation pairs

    Returns:
        - A: dict => final dict of A probabilities
        - E: dict => final dict of E probabilities
        - pi: dict => final dict of pi probabilities
    """
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


def write_probs(A: dict, E: dict, pi: dict, transmission_path: str):
    """
    Function that writes estimate probabilities into storage paths.

    Arguments:
        - A: dict => dict of A probabilities
        - E: dict => dict of E probabilities
        - pi: dict => dict of pi probabilities
        - transmission_path: str => path for storage
        of estimated values
    """
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