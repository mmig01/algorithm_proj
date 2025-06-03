# fft_convolution.py

import numpy as np


if __name__ == "__main__":
    # Example usage:
    # A(x) = 1 + 2x + 3x^2  →  coefficients [1, 2, 3]
    # B(x) = 4 + 5x + 6x^2 + 7x^3  →  coefficients [4, 5, 6, 7]

   # 20 비트 소수
    p = 1048573
    n = 4
    w = 683314

    test_read = ["ATGCATG"]
    test_ref = "ATGCATGTTTTTTTTTTTTTTTTTTTTTTTTTT"

    inv_test_vec = []
    for i in range(len(test_read[0])):
        if test_read[0][i] == 'A':
            inv_test_vec.append(w**(n - 1) % p)
        elif test_read[0][i] == 'T':
            inv_test_vec.append(w**(n - 2) % p)
        elif test_read[0][i] == 'G':
            inv_test_vec.append(w**(n - 3) % p)
        elif test_read[0][i] == 'C':
            inv_test_vec.append(w**(n - 4) % p)
    inv_test_vec.reverse()
    
    test_vec = []
    for i in range(len(test_ref)):
        if test_ref[i] == 'A':
            test_vec.append(w)
        elif test_ref[i] == 'T':
            test_vec.append(w**2 % p)
        elif test_ref[i] == 'G':
            test_vec.append(w**3 % p)
        elif test_ref[i] == 'C':
            test_vec.append(w**4 % p)
            
    print("test_vec:", test_vec)
    print("inv_test_vec:", inv_test_vec)

    A = np.array(inv_test_vec, dtype=np.int64)
    B = np.array(test_vec, dtype=np.int64)

    N = len(A) + len(B)
    print("A:", A)
    print("B:", B)

    C = np.polymul(A, B)
    C = C % p  # Apply modulo operation
    print("Convolution result C:", C)

    # Expected output: [ 4 13 28 34 33 21 ]