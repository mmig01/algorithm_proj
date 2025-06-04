# fft_convolution.py

import numpy as np
import csv

file_path = "/Users/mingikim/algorithm_proj/Original.txt"
csv_path = "/Users/mingikim/algorithm_proj/sequence_list.csv"

def make_sequence_list(file_path, csv_path):
    """
    주어진 텍스트 파일에서 염기서열을 읽어와 CSV 파일로 저장합니다.
    
    :param file_path: 원본 텍스트 파일 경로
    :param csv_path: 저장할 CSV 파일 경로
    """
    with open(file_path, "r") as f:
        seq_str = f.read().replace("\n", "").strip()

    seq_list = list(seq_str)  # 예: ["C", "A", "G", … ]

    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        for base in seq_list:
            writer.writerow([base])
    print(f"sequence_list.csv 파일이 생성되었습니다: {csv_path}")


def read_sequence_list(csv_path):
    """
    CSV 파일에서 염기서열 리스트를 읽어옵니다.
    
    :param csv_path: CSV 파일 경로
    :return: 염기서열 리스트
    """
    seq_list = []
    with open(csv_path, "r", newline="") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            seq_list.append(row[0])
    return seq_list


if __name__ == "__main__":
    # Example usage:
    # A(x) = 1 + 2x + 3x^2  →  coefficients [1, 2, 3]
    # B(x) = 4 + 5x + 6x^2 + 7x^3  →  coefficients [4, 5, 6, 7]
      
    """
    1. make_sequence_list 함수를 사용하여 염기서열을 CSV 파일로 저장합니다.
    """
    make_sequence_list(file_path, csv_path)

    """
    2. read_sequence_list 함수를 사용하여 CSV 파일에서 염기서열을 읽어옵니다.
    """
    seq_list_loaded = read_sequence_list(csv_path)
    print("30001번째 염기 : ", seq_list_loaded[30000])  # 30000번째 염기 출력

 

   # 20 비트 소수
    p = 1048573
    n = 4
    w = 683314

    test_read = ["ATGCATG"]
    test_ref = seq_list_loaded

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
            
    # print("test_vec:", test_vec)
    # print("inv_test_vec:", inv_test_vec)

    A = np.array(inv_test_vec, dtype=np.int64)
    B = np.array(test_vec, dtype=np.int64)
    
    N = len(A) + len(B)
    print("A:", A)
    print("B:", B)

    C = np.polymul(A, B)
    C = C % p  # Apply modulo operation
    for i in range(len(C)):
        if C[i] == len(A):
            print(f"염기서열 {i - len(A) + 2} 에서 발견")
    

   