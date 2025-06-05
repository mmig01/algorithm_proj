import csv

import numpy as np

# --- NumPy FFT 기반 다항식 곱 ---
def conv_fft(a, b):
    """
    NumPy FFT를 이용한 실수 다항식 컨볼루션.
    * a, b : 정수 리스트 (또는 배열)
    * 반환  : 정수 리스트 (길이 len(a)+len(b)-1)
    """
    # 필요 길이
    need = len(a) + len(b) - 1
    n = 1 << (need - 1).bit_length()  # 2^k >= need

    # float64 배열 패딩
    a_pad = np.zeros(n, dtype=np.float64)
    b_pad = np.zeros(n, dtype=np.float64)
    a_pad[: len(a)] = a
    b_pad[: len(b)] = b

    # FFT → 곱 → IFFT
    fa = np.fft.fft(a_pad)
    fb = np.fft.fft(b_pad)
    fc = fa * fb
    c_pad = np.fft.ifft(fc)

    # 반올림 후 int64 캐스팅
    c_int = np.rint(c_pad.real).astype(np.int64)
    return c_int[:need].tolist()


"""
메소드를 모아 놓은 파일
"""

def brute_force(reference_list, reads, M, L):
    """
    brute-force 방식으로 염기서열 복구하는 함수

    reference_list: 레퍼런스 염기 서열 리스트
    reads: 읽어온 read 염기 서열 리스트
    M: read 개수
    L: 각 read의 길이
    """
    
    # 레퍼런스 길이의 0 으로 초기화된 배열 생성
    result_vec = ["0" for _ in range(len(reference_list))]
    
    # read 개수 M 만큼의 brute-force 연산 수행
    for i in range(M):
        read = reads[i]
        for j in range(len(reference_list) - L + 1):
            is_match = True
            for k in range(L):
                if reference_list[j + k] != read[k]:
                    is_match = False
            if is_match:
                # 일치하는 부분이 발견되면 해당 위치에 read를 삽입
                for k in range(L):
                    result_vec[j + k] = read[k]
                break
    return result_vec

def packing(reference_list, reads, M):
    """
    root of unity 를 이용하여 각 염기를 다항식 계수로 packing 하는 함수

    1. 레퍼런스 염기 서열
    - A : w, T: w^2, G: w^3, C: w^4
    2. 읽어온 read 염기 서열
    - A : w^(n-1), T: w^(n-2), G: w^(n-3), C: w^(n-4)
    """

    # 20 비트 소수로 설정
    p = 1048573
    n = 4
    w = 683314


    packed_reference = []
    for i in range(len(reference_list)):
        if reference_list[i] == 'A':
            packed_reference.append(w)
        elif reference_list[i] == 'T':
            packed_reference.append(w**2 % p)
        elif reference_list[i] == 'G':
            packed_reference.append(w**3 % p)
        elif reference_list[i] == 'C':
            packed_reference.append(w**4 % p)

    # 레퍼런스 길이의 0 으로 초기화된 배열 생성
    result_vec = ["0" for _ in range(len(reference_list))]
    
    # read 개수 M 만큼의 fft 연산 수행
    for i in range(M):
        read_vec_represented_by_root_of_unity = []
        for j in range(len(reads[0])):
            if reads[i][j] == 'A':
                read_vec_represented_by_root_of_unity.append(w**(n - 1) % p)
            elif reads[i][j] == 'T':
                read_vec_represented_by_root_of_unity.append(w**(n - 2) % p)
            elif reads[i][j] == 'G':
                read_vec_represented_by_root_of_unity.append(w**(n - 3) % p)
            elif reads[i][j] == 'C':
                read_vec_represented_by_root_of_unity.append(w**(n - 4) % p)
        
        # read 다항식의 경우 역순으로 뒤집음
        inv_read_vec = read_vec_represented_by_root_of_unity[::-1]  # 역순으로 뒤집기

        # NumPy FFT 기반 다항식 곱
        C = conv_fft(packed_reference, inv_read_vec)
        C = [coef % p for coef in C]  # 모듈로 연산
        read_length = len(inv_read_vec)
        for k in range(len(C)):
            if C[k] == read_length:
                target_index = k - read_length + 1
                result_vec[target_index : target_index + read_length] = reads[i]
        
    return result_vec

def make_sequence_list(file_path, csv_path):
    """
    주어진 reference 텍스트 파일에서 염기서열을 읽어와 CSV 파일로 저장하는 함수
    
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


def get_reference_list(csv_path):
    """
    CSV 파일에서 염기서열 리스트를 읽어오는 메소드
    
    :param csv_path: CSV 파일 경로
    :return: 염기서열 리스트
    """
    seq_list = []
    with open(csv_path, "r", newline="") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            seq_list.append(row[0])
    return seq_list

def generate_reads(file_path, csv_path, L, M):
    """
    reference 를 균등하게 잘라서 read 를 생성하는 함수

    file_path: Original.txt 파일 경로
    csv_path: 생성된 read를 저장할 CSV 파일 경로
    read_length: 각 read의 길이 L
    num_reads: 생성할 read 개수 M
    returns: list of reads (문자열)
    """
    # 1. Original.txt에서 FASTA 형식 없이 단일 시퀀스로 읽어옴
    with open(file_path, "r") as f:
        seq = "".join([line.strip() for line in f if not line.startswith(">")])
    seq_len = len(seq)
    if M < 2:
        raise ValueError("M 은 최소 2 이상이어야 합니다.")
    if L > seq_len:
        raise ValueError("read_length가 시퀀스 길이보다 클 수 없습니다.")
    
    # 2. M개의 시작 위치를 균등 분포로 계산 (첫 위치=0, 마지막 위치=seq_len - read_length)
    positions = [
        int(round(i * (seq_len - L) / (M - 1)))
        for i in range(M)
    ]

    # 3. 각 위치에서 read_length 길이만큼 잘라서 read 생성
    reads = [seq[p : p + L] for p in positions]

    
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["read_id", "sequence"])
        for idx, r in enumerate(reads, start=1):
            writer.writerow([idx, r])
    print(f"generated_reads.csv 파일이 생성되었습니다: {csv_path}")


    return reads
