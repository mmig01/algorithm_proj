import os
import time
import random
from method import brute_force, generate_reads, get_reference_list, make_sequence_list, packing

# 스크립트가 위치한 경로
script_dir = os.path.dirname(os.path.realpath(__file__))

# 그 디렉터리를 기준으로 상대경로 지정
original_ref_path = os.path.join(script_dir, "parameters/Original.txt")
reference_list_csv_path  = os.path.join(script_dir, "parameters/sequence_list.csv")
reads_list_csv_path = os.path.join(script_dir, "parameters/reads_list.csv")

if __name__ == "__main__":
      
    """
    1. make_sequence_list 함수를 사용하여 염기서열을 CSV 파일로 저장합니다.
    """
    make_sequence_list(original_ref_path, reference_list_csv_path)

    """
    2. get_sequence_list 함수를 사용하여 CSV 파일에서 염기서열을 읽어옵니다.
    """
    reference_list = get_reference_list(reference_list_csv_path)
    

    # 예시: L과 M을 입력 받아 read 생성 및 CSV 저장
    L = 200  # 예시 값, 실제 실행 시 바꿔 입력 가능합니다.
    M = 500  # 예시 값, 실제 실행 시 바꿔 입력 가능합니다.
    reads = generate_reads(file_path=original_ref_path, csv_path=reads_list_csv_path, L=L, M=M)
    # read 순서를 무작위로 배치
    random.shuffle(reads)
   
    """
    3. 
        1) brute-force 방식으로 염기서열 복구
        2) FFT를 사용하여 다항식 곱 계산 후 염기서열 복구
    """
    user_input = input("brute-force 방식으로 염기서열 복구를 원하시면 1, FFT를 사용하여 다항식 곱 계산 후 염기서열 복구를 원하시면 2 를 입력하세요: ")
    
    try:
        num = int(user_input)
        if num == 1:
            print("brute-force 방식으로 염기서열 복구를 시작합니다...")
            start_time = time.perf_counter()
            result_vec = brute_force(reference_list=reference_list, reads=reads, M=M, L=L)
            end_time = time.perf_counter()
            elapsed_ms = (end_time - start_time)
            print(f"실행 시간: {elapsed_ms:.2f} s")
        elif num == 2:
            print("FFT를 사용하여 다항식 곱 계산 후 염기서열 복구를 시작합니다...")
            start_time = time.perf_counter()
            result_vec = packing(reference_list=reference_list, reads=reads, M=M)
            end_time = time.perf_counter()
            elapsed_ms = (end_time - start_time)
            print(f"실행 시간: {elapsed_ms:.2f} s")
        else:
            print("잘못된 입력입니다. 1 또는 2를 입력해주세요.")
            exit(1)
    except ValueError:
        print("유효한 정수가 아닙니다.")
   
    """
    4. 특정 인덱스만 출력
    index : 30000 30036 30046 30877 30923 30976 31015 31273 31276 31294
    정상유전자 : G     C     T     T     G     C     C     G     G     C

    위 인덱스에 해당하는 유전자가 정상 유전자가 아니라면 알츠하이머 발병 가능성이 높음.
    """
    target_index = [30000, 30036, 30046, 30877, 30923, 30976, 31015, 31273, 31276, 31294]
    for idx in target_index:
        print(f"Index: {idx}, 결과: {result_vec[idx]}, 정상 유전자: {reference_list[idx]}")
        if result_vec[idx] != reference_list[idx]:
            print(f"알츠하이머 발병 가능성이 높습니다.")
            break

   