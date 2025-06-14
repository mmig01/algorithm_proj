# Algorithm project
알고리즘 프로젝트 - 2020111995 김민기

# 실행법

### 가상환경 설치

1. 맥북

```py
# 가상환경 설치
python3 -m venv venv

# 가상환경 실행
source venv/bin/activate
```

2. 윈도우

```py
# 가상환경 설치 
python -m venv venv

# 가상환경 실행
source venv/Scripts/activate
```



### 필요 프로그램 설치

```py
pip install -r requirements.txt
```

### 프로그램 실행

```python
python -u "string_match.py 의 경로"

# ex) python -u "/Users/mingikim/algorithm_proj/string_match.py"
```



# 파일구조

1. parameters 폴더
   - Original.txt : 50000 개의 염기서열이 저장 되어 있는 파일
   - sequence_list.csv : Original.txt 내부에 있는 레퍼런스를 list 형태로 변환 후 저장하는 파일
   - reads_list.csv : 패턴 매칭 알고리즘 실행 시, reference 로부터 생성된 read 가 저장되는 파일
   
2. projcet 폴더
   - method.py : 알고리즘 실행 시 필요한 함수들을 모아 놓은 파이썬 파일

3. string_match.py : main 함수가 위치하는 파일. 해당 파일을 실행하면 대화형 프롬프트가 실행 되고, 염기 서열 복구 및 알츠하이머 판단 결과를 보여줌.
