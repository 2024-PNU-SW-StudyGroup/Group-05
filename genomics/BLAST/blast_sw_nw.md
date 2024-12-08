# 서열 비교 알고리즘 및 BLAST 가이드
작성일: 2024-11-07

## 목차
1. [소개](#소개)
2. [기본 알고리즘](#기본-알고리즘)
   - [Needleman-Wunsch](#needleman-wunsch)
   - [Smith-Waterman](#smith-waterman)
   - [FM-Index](#fm-index)
3. [BLAST (Basic Local Alignment Search Tool)](#blast)
4. [실습 가이드](#실습-가이드)
5. [결과 해석](#결과-해석)

## 소개
서열 비교(Sequence alignment)는 생물정보학의 기본이 되는 분석 방법입니다. 본 문서에서는 주요 서열 비교 알고리즘과 BLAST의 사용법을 다룹니다.

## 기본 알고리즘

### Needleman-Wunsch
전역 정렬(Global alignment) 알고리즘으로, 두 서열의 전체를 비교합니다.

#### 특징
- 동적 프로그래밍 기반
- 최적의 전역 정렬 보장
- 시간 복잡도: O(mn)
- 공간 복잡도: O(mn)

#### 점수 행렬
```
match score: +1
mismatch penalty: -1
gap penalty: -2

      A  C  G  T
   0 -2 -4 -6 -8
A -2  1 -1 -3 -5
C -4 -1  2  0 -2
G -6 -3  0  3  1
T -8 -5 -2  1  4
```

### Smith-Waterman
국소 정렬(Local alignment) 알고리즘으로, 서열의 일부분을 비교합니다.

#### 특징
- 동적 프로그래밍 기반
- 최적의 국소 정렬 보장
- 음수 점수를 0으로 처리
- 시간 복잡도: O(mn)

#### 주요 차이점 (vs Needleman-Wunsch)
1. 행렬 초기화: 첫 행과 열을 0으로 초기화
2. 점수 계산: 음수 점수를 0으로 대체
3. 트레이스백: 최대 점수에서 시작

### FM-Index
접미사 배열(Suffix Array)과 BWT(Burrows-Wheeler Transform)를 기반으로 한 텍스트 검색 자료구조입니다.

#### 특징
- 빠른 검색 속도
- 메모리 효율적
- 전체 텍스트 압축 지원
- 정확한 매칭 검색에 최적화

#### BWT 예시
```
Original: BANANA$
Rotations:
BANANA$
A$BANAN
ANA$BAN
ANANA$B
BANANA$
NA$BANA
NANA$BA

Sorted:
A$BANAN
ANA$BAN
ANANA$B
BANANA$
NA$BANA
NANA$BA

BWT: BNN$AAA
```

## BLAST
BLAST는 서열 데이터베이스 검색을 위한 휴리스틱 알고리즘입니다.

### BLAST 종류
- blastn: 뉴클레오티드 vs 뉴클레오티드
- blastp: 단백질 vs 단백질
- blastx: 번역된 뉴클레오티드 vs 단백질
- tblastn: 단백질 vs 번역된 뉴클레오티드
- tblastx: 번역된 뉴클레오티드 vs 번역된 뉴클레오티드

### 알고리즘 단계
1. Seeding: 짧은 일치 서열 찾기
2. Extension: seed 주변 확장
3. Evaluation: 통계적 유의성 평가

### 주요 매개변수
- -query: 질의 서열 파일
- -db: 검색할 데이터베이스
- -evalue: 기대값 임계치
- -outfmt: 출력 형식
- -num_threads: 사용할 스레드 수

## 실습 가이드

### BLAST 설치
```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# Conda
conda install -c bioconda blast
```

### 데이터베이스 생성
```bash
makeblastdb -in sequences.fasta -dbtype nucl -out mydb
```

### BLAST 실행
```bash
# 기본 검색
blastn -query query.fa -db mydb -out results.txt

# 표 형식 출력
blastn -query query.fa -db mydb -outfmt 6 -out results.tsv

# XML 형식 출력
blastn -query query.fa -db mydb -outfmt 5 -out results.xml
```

## 결과 해석

### BLAST 결과 필드 (표 형식)
1. query id
2. subject id
3. % identity
4. alignment length
5. mismatches
6. gap opens
7. query start
8. query end
9. subject start
10. subject end
11. evalue
12. bit score

### E-value 해석
- E-value < 1e-50: 매우 강한 일치
- E-value < 1e-10: 강한 일치
- E-value < 0.01: 유의미한 일치
- E-value > 0.01: 약한 일치

### Bit Score 해석
- 높을수록 더 좋은 정렬
- 데이터베이스 크기와 독립적
- 일반적으로 50 이상이 유의미