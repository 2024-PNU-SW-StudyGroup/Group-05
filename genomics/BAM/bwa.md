# SAM/BAM 파일 형식 및 BWA-MEM2 가이드
작성일: 2024-11-07

## 목차
1. [소개](#소개)
2. [SAM/BAM 파일 형식](#samBam-파일-형식)
3. [BWA-MEM2](#bwa-mem2)
4. [분석 파이프라인](#분석-파이프라인)
5. [유용한 도구들](#유용한-도구들)
6. [결과 해석](#결과-해석)

## 소개
SAM(Sequence Alignment/Map)과 BAM(Binary Alignment/Map)은 시퀀싱 리드의 레퍼런스 지놈에 대한 정렬 결과를 저장하는 표준 파일 형식입니다. BWA-MEM2는 차세대 시퀀싱 데이터를 레퍼런스 지놈에 매핑하는 고성능 도구입니다.

## SAM/BAM 파일 형식

### 헤더 섹션
SAM 파일은 '@'로 시작하는 헤더 라인들로 시작합니다:
```
@HD     VN:1.6  SO:coordinate
@SQ     SN:chr1 LN:248956422
@RG     ID:S1   SM:sample1
@PG     ID:bwa  PN:bwa  VN:0.7.17
```

주요 헤더 태그:
- @HD: 헤더 라인
- @SQ: 레퍼런스 시퀀스 정보
- @RG: 리드 그룹
- @PG: 프로그램
- @CO: 주석

### 정렬 섹션
각 정렬은 탭으로 구분된 11개의 필수 필드를 포함합니다:

1. QNAME: Query 이름
2. FLAG: 비트 플래그
3. RNAME: 레퍼런스 이름
4. POS: 1-based 레퍼런스 위치
5. MAPQ: 매핑 품질
6. CIGAR: CIGAR 문자열
7. RNEXT: mate의 레퍼런스 이름
8. PNEXT: mate의 위치
9. TLEN: 템플릿 길이
10. SEQ: 시퀀스
11. QUAL: 품질 점수

예시:
```
read1   99  chr1    1000    60  100M    =   1100    200 ATCG... @@@@...
read2   147 chr1    1100    60  100M    =   1000    -200    GCTA... ####...
```

### FLAG 필드 해석
자주 사용되는 FLAG 값:
- 0x1 (1): paired-end
- 0x2 (2): properly paired
- 0x4 (4): unmapped
- 0x8 (8): mate unmapped
- 0x10 (16): reverse strand
- 0x20 (32): mate reverse strand
- 0x40 (64): first in pair
- 0x80 (128): second in pair
- 0x100 (256): secondary alignment
- 0x400 (1024): duplicate
- 0x800 (2048): supplementary alignment

### CIGAR 문자열
정렬 상태를 나타내는 문자열:
- M: 매치/미스매치
- I: 삽입
- D: 삭제
- N: 스킵
- S: 소프트 클리핑
- H: 하드 클리핑
- P: 패딩
- =: 시퀀스 매치
- X: 시퀀스 미스매치

## BWA-MEM2

### 특징
- 빠른 속도와 높은 정확도
- SIMD 가속화 지원
- 긴 리드 지원
- 키메릭 리드 처리
- 멀티스레딩 지원

### 설치 방법
```bash
# Github에서 설치
git clone https://github.com/bwa-mem2/bwa-mem2.git
cd bwa-mem2
make

# Conda를 이용한 설치
conda install -c bioconda bwa-mem2
```

### 인덱스 생성
```bash
bwa-mem2 index reference.fa
```

### 매핑 실행
```bash
# Single-end 리드
bwa-mem2 mem reference.fa reads.fq > aln.sam

# Paired-end 리드
bwa-mem2 mem reference.fa read1.fq read2.fq > aln.sam

# 멀티스레딩 사용
bwa-mem2 mem -t 16 reference.fa reads.fq > aln.sam
```

### 주요 매개변수
- -t: 스레드 수
- -k: 시드 길이 (default: 19)
- -w: 밴드 폭 (default: 100)
- -d: off-diagonal X-dropoff (default: 100)
- -r: 시드 간격 (default: 1.5)
- -A: 매칭 점수 (default: 1)
- -B: 미스매치 페널티 (default: 4)
- -O: gap open 페널티 (default: 6)
- -E: gap extension 페널티 (default: 1)

## 분석 파이프라인

### 기본 워크플로우
1. 레퍼런스 인덱스 생성
2. BWA-MEM2로 매핑
3. SAM to BAM 변환
4. BAM 정렬
5. 중복 마킹
6. 인덱스 생성

## 유용한 도구들

### SAMtools
```bash
# SAM to BAM 변환
samtools view -bS aln.sam > aln.bam

# BAM 정렬
samtools sort aln.bam -o aln.sorted.bam

# 인덱스 생성
samtools index aln.sorted.bam

# 통계 확인
samtools flagstat aln.sorted.bam
```

### Picard
```bash
# 중복 마킹
java -jar picard.jar MarkDuplicates \
    I=input.bam \
    O=marked_duplicates.bam \
    M=marked_dup_metrics.txt
```

## 결과 해석

### 매핑 품질 평가
- MAPQ 점수 분포 확인
- 매핑율 계산
- 중복율 확인
- 커버리지 분석

### 주요 품질 지표
1. 매핑율: 총 리드 중 매핑된 리드의 비율
2. 고품질 매핑: MAPQ ≥ 30인 리드의 비율
3. 중복율: 중복 마킹된 리드의 비율
4. 커버리지: 각 위치별 평균 리드 깊이

### 문제 해결
1. 낮은 매핑율
   - 리드 품질 확인
   - 어댑터 오염 검사
   - 레퍼런스 적합성 확인

2. 높은 중복율
   - 라이브러리 복잡도 확인
   - PCR 사이클 수 조정
   - 입력 DNA 양 확인

3. 불균일한 커버리지
   - GC 편향 검사
   - 캡처 효율성 확인
   - 시퀀싱 깊이 조정