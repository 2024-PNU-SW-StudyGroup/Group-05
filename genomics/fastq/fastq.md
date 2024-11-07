# FastQ 품질 관리 및 전처리 가이드
작성일: 2024-11-07

## 목차
1. [소개](#소개)
2. [FastQ 파일 형식](#fastq-파일-형식)
3. [FastQC](#fastqc)
4. [Cutadapt](#cutadapt)
5. [분석 파이프라인](#분석-파이프라인)
6. [결과 해석](#결과-해석)

## 소개
본 문서는 NGS(Next Generation Sequencing) 데이터의 품질 관리와 전처리 과정을 설명합니다. FastQ 파일의 품질을 확인하는 FastQC와 어댑터 및 저품질 시퀀스를 제거하는 Cutadapt의 사용법을 다룹니다.

## FastQ 파일 형식
FastQ는 DNA/RNA 시퀀싱 데이터를 저장하는 텍스트 기반 형식입니다.

### 구조
각 리드는 4줄로 구성됩니다:
```
@SEQ_ID                    # 시퀀스 식별자 (@ 로 시작)
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT    # 염기서열
+                         # 구분자 (+ 로 시작)
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65    # 품질 점수
```

### 품질 점수 (Phred Score)
- ASCII 문자로 인코딩된 품질 점수
- Phred+33 또는 Phred+64 인코딩 사용
- Q = -10 log₁₀ P (P: 에러 확률)
  - Q30: 99.9% 정확도
  - Q20: 99% 정확도

## FastQC
FastQC는 시퀀싱 데이터의 품질을 평가하는 도구입니다.

### 설치 방법
```bash
# Conda를 이용한 설치
conda create -n qc_env
conda activate qc_env
conda install -c bioconda fastqc

# Ubuntu/Debian
sudo apt-get install fastqc
```

### 주요 기능
1. 기본 통계
2. 시퀀스 품질 점수
3. GC 함량 분포
4. 시퀀스 길이 분포
5. 중복 시퀀스 분석
6. 오버리프리젠티드 시퀀스
7. K-mer 분석

### 실행 방법
```bash
# 단일 파일 분석
fastqc sample.fastq

# 여러 파일 동시 분석
fastqc sample1.fastq sample2.fastq

# 결과 저장 디렉토리 지정
fastqc -o output_dir sample.fastq

# 압축 파일 직접 분석
fastqc sample.fastq.gz
```

## Cutadapt
Cutadapt는 어댑터 제거 및 품질 필터링을 수행하는 도구입니다.

### 설치 방법
```bash
# pip를 이용한 설치
pip install cutadapt

# Conda를 이용한 설치
conda install -c bioconda cutadapt
```

### 주요 기능
1. 어댑터 시퀀스 제거
2. 품질 기반 트리밍
3. 길이 기반 필터링
4. Paired-end 리드 처리

### 기본 사용법
```bash
# 3' 어댑터 제거
cutadapt -a AGATCGGAAGAG -o output.fastq input.fastq

# 품질 점수 기반 트리밍
cutadapt -q 20 -o output.fastq input.fastq

# 최소 길이 필터링
cutadapt -m 50 -o output.fastq input.fastq

# Paired-end 데이터 처리
cutadapt -a ADAPT1 -A ADAPT2 -o out1.fastq -p out2.fastq in1.fastq in2.fastq
```

## 분석 파이프라인
일반적인 전처리 파이프라인은 다음 단계로 구성됩니다:

1. 원시 데이터 품질 확인 (FastQC)
2. 어댑터 제거 및 품질 필터링 (Cutadapt)
3. 처리된 데이터 품질 재확인 (FastQC)

## 결과 해석

### FastQC 결과 해석
- Per base sequence quality: 대부분의 베이스가 Q30 이상
- Per sequence GC content: 정규 분포 형태
- Sequence duplication levels: 낮은 중복율
- Overrepresented sequences: 어댑터 오염 확인

### Cutadapt 결과 해석
- 처리된 리드 수
- 제거된 어댑터 수
- 필터링된 리드 수
- 최종 데이터의 통계