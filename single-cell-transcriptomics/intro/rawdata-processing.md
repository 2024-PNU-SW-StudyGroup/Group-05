문서가 매우 길어서 중요한 섹션별로 나누어 번역하도록 하겠습니다.



# Raw Data Processing (원시 데이터 처리)

## 개요 및 소개

이 섹션에서는 단일세포 및 단일핵 RNA 시퀀싱(sc/snRNA-seq) 데이터의 "전처리(preprocessing)"라고 일반적으로 불리는 기본적인 문제들을 다룹니다. 이는 일반적인 용어이지만, 데이터를 처리하고 표현하는 방법에 대한 중요한 결정을 포함하는 여러 단계로 구성되어 있어 약간의 오해의 소지가 있을 수 있습니다.

여기서는 이 단계를 "raw data processing(원시 데이터 처리)"라고 부르며, 다음과 같은 단계들을 포함합니다:

1. Lane-demultiplexed FASTQ 파일에서 시작
2. 각 세포 내의 각 유전자로부터 발생하는 고유 분자의 추정 수를 나타내는 count matrix 생성
3. 각 분자의 스플라이싱 상태에 따른 분류 

## Raw Data Quality Control (원시 데이터 품질 관리)

FASTQ 파일을 얻은 후에는 [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)와 같은 QC 도구를 사용하여 리드 품질을 빠르게 진단할 수 있습니다. 주요 평가 항목:

- Basic Statistics (기본 통계)
- Per Base Sequence Quality (염기별 시퀀스 품질)
- Per Tile Sequence Quality (타일별 시퀀스 품질) 
- Per Sequence Quality Scores (시퀀스별 품질 점수)
- Per Base Sequence Content (염기별 시퀀스 구성)
- Per Sequence GC Content (시퀀스별 GC 함량)
- Per Base N Content (염기별 N 함량)
- Sequence Length Distribution (시퀀스 길이 분포)
- Sequence Duplication Levels (시퀀스 중복 수준)
- Overrepresented Sequences (과다 표현된 시퀀스)
- Adapter Content (어댑터 함량)

## Alignment and Mapping (정렬 및 매핑)

매핑 또는 정렬은 단일세포 원시 데이터 처리의 기본 단계입니다. 시퀀싱된 각 단편(예: 유전체 또는 전사체 loci와 유사한 리드)의 잠재적 위치를 결정하는 과정을 의미합니다.

시퀀싱 프로토콜에 따라 결과 시퀀스 파일은 다음을 포함합니다:
- Cell barcodes (CB, 세포 바코드)
- Unique molecule identifier (UMI, 고유 분자 식별자)
- Raw cDNA sequence (원시 cDNA 시퀀스)

네, 계속해서 번역하겠습니다.

## 매핑 도구 및 접근법

현재 가장 널리 사용되는 단일세포 RNA-seq 데이터 처리 전용 도구들:
- `Cell Ranger` (10x Genomics의 상용 소프트웨어)
- `zUMIs` 
- `alevin`
- `RainDrop`
- `kallisto|bustools`
- `STARsolo`
- `alevin-fry`

이러한 도구들은 크게 두 가지 기준으로 분류할 수 있습니다:
1. 수행하는 매핑의 유형
2. 리드를 매핑하는 참조 시퀀스의 유형

### 매핑의 유형

세 가지 주요 매핑 알고리즘 유형:

1. "Spliced alignment(스플라이스된 정렬)"
   - 리드가 여러 스플라이스 접합부를 걸쳐 정렬될 수 있음
   - 게놈에 매핑할 때 주로 사용

2. "Contiguous alignment(연속 정렬)"
   - 참조의 연속된 부분 문자열과의 정렬 탐색
   - 작은 삽입과 결실은 허용되나 큰 간격은 허용되지 않음

3. "Lightweight mapping(경량 매핑)"
   - k-mer나 다른 유형의 정확한 매치만을 기반으로 함
   - 속도가 빠르지만 정렬 품질 평가가 어려움

### 참조 시퀀스 유형

매핑에 사용되는 세 가지 주요 참조 시퀀스 유형:

1. **전체 게놈 (Full genome)**
- 장점:
  - 모든 게놈 위치의 리드 확인 가능
  - 인트론과 주석되지 않은 영역의 리드도 포함
- 단점:
  - 큰 메모리 요구사항 (인간 게놈의 경우 32GB 이상)
  - 계산 비용이 높음

2. **주석된 전사체 (Annotated transcriptome)**
- 장점:
  - 더 작은 참조 크기로 계산 자원 절약
  - 스플라이스된 정렬이 필요 없음
- 단점:
  - 스플라이스된 전사체 외부의 리드 캡처 불가
  - Single-nucleus 데이터에는 부적합

3. **증강된 전사체 (Augmented transcriptome)**
- 장점:
  - 게놈보다 작은 인덱스 크기
  - 연속 정렬만 필요
  - 스플라이스된 전사체 외부의 리드도 포함
- 단점:
  - 전체 게놈 매핑보다는 덜 포괄적

## Cell Barcode Correction (세포 바코드 보정)

Droplet 기반 단일세포 분리 시스템의 주요 오류 원인:
1. Doublet/Multiplet: 하나의 바코드가 두 개 이상의 세포와 연관
2. Empty Droplet: 세포가 없는 방울에 ambient RNA가 태그됨
3. Sequence error: PCR 증폭이나 시퀀싱 과정의 오류

바코드 보정을 위한 일반적인 전략:
1. 알려진 바코드 목록과 대조 보정
2. Knee/elbow 기반 방법
3. 사용자가 제공한 예상 세포 수 기반 필터링
4. 강제된 유효 세포 수 기반 필터링
