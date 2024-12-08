(introduction:analysis-frameworks)=

# 분석 프레임워크와 도구

## 단일세포 분석 프레임워크와 컨소시엄

앞서 설명한 대로 카운트 매트릭스를 얻은 후에는 탐색적 데이터 분석 단계가 시작됩니다. 데이터의 크기와 복잡성으로 인해 특수한 도구가 필요합니다. 초기 단일세포 분석에서는 사용자 정의 스크립트로 데이터를 분석했지만, 이제는 이를 위한 전문 프레임워크가 존재합니다. 가장 인기 있는 3가지 옵션은 R 기반의 Bioconductor와 Seurat 생태계, 그리고 Python 기반의 scverse 생태계입니다. 이들은 사용하는 프로그래밍 언어뿐만 아니라 기본 데이터 구조와 사용 가능한 특수 분석 도구에서도 차이가 있습니다.

Bioconductor는 단일세포를 포함한 다양한 생물학적 분석을 위한 무료 오픈소스 소프트웨어를 개발, 지원 및 공유하는 프로젝트입니다. 일관된 개발자/사용자 경험과 사용자 친화적인 매뉴얼이 포함된 광범위한 문서화가 Bioconductor의 가장 큰 장점입니다. Seurat은 단일세포 데이터 분석을 위해 특별히 설계된 잘 알려진 R 패키지입니다. 다중모달 및 공간 데이터를 포함한 모든 분석 단계에 대한 도구를 제공합니다. 잘 작성된 매뉴얼과 큰 사용자 기반이 Seurat의 특징입니다.

하지만 두 R 옵션 모두 매우 큰 데이터셋(50만 개 이상의 세포)에서는 확장성 문제가 발생할 수 있어, Python 기반 커뮤니티가 scverse 생태계를 개발하게 되었습니다. scverse는 생명과학 분야의 기초 도구 개발에 중점을 둔 조직이자 생태계입니다. 확장성, 확장 가능성, 기존 Python 데이터 및 머신러닝 도구와의 강력한 상호운용성이 scverse 생태계의 장점입니다.

세 생태계 모두 프레임워크 간 상호운용성을 위한 많은 노력을 기울이고 있습니다. 이는 "상호운용성" 장에서 논의될 것입니다. 이 책은 항상 해당 질문에 가장 적합한 도구에 초점을 맞추므로 특정 문제에 대해 위에서 언급한 생태계를 혼합하여 사용할 것입니다. 하지만 모든 분석의 기초는 다음 두 가지 이유로 scverse 생태계가 될 것입니다:

1. 이 책에서는 생태계와 프로그래밍 언어를 정기적으로 전환하지만, 일관된 데이터 구조와 도구 사용은 독자가 구현 세부사항이 아닌 개념에 집중하는 데 도움이 됩니다.
2. Bioconductor 생태계만을 다루는 훌륭한 책이 이미 존재합니다. Bioconductor로 단일세포 분석만 배우고 싶은 사용자는 해당 책을 읽어보시기 바랍니다.

다음 섹션에서는 scverse 생태계를 자세히 소개하고 가장 중요한 데이터 구조를 중심으로 주요 개념을 설명합니다. 이 소개는 데이터 구조와 프레임워크의 모든 측면을 다룰 수는 없습니다. 사용 가능한 모든 분석 함수를 소개하는 것은 범위를 벗어납니다. 따라서 필요한 경우 각 프레임워크의 튜토리얼과 문서를 참조하시기 바랍니다.

## AnnData로 단일모달 데이터 저장하기

앞서 논의한 대로, 정렬과 유전자 주석 후에 유전체 데이터는 일반적으로 특성 행렬로 요약됩니다. 이 행렬은 `관측치 수 x 변수 수` 형태를 가지며, scRNA-seq의 경우 관측치는 세포 바코드이고 변수는 주석이 달린 유전자입니다. 분석 과정에서 이 행렬의 관측치와 변수에는 계산된 측정값(예: 품질 관리 메트릭스 또는 잠재 공간 임베딩)과 사전 지식(예: 공여자 출처 또는 대체 유전자 식별자)이 주석으로 달립니다. scverse 생태계에서는 AnnData를 사용하여 데이터 행렬을 이러한 주석과 연결합니다. 빠르고 메모리 효율적인 변환을 위해 AnnData는 희소 행렬과 부분 읽기도 지원합니다.

AnnData는 R 생태계의 데이터 구조(예: Bioconductor의 SummarizedExperiment 또는 Seurat의 객체)와 광범위하게 유사하지만, R 패키지는 전치된 특성 행렬을 사용합니다.

핵심적으로 AnnData 객체는 `X`에 희소 또는 밀집 행렬(scRNA-Seq의 경우 카운트 행렬)을 저장합니다. 이 행렬은 `obs_names x var_names` 차원을 가지며, obs(=관측치)는 세포 바코드에 해당하고 var(=변수)는 유전자 식별자에 해당합니다. 이 행렬 `X`는 각각 세포와 유전자의 주석을 저장하는 Pandas DataFrame인 `obs`와 `var`로 둘러싸여 있습니다. 또한 AnnData는 관측치(`obsm`) 또는 변수(`varm`)에 대한 전체 계산 행렬을 해당 차원과 함께 저장합니다. 세포와 세포 또는 유전자와 유전자를 연결하는 그래프와 같은 구조는 일반적으로 `obsp`와 `varp`에 저장됩니다. 다른 슬롯에 맞지 않는 비구조화된 데이터는 `uns`에 비구조화된 데이터로 저장됩니다. `layers`에 `X`의 더 많은 값을 저장하는 것도 가능합니다. 예를 들어 원시 정규화되지 않은 카운트 데이터를 `counts` 레이어에 저장하고 정규화된 데이터를 이름 없는 기 레이어에 저장하는 것이 사용 사례가 될 수 있습니다.
AnnData는 주로 단일모달(예: scRNA-Seq만) 데이터용으로 설계되었습니다. 하지만 이 장에서 나중에 다룰 MuData와 같은 AnnData의 확장은 다중모달 데이터의 효율적인 저장과 접근을 가능하게 합니다.

## AnnData 객체 초기화하기
이 섹션은 AnnData의 "시작하기" 튜토리얼에서 영감을 받았습니다: https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html
희소 카운트 정보가 포함된 간단한 AnnData 객체를 만들어보겠습니다. 이는 예를 들어 유전자 발현 카운트를 나타낼 수 있습니다. 먼저 필요한 패키지들을 임포트합니다.

```python
import numpy as np
import pandas as pd
import anndata as ad
```

이제 무작위 데이터로 AnnData 객체를 만들어보겠습니다:
```python
# 관측치(세포)와 변수(유전자)의 수를 정의합니다
n_obs = 1000
n_vars = 100

# 무작위 데이터 생성
data = np.random.binomial(100, 0.005, (n_obs, n_vars))
obs_names = [f'cell_{i:d}' for i in range(n_obs)]
var_names = [f'gene_{i:d}' for i in range(n_vars)]

# AnnData 객체 생성
adata = ad.AnnData(
    X=data,
    obs=dict(obs_names=obs_names),
    var=dict(var_names=var_names),
    dtype='int32',
)
```

## 주석 추가하기
AnnData 객체에 주석을 추가하는 것은 매우 간단합니다. 관측치(세포)나 변수(유전자)에 대한 주석을 추가할 수 있습니다:

```python
# 관측치에 대한 그룹 정보 추가
adata.obs['group'] = np.random.choice(['A', 'B', 'C'], size=n_obs)

# 변수에 대한 추가 정보
adata.var['highly_variable'] = np.random.choice([True, False], size=n_vars)
```

## 데이터 저장 및 로드
AnnData 객체는 h5ad 형식으로 저장할 수 있습니다.

```python
adata.write('my_data.h5ad')

# 데이터 다시 로드하기
adata = ad.read_h5ad('my_data.h5ad')
```

# MuData로 다중모달 데이터 저장하기
최근 단일세포 기술의 발전으로 동일한 세포에서 여러 가지 특성을 동시에 측정할 수 있게 되었습니다. 예를 들어, CITE-seq은 전사체와 단백질 수준을 동시에 측정할 수 있습니다. 이러한 다중모달 데이터를 효율적으로 저장하고 분석하기 위해 MuData가 개발되었습니다.
MuData는 여러 개의 AnnData 객체를 하나의 컨테이너에 저장합니다. 각 AnnData 객체는 서로 다른 모달리티(예: RNA, ATAC, 단백질 등)를 나타냅니다. 이를 통해 각 모달리티의 데이터를 독립적으로 처리하면서도 모달리티 간의 관계를 분석할 수 있습니다.
## 설치
MuData는 다음과 같이 설치할 수 있습니다:
```bash
pip install mudata
```
## MuData 객체 생성하기
RNA-seq과 ATAC-seq 데이터를 포함하는 MuData 객체를 만드는 예시를 살펴보겠습니다:

```python
import mudata as md

# RNA와 ATAC 데이터에 대한 별도의 AnnData 객체 생성
rna = ad.AnnData(np.random.binomial(100, 0.005, (n_obs, n_vars)))
atac = ad.AnnData(np.random.binomial(100, 0.005, (n_obs, n_vars)))

# MuData 객체 생성
mdata = md.MuData({'rna': rna, 'atac': atac})
```

더 자세한 내용은 muon API 참조 문서(https://muon.readthedocs.io/en/latest/api/index.html)와 muon 튜토리얼(https://muon-tutorials.readthedocs.io/en/latest/)을 참조하시기 바랍니다.