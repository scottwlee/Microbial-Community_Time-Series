# Microbial-Community_Time-Series (KR)
미생물 군집 시계열 데이터(Microbial Community Time-Series Data)의 전처리, 환경 변수 분석 및 시계열 모델링을 위한 R 스크립트 저장소입니다.

## 📂 저장소 구조 및 파이프라인 (Workflow & File Structure)

### 1. 데이터 선별 및 전처리 (Data Selection & Preprocessing)
- `01_ENT&AA-Selection.R`: 타겟 타겟군(ENT, AA 등) 및 데이터 선별 스크립트
- `02_Temperature-Trend.R`: 수온 및 기온 등 주요 환경 변수의 트렌드 분석
- `04_Data Integration&Processing.R`: 이질적 데이터셋 통합 및 표준화/전처리 파이프라인

### 2. 시계열 분석 모듈 (Time-Series Analysis)
- `05_Time-Series_3-1.R` ~ `05_Time-Series_3-14.R`: 미생물 군집 동태 분석, 패턴 모델링 및 시각화를 단계별로 수행하는 모듈 스크립트 (3-1부터 3-14까지 순차적 구성)

### 3. 헬퍼 함수 및 공통 모듈 (Helper Functions)
- `F2_HelperFunction_MDR_Block.R`: 시계열 블록 연산 및 주요 데이터 처리에 필요한 사용자 정의 함수 모음

================================================================================

# Microbial-Community_Time-Series (US)
A collection of R scripts designed for preprocessing, environmental trend analysis, and time-series modeling of microbial community datasets.

## 📂 Workflow & File Structure

### 1. Data Selection & Preprocessing
- `01_ENT&AA-Selection.R`: Data filtering and selection script for target groups (e.g., ENT, AA).
- `02_Temperature-Trend.R`: Analysis of environmental temperature trends.
- `04_Data Integration&Processing.R`: Data integration and preprocessing pipeline across datasets.

### 2. Time-Series Analysis Modules
- `05_Time-Series_3-1.R` ~ `05_Time-Series_3-14.R`: Step-by-step modular scripts (3-1 to 3-14) performing time-series analysis, dynamic modeling, and visualization.

### 3. Helper Functions
- `F2_HelperFunction_MDR_Block.R`: Custom utility functions for MDR block operations and statistical calculations.
