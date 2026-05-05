# pr0003-OralProbPhage 实施方案

## 项目概述
**目标**：构建"口腔细菌组–病毒组–宿主"知识层，输出标准化候选菌株和噬菌体关系信息。

**应用方向**：口腔益生菌（S. salivarius K12/M18）与噬菌体联合疗法，重点针对牙周病、龋病、口臭。

---

## 实施阶段（基于 Notion 方案）

### 阶段 1：基础设施与 Schema 设计（Week 1-2）

#### 任务清单
1. **创建 qproj 分析模块**
   - `010-schema-design`：定义 5 个核心对象（Sample, TaxonProfile, CandidateMicrobe, PhageHostLink, EvidencePacket）
   - `020-data-source-registry`：登记公共队列与数据源
   - `030-preprocessing-pipeline`：统一预处理流程

2. **建立数据模型**
   - 参考 eHOMD/HOMD 设计 stable ID 方案
   - 为 vOTU 设计稳定 ID 系统（借鉴 HMT 思路）
   - 输出 YAML/JSON schema 原型

3. **初始化 SQLite 数据库**
   - 创建 `pr0003.sqlite`
   - 建表：`step01_sources`, `step02_samples`, `step03_taxon_profiles` 等

#### 执行命令
```r
# 在 R 中执行
qproj::proj_use_workflow("analyses")
qproj::use_qmd("010-schema-design", "analyses")
qproj::use_qmd("020-data-source-registry", "analyses")
qproj::use_qmd("030-preprocessing-pipeline", "analyses")
```

---

### 阶段 2：细菌组信号分析（M1a）（Week 3-4）

#### 任务清单
1. **数据抓取与重处理（Task M1a-1）**
   - 输入：论文补充数据中的 SRA 号
   - 数据源：
     - Periodontitis 通用标志物研究
     - Gingivitis shotgun 队列
     - 现代与历史口腔微生物组大合集
     - 口臭–牙周病–VSC 综述
   - 动作：下载 fastq、质控、分类学分析（使用 vpipe 或类似工具）

2. **跨队列差异分析（Task M1a-2）**
   - 对每个疾病表型分别做：
     - periodontitis（牙周炎）
     - gingivitis（牙龈炎）
     - caries（龋病）
     - halitosis（口臭）
   - 识别健康富集菌（health_enriched）和疾病富集菌（disease_enriched）
   - 输出：CandidateMicrobe 对象

3. **VSC/临床指标关联（Task M1a-3）**
   - 对带 VSC（挥发性硫化合物）或炎症指标的队列做关联分析
   - 重点关注：
     - S. salivarius K12/M18 群
     - 口臭相关病原（P. gingivalis 等）

#### 创建 qmd 模块
```r
qproj::use_qmd("040-m1a-data-processing", "analyses")
qproj::use_qmd("050-m1a-differential-analysis", "analyses")
qproj::use_qmd("060-m1a-vsc-association", "analyses")
```

---

### 阶段 3：病毒–宿主层分析（M2a）（Week 5-6）

#### 任务清单
1. **口腔 virome 重处理（Task M2a-1）**
   - 输入：RA 研究和健康队列的牙菌斑/唾液 metagenome
   - 工具：ViWrap 或 ViOTUcluster
   - 动作：病毒识别、vOTU clustering、质量过滤
   - 输出：vOTU catalog

2. **病毒–宿主预测（Task M2a-2）**
   - 对每个 vOTU 综合多种证据：
     - 序列相似性（BLAST/DIAMOND）
     - CRISPR spacer 匹配
     - 代谢互补分析
   - 输出：PhageHostLink 对象

3. **与 M1a 结果映射（Task M2a-3）**
   - 将 CandidateMicrobe（health/disease enriched）映射到 PhageHostLink
   - 分析噬菌体-宿主兼容性
   - 输出：兼容性矩阵（供后续 M6+ 使用）

#### 创建 qmd 模块
```r
qproj::use_qmd("070-m2a-virome-processing", "analyses")
qproj::use_qmd("080-m2a-host-prediction", "analyses")
qproj::use_qmd("090-m2a-mapping", "analyses")
```

---

### 阶段 4：整合与迭代（Week 7-8）

#### 任务清单
1. **构建统一知识库**
   - 整合 M1a 和 M2a 结果
   - 生成 EvidencePacket 对象
   - 更新 SQLite 数据库

2. **第一次迭代评估**
   - 检查数据覆盖度
   - 识别缺失数据源
   - 优化预处理流程

3. **准备数据论文**
   - 撰写方法部分
   - 生成图表
   - 开始撰写 manuscript（创建 `100-manuscript.qmd`）

#### 创建 qmd 模块
```r
qproj::use_qmd("100-knowledge-integration", "analyses")
qproj::use_qmd("110-evaluation", "analyses")
qproj::use_qmd("990-manuscript", "analyses")
```

---

## 工具与技术栈

### 核心工具
- **数据预处理**：vpipe（病毒组）、metaSPAdes/MEGAHIT（组装）
- **分类学分析**：BLAST、DIAMOND、Kraken2
- **病毒识别**：ViWrap、ViOTUcluster、CheckV
- **宿主预测**：多种工具集成（CRISPR、相似性、代谢互补）
- **统计分析**：R（tidyverse、DESeq2、phyloseq）

### 数据库
- **eHOMD**：口腔细菌参考数据库
- **RefSeq/GenBank**：噬菌体基因组
- **IMG/VR**：病毒参考数据库

---

## 输出产物

### 数据产物
1. `pr0003.sqlite`：包含所有核心对象的关系型数据库
2. `data/010-schema-design/`：Schema 定义文件（YAML/JSON）
3. `data/020-data-source-registry/`：数据源登记表
4. `data/040-m1a-data-processing/`：处理后的细菌组 profiles
5. `data/070-m2a-virome-processing/`：vOTU catalog
6. `data/080-m2a-host-prediction/`：PhageHostLink 对象

### 可发表产出
1. **数据论文**：口腔细菌–病毒公共资源整合（目标期刊：Scientific Data、GigaScience）
2. **知识库原型**：公开的口腔微生物–病毒知识库（可为后续专利打基础）
3. **方法学论文**：统一预处理 pipeline（目标期刊：Microbiome、Cell Reports Methods）

---

## 时间线总览

| 周次 | 阶段 | 主要产出 |
|------|------|---------|
| 1-2 | 基础设施 | Schema 设计、数据源登记、预处理流程 |
| 3-4 | M1a 细菌组 | CandidateMicrobe 对象、差异分析结果 |
| 5-6 | M2a 病毒组 | vOTU catalog、PhageHostLink 对象 |
| 7-8 | 整合迭代 | 知识库原型、第一次评估、manuscript 初稿 |
| 9-12 | 扩展迭代 | 整合更多数据源、完善分析、准备投稿 |

---

## 关键成功指标（KPI）

1. **数据覆盖度**：整合 ≥5 个公共队列，覆盖 ≥1000 个样本
2. **稳健信号**：在 ≥3 个独立队列中重复的健康/疾病富集菌
3. **宿主预测质量**：≥70% 的 vOTU 有 ≥2 种证据支持
4. **知识库完整性**：5 类核心对象均有数据，且关系完整
5. **可重现性**：从原始数据到最终结果的完整 pipeline

---

**更新时间**：2026-05-05
