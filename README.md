# OralProbPhage

构建"口腔细菌组–病毒组–宿主"知识层，为口腔益生菌噬菌体联合疗法提供数据基础。

## 研究背景

本项目针对口腔疾病（牙周病、龋病、口臭）的益生菌疗法瓶颈：
- **益生菌效果不佳**：S. salivarius K12/M18 在部分人群中无效
- **核心假说**：口腔内源噬菌体（尤其是针对益生菌的噬菌体）可能干扰益生菌定植
- **联合策略**：筛选"益生菌 + 其噬菌体"组合，实现协同定植与病原抑制

## Quick start

```r
# 设置工作流目录
qproj::proj_use_workflow("analyses")

# 创建分析模块（已创建）
qproj::use_qmd("010-schema-design", "analyses")
qproj::use_qmd("020-data-source-registry", "analyses")
qproj::use_qmd("030-preprocessing-pipeline", "analyses")
qproj::use_qmd("040-m1a-data-processing", "analyses")
qproj::use_qmd("050-m1a-differential-analysis", "analyses")
qproj::use_qmd("060-m1a-vsc-association", "analyses")
qproj::use_qmd("070-m2a-virome-processing", "analyses")
qproj::use_qmd("080-m2a-host-prediction", "analyses")
qproj::use_qmd("090-m2a-mapping", "analyses")
qproj::use_qmd("100-knowledge-integration", "analyses")
qproj::use_qmd("110-evaluation", "analyses")
qproj::use_qmd("990-manuscript", "analyses")
```

## Project structure

- `analyses/`: qproj 分析模块
  - `010-schema-design.qmd`: 统一对象层（M0a）设计，定义 5 个核心对象
  - `020-data-source-registry.qmd`: 公共队列与数据源登记
  - `030-preprocessing-pipeline.qmd`: 统一预处理流程生成
  - `040-m1a-data-processing.qmd`: M1a 细菌组信号层数据抓取与重处理
  - `050-m1a-differential-analysis.qmd`: M1a 跨队列差异分析与 CandidateMicrobe 生成
  - `060-m1a-vsc-association.qmd`: M1a 与 VSC/临床指标关联分析
  - `070-m2a-virome-processing.qmd`: M2a 口腔 virome 重处理与 vOTU catalog
  - `080-m2a-host-prediction.qmd`: M2a 病毒–宿主预测与 PhageHostLink 构建
  - `090-m2a-mapping.qmd`: M2a 与 M1a 结果映射
  - `100-knowledge-integration.qmd`: 知识库整合
  - `110-evaluation.qmd`: 内部评估与文档化
  - `990-manuscript.qmd`: 手稿（数据论文）
- `analyses/data/`: 数据目录
  - `00-raw/`: 原始数据（只读）
  - 其他目录为各模块产出
- `R/utils.R`: 辅助函数（load_db/close_db/write_table_db）
- `analyses/pr0003.sqlite`: 项目 SQLite 数据库（git 忽略）
- `analyses/background.md`: Notion 研究方案背景
- `IMPLEMENTATION_PLAN.md`: 详细实施方案

## Current status

**Phase**: 初始化完成，准备开始阶段 1（基础设施与 Schema 设计）

**已完成**:
- ✅ qproj 框架初始化
- ✅ 12 个分析模块已创建
- ✅ SQLite 数据库准备（pr0003.sqlite）
- ✅ 背景文档已保存（analyses/background.md）
- ✅ 实施方案已制定（IMPLEMENTATION_PLAN.md）

**下一步**:
1. 执行 `010-schema-design.qmd`：定义 5 个核心对象
2. 执行 `020-data-source-registry.qmd`：登记公共队列
3. 开始 M1a 细菌组信号分析

## Core Objectives

### P0: 总体设计
构建一个可迭代更新的"口腔细菌组–病毒组–宿主"知识层，输出标准化的候选菌株和噬菌体关系信息。

### M0a: 统一对象层
定义 5 个核心对象：Sample、TaxonProfile、CandidateMicrobe、PhageHostLink、EvidencePacket

### M1a: 细菌组信号层
在多个队列中识别口腔疾病（牙周炎、牙龈炎、龋病、口臭）的健康/疾病富集菌

### M2a: 病毒–宿主层
构建口腔 virome catalog，预测噬菌体–宿主关系，与 M1a 结果映射

## Key Datasets

### 细菌组
- Periodontitis 通用标志物整合研究
- Gingivitis vs 健康 shotgun 队列
- 现代与历史口腔微生物组大合集
- 口臭–牙周病–VSC 综述

### 病毒组
- RA 研究中的口腔–肠道联合 virome
- 健康人群口腔 DNA virome
- Oral virome atlas 计划

## Tools & Resources

- **数据库**: eHOMD（口腔细菌）、IMG/VR（病毒）
- **预处理**: vpipe（病毒组）、metaSPAdes（组装）
- **分析**: ViWrap、ViOTUcluster（病毒组）、DESeq2（差异分析）
- **宿主预测**: 多工具集成（CRISPR、相似性、代谢互补）

## Expected Outputs

1. **数据论文**: 口腔细菌–病毒公共资源整合（目标：Scientific Data）
2. **知识库原型**: 口腔微生物–病毒知识库（可为专利打基础）
3. **方法学论文**: 统一预处理 pipeline（目标：Microbiome）

## Citations

研究方案来源：[Notion 页面](https://www.notion.so/p0-P0-M0-M2--35725833f8a781b799d2ee2bf8413161)

核心参考文献：
- eHOMD: Expanded Human Oral Microbiome Database
- Universal oral microbiome signature for periodontitis (PMC11316925)
- Dysbiotic oral and gut viromes in RA (biorxiv 2021.03.05.434018)
- ViWrap: Modular pipeline for virome analysis (iMeta 2023)
- ViOTUcluster: High-throughput viromics pipeline (iMeta 2024)

---

**Last updated**: 2026-05-05
