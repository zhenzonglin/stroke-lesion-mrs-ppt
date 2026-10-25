# Stroke lesion ROI and longitudinal mRS summaries

只读提取SAS临床队列的长期mRS分布，并提供独立的急性梗死病灶ROI统计入口。输出用于后续交互式PowerPoint制作；本仓库目前不包含PPT生成器，也不包含患者数据或真实统计结果。

## 当前轨迹图：直接按SAS构建两组缺血性卒中人群

轨迹图的两个人群均直接来自SAS，不检查病灶文件是否存在、不要求影像匹配、不根据空掩膜排除患者：

- **全体缺血性卒中**：`D_DIAG=1`。出院mRS缺失或0–1分仍保留在该总体中。
- **出院时致残性缺血性卒中**：上述人群中，`D_MRS`为有效整数2–6者。

```bash
git pull --ff-only
conda activate wmh-hcy
python run_clinical_trajectories.py --config config.v2.local.json
```

此入口仅使用配置中的两份SAS路径、ID宽度和SAS编码设置，不使用`lesion_root`、图谱或既往影像运行目录。早期与长期记录按`code_n`连接并保留ID并集，再按`D_DIAG=1`筛选；没有诊断依据的记录不推定为缺血性卒中。各患者不需六次随访完整，缺失记录不插补。

输出新建`outputs/clinical_trajectories_<UTC>/share.zip`。其中`clinical_trajectories.json`提供筛选流程、出院mRS频数、两个人群的六次精确分组计数、比例和有效分母。终端PAGE及`share/screenshots.txt`可供截图回传。患者ID和逐人状态只留在工作站的`private/`，分享包不包含它们。

分组仍为mRS 0–1恢复良好、2–3功能残疾、4–6功能丧失；死亡仍仅按mRS 6定义。每次曲线的分母为该人群当次有效mRS 0–6记录，三组未舍入比例在有效分母大于0时合计100%。出院6分后的0–5分访视标为冲突，不进入该次有效分母；缺失不前填。

后续两页PPT分别展示上述两组轨迹，每页配同一张全影像队列的脑区累及图。脑图沿用全影像队列的分母，轨迹采用本次SAS队列的分母，二者分别标注。基底节与内囊合并为一个脑区时，累及人数仍须使用患者去重并集，不能相加或用四区联合人数替代。以下影像统计入口保留，可用于脑区图；无需先运行它们才能生成临床轨迹。

## 影像队列亚组：出院时致残性缺血性卒中

在原v2运行的有效、非空且匹配临床的影像队列中，选择 **D_DIAG=1 且 D_MRS为有效整数2–6** 的患者。D_MRS是出院mRS，不用3个月mRS代替；TIA、出院mRS 0–1、出院mRS缺失或非法值不纳入。按`D_MRS≥2`的原条件保留6分，并单独报告出院6分人数；不额外增加出院存活条件。

```bash
git pull --ff-only
conda activate wmh-hcy
python supplement_discharge_cohort.py --run-dir outputs/<原v2运行时间戳>
```

这里的来源目录应含原始`share/ppt_data.json`和`private/`，不是`deep_union_...`补统计目录。命令从原配置读取早期SAS的`code_n、D_DIAG、D_MRS`；读取前后均核对SAS与原运行的SHA256。缺少指定字段或ID不合格时停止，不猜测替代字段。

原运行未缓存逐人ROI归属，因此需重新读取原匹配掩膜一次。程序同时计算13个独立ROI与四区体素并集，复核原有效患者集合、排除原因及总体/13区全部六次mRS频数，随后对入选患者重新汇总。其他13区保留原显示编号和顺序；这些编号代表原队列排序，不宣称是新亚组的累及率排序。放射冠、基底节、内囊、半卵圆中心联合人群仍单独统计，每名患者只计一次。

新输出位于`outputs/discharge_<UTC>/`，不会改写原运行。分享包包含筛选流程、出院mRS频数、总体和联合人群曲线、13区累及比例及六次残疾/死亡率。`private/`保留逐人筛选决定和派生状态，不进入分享包。优先回传`share.zip`；也可发送全部PAGE截图。

新亚组的累及比例分母为入选患者总数；每次功能分组及残疾/死亡率的分母为该人群当次有效mRS 0–6人数。缺失不作为曲线类别；缺失者仍在亚组总数中，不进行完整病例筛选。三组仍为0–1、2–3、4–6，有效分母大于0时未舍入比例合计100%。出院mRS 6后出现随访0–5的记录另标冲突，不解释为恢复；缺失不前填为死亡。

此步骤只准备新队列数据。后续PPT沿用已确认的版式与逐时间点三组同步延伸动画，替换队列标题、人数和比例。卒中前及早期轨迹仍是无实测数值的示意，3个月起曲线高度对应实际比例。

## 全体队列：保留13区，另补算四区联合人群

保留全部13个独立ROI，按累及人数降序编号。另取**放射冠＋基底节＋内囊＋半卵圆中心**的患者并集，单独展示其功能分组随访曲线；该联合人群不替换原ROI。

```bash
git pull --ff-only
conda activate wmh-hcy
python supplement_deep_union.py --run-dir outputs/<原v2运行时间戳>
```

此命令不重读SAS，复用原运行保存的临床状态，并重新核查四个ROI及其体素并集。只有有效患者集合、排除原因、四个组成区和总体的六次mRS频数全部复现时才输出完成结果。每个患者在联合人群中仅计一次，检查包括重叠分布及四集合包含排除关系。

输出新建`outputs/deep_union_<UTC>/share.zip`。其中`deep_union_supplement.json`的`union_curve`给出六次0–1、2–3、4–6三组计数、比例及有效分母，`ranked_rois_13`和`ranked_outcome_rates_13`完整保留原13区。分享包只含群体汇总。`private/`需留在工作站，不能仅凭各ROI边际人数或截图推算联合结局。

曲线纵轴若表示百分比，应按实际占比升降，三类合计100%；卒中前与急性期没有实测mRS时，只能另标无数值示意。群体比例连线不代表个体恢复或转移。

## 可选补统计：基底节-内囊并集与12区排序

基底节和内囊可涉及同一患者，不能相加人数或平均两区结局比例。补算命令从原运行的`private/derived_states.csv`复用已经核对的六次mRS，不重读SAS；从`private/mask_inventory.csv`定位原掩膜，只核查基底节、内囊及其体素并集。

```bash
git pull --ff-only
conda activate wmh-hcy
python supplement_bgic.py --run-dir outputs/<原v2运行时间戳>
```

仅当原有效患者集合、排除原因、两组成区及总体的全部六次mRS分布均复现时，输出才标为完成。逐个患者在并集中仅计一次，0.1 mL敏感性阈值也作用于体素并集，不能合并两区各自达到阈值的人群。

补算结果新建`outputs/bgic_<UTC>/`，不覆盖旧运行。`share.zip`含`bgic_supplement.json`、合并区六次结局、12区累及排序及对应结局CSV。两页核对分别显示并集/交集人数和六次结局、全部12区排序。分享包不含患者ID。其他11区沿用原运行汇总。

原运行没有记录每张影像的SHA256，因此补算核对文件大小、有效患者集合和组成区频数，不宣称已证明所有影像逐字节未变。`private/`必须仍保留在工作站；仅有截图或分享包无法恢复并集患者名单。

## 已有工作站：更新13区与六次残疾/死亡率

v2将尾状核、壳核和苍白球合并为基底节，丘脑单列，增加内囊、放射冠、半卵圆中心和胼胝体。六次随访为3、12、24、36、48、60个月。脑图旁的**残疾/死亡率 = 该ROI人群该次mRS 2–6人数 / 该次有效mRS 0–6人数**，各次有效人数可以不同。mRS 6仍为死亡。

```bash
cd /data/usersdir/linzhenzong/stroke-lesion-mrs-ppt
git pull --ff-only
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate wmh-hcy
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
python configure_workstation.py --data-root /data/usersdir/linzhenzong --atlas-version v2
python run_summary.py --config config.v2.local.json
```

配置写入独立的`config.v2.local.json`，保留旧配置及旧结果。如果v2配置已经生成，只需直接运行最后一条命令。每次统计新建输出目录；新版汇总要重新读取掩膜，不能把尾状核与壳核人数相加得到基底节人数。

默认生成六页“脑图标注”表，分别列出ROI累及比例、mRS 2–6人数、有效mRS人数及残疾/死亡率。回传当次`share.zip`即可补全PPT；仅截图时请包含六页脑图标注、影像排除与图谱覆盖检查。

### v2脑区定义

| 脑区 | 定义 |
| --- | --- |
| 基底节 | 原MNI尾状核、壳核与Harvard–Oxford双侧苍白球的体素并集 |
| 内囊 | JHU双侧前肢、后肢及豆状核后部 |
| 放射冠 | JHU双侧前部、上部及后部 |
| 胼胝体 | JHU膝部、体部及压部 |
| 半卵圆中心 | 操作性近似：Harvard–Oxford双侧大脑白质中，MNI z高于该图谱侧脑室最高体素中心（z > 28 mm），排除JHU胼胝体 |

其余8区（额叶、顶叶、颞叶、枕叶、岛叶、丘脑、小脑、脑干）与v1逐体素一致。MNI脑叶分区是宏观结构标签，不能称为只含灰质的组织掩膜。半卵圆中心不是JHU原生标签，上述定义不是经过独立验证的标准分割；空间依据与完整构建记录见`atlas/v2/provenance.json`。

统计使用13个独立二值通道，允许不同ROI体素和患者重叠。全图谱覆盖率按所有通道的并集计算，区域比例不可相加。`display_MNI_HO_JHU_2mm.nii.gz`仅用于画图，其覆盖优先级不会改变统计掩膜。v1默认配置及资源保留，显式`--atlas-version v2`才选择新图谱。

### 新增分享表

- `roi_outcome_rates.csv`：每个ROI、每次随访的2–5、2–6、2–3三种结局定义，均附分子、有效分母、总N与mRS 6人数；PPT选择2–6。
- `primary_recovery_distributions.csv`及敏感性对应表：0–1恢复良好、2–3功能残疾、4–6功能丧失，比例分母为当次有效mRS记录。
- `ppt_data.json`保留原始mRS计数、规则、图谱指纹及新增聚合表。旧`state3`的0–2/3–5/6定义不改名混用。

上述类别均为本次展示约定。缺失、非法或冲突记录不进入该次结局比例分母，不作插补；在原始核对表中保留。

## 工作站首次运行（v1原流程）

在已有 `wmh-hcy` 环境的工作站执行：

```bash
cd /data/usersdir/linzhenzong
git clone https://github.com/zhenzonglin/stroke-lesion-mrs-ppt.git
cd stroke-lesion-mrs-ppt

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate wmh-hcy
export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8

python configure_workstation.py --data-root /data/usersdir/linzhenzong
python run_summary.py --config config.local.json
```

配置命令按原项目结构生成以下本地路径：

| 输入 | 相对 `--data-root` 的位置 |
| --- | --- |
| 急性梗死掩膜 | `Substain/C3_LesionMNI2mm` |
| 3、12个月mRS | `WMH/sas/c3_bl_v1y_210714.sas7bdat` |
| 24、36、48、60个月mRS | `WMH/sas/cnsr3_lzz_20260130.sas7bdat` |

`config.local.json` 被Git忽略。配置命令拒绝覆盖已有文件；修改路径或阈值时编辑该文件。也可通过 `--lesion-root`、`--early-sas`、`--late-sas` 显式指定其他输入路径。

使用Python 3.11或更新版本。运行沿用现有环境，不自动安装依赖；新建独立环境时可执行 `python -m pip install -r requirements.txt`。图谱已随代码附带，无需另行下载。

## 数据规则

- 对早期与长期SAS的ID取并集，再与全部病灶掩膜匹配。保留仅有某一阶段随访者，不套用既往3个月独立队列的限制。
- 所有时间点仅以该次 **mRS=6** 定义死亡，不读取独立死亡字段，不填补缺失。3、12个月源字段若没有6分，零计数不能解释为没有真实死亡。
- 既往6分后出现0–5分，仅将后一次访视标为 `conflict`。保留固定人数及原始频数，不作为死亡后恢复。
- ROI主定义为任意正重叠；同时导出重叠至少0.1 mL的敏感性版本。所有ROI使用相同有效匹配分母，一个患者可属于多个ROI。
- 图谱是MNI结构图谱9区与Harvard–Oxford脑干的10区展示组合，双侧合并。未覆盖病灶比例单独报告。它不是新的已验证标准图谱，也不宣称完整脑覆盖。
- 只接受与图谱网格相容的3D二值掩膜；允许保留物理坐标的轴重排，不配准或插值。缺少空间单位但网格完全相容时按配置解释为mm，并单列人数。
- 未知文件命名、重复ID或扫描不完整会阻断统计。逐例影像格式/空间问题有明确排除记录，不把空掩膜当作零病灶负担。

具体规则与状态解释见 [运行说明](运行说明.txt)。

## 输出和后续查看

每次运行新建 `outputs/<UTC时间戳>/`，不覆盖原数据或旧结果。

- `share.zip`：推荐回传文件，仅含群体汇总、来源、图谱指纹和检查记录。
- `share/ppt_data.json`：总体及各ROI的精确频数、相邻转移和敏感性结果。
- `share/*.csv`：ROI、mRS、三状态与变化分类表。
- `share/screenshots.txt`、`share/page_XX.txt`：适合逐页截图。
- `private/`：含患者ID、输入路径和逐人派生状态，仅留工作站，不进入 `share.zip`。

终端会打印当次输出路径。重新查看已有结果时，使用实际时间戳替换示例：

```bash
python run_summary.py --review outputs/<运行时间戳>/share/ppt_data.json --page 1
```

此命令只读已有汇总，不重新扫描数据。再次统计时直接运行 `python run_summary.py --config config.local.json`。

## 验证

```bash
python -m pip install -r requirements-test.txt
python -m pytest -q
```

测试使用合成NIfTI和模拟SAS读取器，覆盖ID连接、死亡定义、ROI重叠、人数守恒、输入不变和分享包范围。软件测试不代表真实患者结果或配准质量已经验证。

## 图谱来源和许可

图谱来源、变换记录、SHA256与原始文件见 `atlas/provenance.json` 和 `atlas/sources/`。构建脚本为 `build_display_atlas.py`。第三方图谱继续适用原许可，不能把整个仓库视为不受限制的商业许可。

- [FSL标准图谱说明](https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html)
- [FSL及图谱许可](https://fsl.fmrib.ox.ac.uk/fsl/docs/license.html)
- [Harvard–Oxford CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
- [随附来源与修改说明](atlas/来源与许可.txt)
