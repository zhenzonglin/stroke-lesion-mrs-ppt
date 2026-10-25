# Stroke lesion ROI and longitudinal mRS summaries

只读提取急性梗死病灶的ROI累及比例，以及全体可匹配患者和各ROI人群的长期mRS分布、相邻随访转移。输出用于后续交互式PowerPoint制作；本仓库目前不包含PPT生成器，也不包含患者数据或真实统计结果。

## 工作站首次运行

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
