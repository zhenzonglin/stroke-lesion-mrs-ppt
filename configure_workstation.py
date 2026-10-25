#!/usr/bin/env python3
"""生成不纳入Git的工作站配置，不读取或修改患者数据。"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def create_config(data_root: Path, output: Path, lesion_root: Path | None = None,
                  early_sas: Path | None = None, late_sas: Path | None = None,
                  atlas_version: str = "v1", display_month: int | None = None,
                  display_definition: str | None = None) -> dict:
    base = Path(__file__).resolve().parent
    cfg = json.loads((base / "config.example.json").read_text(encoding="utf-8"))
    if atlas_version not in {"v1", "v2"}:
        raise ValueError("atlas_version must be v1 or v2")
    if (display_month is None) != (display_definition is None):
        raise ValueError("display_month and display_definition must be selected together")
    if display_month is not None:
        if display_month not in (3, 12, 24, 36, 48, 60):
            raise ValueError("Unsupported display month")
        if display_definition not in ("mrs_2_5", "mrs_2_6", "mrs_2_3"):
            raise ValueError("Unsupported display outcome definition")
        cfg["display_outcome"] = {"month": display_month, "definition": display_definition}
    if atlas_version == "v2":
        cfg.update(atlas_path="atlas/v2/roi_masks_MNI_HO_JHU_2mm.nii.gz",
                   atlas_labels="atlas/v2/labels.tsv", atlas_provenance="atlas/v2/provenance.json",
                   atlas_format="independent_binary_channels")
        if display_month is None:
            cfg["display_outcome"] = {"months": [3, 12, 24, 36, 48, 60], "definition": "mrs_2_6"}
    cfg["atlas_version"] = atlas_version
    data_root = data_root.expanduser().resolve()
    cfg["lesion_root"] = str((lesion_root or data_root/"Substain/C3_LesionMNI2mm").expanduser().resolve())
    cfg["source_files"] = {
        "early": str((early_sas or data_root/"WMH/sas/c3_bl_v1y_210714.sas7bdat").expanduser().resolve()),
        "late": str((late_sas or data_root/"WMH/sas/cnsr3_lzz_20260130.sas7bdat").expanduser().resolve()),
    }
    # 默认配置在仓库根目录，图谱相对该目录。其他输出位置使用绝对图谱路径。
    if output.expanduser().resolve().parent != base:
        for key in ("atlas_path", "atlas_labels", "atlas_provenance"):
            cfg[key] = str(base/cfg[key])
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(cfg, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return cfg


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True,
                        help="包含Substain和WMH两个项目的父目录")
    parser.add_argument("--output", type=Path, help="默认v1为config.local.json，v2为config.v2.local.json")
    parser.add_argument("--atlas-version", choices=("v1", "v2"), default="v1")
    parser.add_argument("--display-month", type=int, choices=(3, 12, 24, 36, 48, 60))
    parser.add_argument("--display-definition", choices=("mrs_2_5", "mrs_2_6", "mrs_2_3"))
    parser.add_argument("--lesion-root", type=Path)
    parser.add_argument("--early-sas", type=Path)
    parser.add_argument("--late-sas", type=Path)
    args = parser.parse_args()
    output = args.output or Path(__file__).resolve().parent / ("config.v2.local.json" if args.atlas_version == "v2" else "config.local.json")
    try:
        cfg = create_config(args.data_root, output, args.lesion_root, args.early_sas, args.late_sas,
                            args.atlas_version, args.display_month, args.display_definition)
    except FileExistsError:
        parser.error("配置已存在，未覆盖。使用run_summary.py --config指定该配置，或选择新的--output。")
    except ValueError as exc:
        parser.error(str(exc))
    print("已生成本地配置:", output.resolve())
    print("图谱版本:", cfg["atlas_version"])
    print("病灶目录:", cfg["lesion_root"])
    print("该配置不会自动上传Git；本命令未读取患者数据。")


if __name__ == "__main__":
    main()
