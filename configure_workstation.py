#!/usr/bin/env python3
"""生成不纳入Git的工作站配置，不读取或修改患者数据。"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def create_config(data_root: Path, output: Path, lesion_root: Path | None = None,
                  early_sas: Path | None = None, late_sas: Path | None = None) -> dict:
    base = Path(__file__).resolve().parent
    cfg = json.loads((base / "config.example.json").read_text(encoding="utf-8"))
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
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent/"config.local.json")
    parser.add_argument("--lesion-root", type=Path)
    parser.add_argument("--early-sas", type=Path)
    parser.add_argument("--late-sas", type=Path)
    args = parser.parse_args()
    try:
        cfg = create_config(args.data_root, args.output, args.lesion_root, args.early_sas, args.late_sas)
    except FileExistsError:
        parser.error("配置已存在，未覆盖。直接运行run_summary.py，或自行编辑现有config.local.json。")
    print("已生成本地配置:", args.output.resolve())
    print("病灶目录:", cfg["lesion_root"])
    print("该配置不会自动上传Git；本命令未读取患者数据。")


if __name__ == "__main__":
    main()
