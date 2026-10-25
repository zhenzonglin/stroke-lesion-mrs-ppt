"""由冻结FSL图谱生成本次展示用分区；不会读取患者数据。"""
from pathlib import Path
import argparse
import hashlib
import json
import tarfile
import xml.etree.ElementTree as ET

import nibabel as nib
import numpy as np

PACKAGE_SHA256 = "3ac89f5f67de84054e1e576c9ad5f7ef6f562cef7b2504a3fa992c079cf72f8c"
PACKAGE_URL = "https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/noarch/fsl-data_atlases-2103.0-0.tar.bz2"
FILES = {
    "MNI9.nii.gz": "data/atlases/MNI/MNI-maxprob-thr25-2mm.nii.gz",
    "HOsub.nii.gz": "data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz",
    "MNI.xml": "data/atlases/MNI.xml",
    "HOsub.xml": "data/atlases/HarvardOxford-Subcortical.xml",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(package: Path, out: Path) -> None:
    if sha(package) != PACKAGE_SHA256:
        raise ValueError("Official package SHA256 mismatch")
    out.mkdir(parents=True, exist_ok=True)
    sources = out / "sources"
    sources.mkdir(exist_ok=True)
    with tarfile.open(package) as archive:
        for name, member in FILES.items():
            (sources / name).write_bytes(archive.extractfile(member).read())
    mni, ho = nib.load(sources / "MNI9.nii.gz"), nib.load(sources / "HOsub.nii.gz")
    if mni.shape != ho.shape or not np.allclose(mni.affine, ho.affine, atol=1e-4, rtol=0):
        raise ValueError("Atlas grids differ; no resampling is performed")
    base = np.asanyarray(mni.dataobj).astype(np.uint8)
    extra = np.asanyarray(ho.dataobj)
    mni_labels = {int(e.attrib["index"]) + 1: e.text.strip()
                  for e in ET.parse(sources / "MNI.xml").findall(".//label")}
    brainstem = [int(e.attrib["index"]) + 1 for e in ET.parse(sources / "HOsub.xml").findall(".//label")
                 if e.text.strip() == "Brain-Stem"]
    assert len(brainstem) == 1
    mask = extra == brainstem[0]
    overlap = {mni_labels[k]: int(np.count_nonzero(mask & (base == k))) for k in mni_labels}
    data = base.copy()
    data[mask] = 10
    header = mni.header.copy()
    header.set_data_dtype(np.uint8)
    header.set_xyzt_units("mm")
    image = nib.Nifti1Image(data, mni.affine, header)
    image.set_qform(mni.affine, code=4)
    image.set_sform(mni.affine, code=4)
    atlas_path = out / "display_MNI9_HObrainstem_2mm.nii.gz"
    nib.save(image, atlas_path)
    names_cn = {1:"尾状核",2:"小脑",3:"额叶",4:"岛叶",5:"枕叶",6:"顶叶",7:"壳核",8:"颞叶",9:"丘脑",10:"脑干"}
    labels = {**mni_labels, 10: "Brainstem"}
    lines = ["label_id\troi_name\troi_name_zh\tsource"]
    for key, label in labels.items():
        source = "HO subcortical maxprob thr25; label8" if key == 10 else f"MNI structural maxprob thr25; label{key}; excluding HO brainstem"
        lines.append(f"{key}\t{label}\t{names_cn[key]}\t{source}")
    labels_path = out / "labels.tsv"
    labels_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    provenance = {"name":"Custom display atlas MNI9 plus HO brainstem", "version":"20261025_v1",
        "is_standard_unmodified_atlas":False, "space":"FSL MNI152 2mm", "bilateral_regions":True,
        "construction":"Keep MNI structural maxprob thr25 labels 1-9; assign all HO subcortical maxprob thr25 Brain-Stem voxels to label10 with priority over MNI9. No interpolation.",
        "overridden_MNI_voxels":overlap, "brainstem_voxels":int(mask.sum()),
        "shape":list(image.shape), "affine":image.affine.tolist(),
        "atlas_sha256":sha(atlas_path), "labels_sha256":sha(labels_path),
        "package_url":PACKAGE_URL, "package_sha256":PACKAGE_SHA256,
        "source_files":{name:sha(sources/name) for name in FILES},
        "source_docs":"https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html",
        "license":"See FSL license and Harvard-Oxford CC BY-SA 4.0; source files and build code retained.",
        "limitations":["Display-level bilateral macro regions, not a new validated atlas.",
                        "Unlabelled lesion voxels are quantified separately; not complete brain coverage.",
                        "MNI grid compatibility does not establish patient anatomical registration quality."]}
    (out / "provenance.json").write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"atlas":str(atlas_path),"sha256":sha(atlas_path),"labels":labels,"override_voxels":sum(overlap.values())},ensure_ascii=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package",type=Path)
    parser.add_argument("--out",type=Path,default=Path(__file__).parent/"atlas")
    args=parser.parse_args()
    build(args.package,args.out)
