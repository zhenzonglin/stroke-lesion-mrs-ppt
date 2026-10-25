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
V2_FILES = {**FILES,
    "JHUlabels.nii.gz": "data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz",
    "JHUlabels.xml": "data/atlases/JHU-labels.xml",
}

V2_REGIONS = [
    (1, "Basal Ganglia", "基底节", "basal_ganglia"),
    (2, "Cerebellum", "小脑", "cerebellum"),
    (3, "Frontal Lobe", "额叶", "frontal_lobe"),
    (4, "Insula", "岛叶", "insula"),
    (5, "Occipital Lobe", "枕叶", "occipital_lobe"),
    (6, "Parietal Lobe", "顶叶", "parietal_lobe"),
    (7, "Internal Capsule", "内囊", "internal_capsule"),
    (8, "Temporal Lobe", "颞叶", "temporal_lobe"),
    (9, "Thalamus", "丘脑", "thalamus"),
    (10, "Brainstem", "脑干", "brainstem"),
    (11, "Corona Radiata", "放射冠", "corona_radiata"),
    (12, "Centrum Semiovale", "半卵圆中心（操作性定义）", "centrum_semiovale"),
    (13, "Corpus Callosum", "胼胝体", "corpus_callosum"),
]


def label_values(xml_path: Path, expected_type: str) -> dict[str, int]:
    """FSL Label图直接用index；Probabilistic的maxprob图才加1。"""
    root = ET.parse(xml_path)
    atlas_type = root.findtext("header/type")
    if atlas_type != expected_type:
        raise ValueError(f"Unexpected atlas type: {atlas_type}, expected {expected_type}")
    offset = 1 if atlas_type == "Probabilistic" else 0
    return {e.text.strip(): int(e.attrib["index"]) + offset for e in root.findall(".//label")}


def make_v2_masks(mni: np.ndarray, ho: np.ndarray, jhu: np.ndarray,
                  affine: np.ndarray, ho_labels: dict[str, int],
                  jhu_labels: dict[str, int]) -> tuple[np.ndarray, dict]:
    """独立ROI可重叠。基底节内部作并集，其余原ROI不被新白质标签覆盖。"""
    if mni.shape != ho.shape or mni.shape != jhu.shape or mni.ndim != 3:
        raise ValueError("Atlas arrays must have identical 3D shape")
    brainstem = ho == ho_labels["Brain-Stem"]
    v1 = mni.copy()
    v1[brainstem] = 10
    pallidum_ids = [ho_labels[side + " Pallidum"] for side in ("Left", "Right")]
    pallidum = np.isin(ho, pallidum_ids)
    ic_names = [part + " of internal capsule " + side for part in
                ("Anterior limb", "Posterior limb", "Retrolenticular part") for side in ("R", "L")]
    cr_names = [part + " corona radiata " + side for part in
                ("Anterior", "Superior", "Posterior") for side in ("R", "L")]
    cc_names = [part + " of corpus callosum" for part in ("Genu", "Body", "Splenium")]
    ic_ids = [jhu_labels[n] for n in ic_names]
    cr_ids = [jhu_labels[n] for n in cr_names]
    cc_ids = [jhu_labels[n] for n in cc_names]
    ic, cr, cc = (np.isin(jhu, ids) for ids in (ic_ids, cr_ids, cc_ids))
    bg = np.isin(v1, (1, 7)) | pallidum
    wm_ids = [ho_labels[side + " Cerebral White Matter"] for side in ("Left", "Right")]
    vent_ids = [ho_labels[side + " Lateral Ventricle"] for side in ("Left", "Right")]
    wm, vent = np.isin(ho, wm_ids), np.isin(ho, vent_ids)
    if not np.any(vent):
        raise ValueError("No lateral ventricle voxels; cannot define centrum semiovale boundary")
    # 用物理z坐标确定脑室顶上方，不把数组轴号误作MNI坐标。
    ijk = np.indices(mni.shape, dtype=np.float64)
    world_z = sum(affine[2, axis] * ijk[axis] for axis in range(3)) + affine[2, 3]
    vent_roof_z = float(world_z[vent].max())
    cso = wm & (world_z > vent_roof_z + 1e-8) & ~cc
    masks = np.stack([bg if key == 1 else ic if key == 7 else cr if key == 11
                      else cso if key == 12 else cc if key == 13 else v1 == key
                      for key, _, _, _ in V2_REGIONS], axis=-1).astype(np.uint8)
    metadata = {
        "basal_ganglia": {"MNI_v1_label_values": [1, 7], "HO_pallidum_label_values": pallidum_ids,
                           "definition": "Union of v1 Caudate and Putamen plus bilateral HO Pallidum"},
        "internal_capsule": {"JHU_label_values": ic_ids, "labels": ic_names},
        "corona_radiata": {"JHU_label_values": cr_ids, "labels": cr_names},
        "corpus_callosum": {"JHU_label_values": cc_ids, "labels": cc_names},
        "centrum_semiovale": {
            "is_operational_approximation": True,
            "definition": "HO bilateral cerebral white matter above the highest HO lateral-ventricle voxel centre, excluding JHU corpus callosum. No erosion or additional atlas labels.",
            "HO_white_matter_label_values": wm_ids, "HO_lateral_ventricle_label_values": vent_ids,
            "MNI_z_strictly_greater_than_mm": vent_roof_z,
            "excluded_JHU_label_values": cc_ids,
            "not_equivalent_to": "JHU Superior Corona Radiata",
            "overlap_with_corona_radiata_voxels": int(np.count_nonzero(cso & cr)),
            "overlap_with_corpus_callosum_voxels": int(np.count_nonzero(cso & cc)),
            "outside_HO_cerebral_white_matter_voxels": int(np.count_nonzero(cso & ~wm)),
            "literature_support_is_anatomical_location_only": True,
            "source_reference": "https://doi.org/10.3389/fneur.2022.888511",
            "limitation": "Custom approximate region, not a validated CSO atlas or a reproduction of the cited study's segmentation. Cortical U-fibres are not specifically removed.",
        },
        "unchanged_v1_label_values": [2, 3, 4, 5, 6, 8, 9, 10],
    }
    return masks, metadata


def _write_image(data: np.ndarray, reference: nib.Nifti1Image, path: Path) -> None:
    header = reference.header.copy()
    header.set_data_dtype(np.uint8)
    header.set_xyzt_units("mm")
    image = nib.Nifti1Image(data.astype(np.uint8), reference.affine, header)
    image.set_qform(reference.affine, code=4)
    image.set_sform(reference.affine, code=4)
    nib.save(image, path)


def build_v2(package: Path, out: Path) -> dict:
    """构建独立13通道统计图、各ROI二值图及单标签可视化图；不读患者数据。"""
    if sha(package) != PACKAGE_SHA256:
        raise ValueError("Official package SHA256 mismatch")
    if (out / "display_MNI9_HObrainstem_2mm.nii.gz").exists():
        raise ValueError("Choose a separate v2 output directory; v1 atlas must remain intact")
    out.mkdir(parents=True, exist_ok=True)
    sources, mask_dir = out / "sources", out / "roi_masks"
    sources.mkdir(exist_ok=True)
    mask_dir.mkdir(exist_ok=True)
    with tarfile.open(package) as archive:
        for name, member in V2_FILES.items():
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError("Missing package member: " + member)
            (sources / name).write_bytes(stream.read())
    mni, ho, jhu = (nib.load(sources / name) for name in ("MNI9.nii.gz", "HOsub.nii.gz", "JHUlabels.nii.gz"))
    for image in (ho, jhu):
        if image.shape != mni.shape or not np.allclose(image.affine, mni.affine, atol=1e-4, rtol=0):
            raise ValueError("Atlas grids differ; no interpolation is permitted")
    masks, definitions = make_v2_masks(np.asanyarray(mni.dataobj), np.asanyarray(ho.dataobj),
        np.asanyarray(jhu.dataobj), mni.affine,
        label_values(sources / "HOsub.xml", "Probabilistic"),
        label_values(sources / "JHUlabels.xml", "Label"))
    if any(not masks[..., i].any() for i in range(masks.shape[-1])):
        raise ValueError("One or more derived ROI channels are empty")
    atlas_path = out / "roi_masks_MNI_HO_JHU_2mm.nii.gz"
    _write_image(masks, mni, atlas_path)
    records = []
    lines = ["label_id\troi_name\troi_name_zh\tsource"]
    for channel, name, chinese, slug in V2_REGIONS:
        file = mask_dir / (slug + ".nii.gz")
        _write_image(masks[..., channel-1], mni, file)
        source = ("MNI Caudate + Putamen union HO Pallidum" if channel == 1 else
                  "JHU ICBM-DTI-81 labels " + ",".join(map(str, definitions[slug]["JHU_label_values"])) if channel in (7, 11, 13) else
                  "Operational HO white matter above lateral ventricle roof excluding JHU CC" if channel == 12 else
                  "Unchanged v1 label " + str(channel))
        lines.append(f"{channel}\t{name}\t{chinese}\t{source}")
        records.append({"label_id": channel, "channel_index": channel-1, "roi_id": f"roi_{channel:04d}",
                        "roi_name": name, "roi_name_zh": chinese, "slug": slug,
                        "mask_path": file.relative_to(out).as_posix(), "mask_sha256": sha(file),
                        "voxel_count": int(masks[..., channel-1].sum()), "definition": source})
    labels = out / "labels.tsv"
    labels.write_text("\n".join(lines) + "\n", encoding="utf-8")
    # 单标签图仅方便PPT显示；后绘制的特异小ROI优先，不可作为统计输入。
    priority = [2, 3, 4, 5, 6, 8, 9, 10, 12, 11, 13, 1, 7]
    display = np.zeros(mni.shape, dtype=np.uint8)
    for key in priority:
        display[masks[..., key-1] > 0] = key
    display_path = out / "display_MNI_HO_JHU_2mm.nii.gz"
    _write_image(display, mni, display_path)
    flat = masks.reshape(-1, masks.shape[-1]).astype(np.int64)
    pairwise = (flat.T @ flat).tolist()
    union_voxels = int(np.any(masks, axis=-1).sum())
    manifest = {"version": "20261025_v2", "atlas_format": "independent_binary_channels",
                "statistical_atlas": atlas_path.name, "labels": labels.name,
                "spatial_shape": list(mni.shape), "affine": mni.affine.tolist(),
                "channels": records, "roi_voxel_overlap_allowed": True,
                "coverage_rule": "OR union of all binary channels; overlapping voxels counted once",
                "display_only_atlas": display_path.name, "display_priority_low_to_high": priority,
                "display_warning": "The priority label image is for visualization only. Use the 4D binary-channel atlas for statistics."}
    (out / "roi_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    provenance = {"name": "Custom independent MNI HO JHU display regions", "version": "20261025_v2",
        "is_standard_unmodified_atlas": False, "space": "FSL MNI152 2mm", "bilateral_regions": True,
        "atlas_format": "independent_binary_channels", "channel_order": [r["roi_name"] for r in records],
        "shape": list(masks.shape), "affine": mni.affine.tolist(), "definitions": definitions,
        "construction": "No resampling. Preserve the eight unchanged v1 regions as binary channels. Merge caudate, putamen and pallidum within one channel. Add independent white-matter channels without overwriting existing channels.",
        "pairwise_overlap_voxels": pairwise, "union_voxels": union_voxels,
        "sum_channel_voxels": int(masks.sum()), "atlas_sha256": sha(atlas_path),
        "labels_sha256": sha(labels), "display_only_atlas_sha256": sha(display_path),
        "package_url": PACKAGE_URL, "package_sha256": PACKAGE_SHA256,
        "source_files": {name: {"package_member": member, "sha256": sha(sources/name)} for name, member in V2_FILES.items()},
        "source_docs": ["https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html", "https://fsl.fmrib.ox.ac.uk/fsl/docs/other/atlas_xml.html"],
        "source_label_index_rules": {"MNI_HO_maxprob": "XML index + 1", "JHU_Label": "XML index directly equals voxel value"},
        "license": "Source atlas terms apply. JHU atlas data may not be used commercially without permission under the supplied FSL atlas terms. Harvard-Oxford is CC BY-SA 4.0. See ../来源与许可.txt, ../FSL_license_original.html and ../CC_BY_SA_4_legalcode.txt.",
        "limitations": ["Custom regions, not a new validated whole-brain atlas.", "Centrum semiovale is an explicitly operational approximate region, not a JHU label.", "White-matter and lobar regions may overlap; counts and proportions must not be added.", "The single-label display image must not be used for ROI statistics.", "MNI grid compatibility does not establish patient anatomical registration quality."]}
    (out / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    notice = """V2 atlas sources and adaptation notice

Frozen source: FSL fsl-data_atlases 2103.0-0, official package URL and SHA256 in provenance.json.
MNI Structural Atlas: Jack Lancaster, UTHSCSA; original labels from the McConnell Brain Imaging Centre, MNI.
Harvard-Oxford Subcortical Atlas: Harvard Center for Morphometric Analysis and contributors credited by FSL.
JHU ICBM-DTI-81 White-Matter Labels: Susumu Mori and colleagues, Johns Hopkins University; ICBM DTI workgroup.
References: Mori et al. MRI Atlas of Human White Matter (2005); FSL dataset documentation.
https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html
https://fsl.fmrib.ox.ac.uk/fsl/docs/license.html

The source atlas license terms remain applicable. Harvard-Oxford is CC BY-SA 4.0.
JHU atlas data may not be used commercially without permission from its copyright holder
under the FSL-distributed atlas terms. This bundle does not grant such permission.
The source FSL license page and CC BY-SA 4.0 legal text are retained in the parent atlas directory.
This academic research bundle does not assert commercial-use rights or endorsement by source institutions.

Modifications: bilateral unions; caudate/putamen/pallidum merged; independent internal capsule,
corona radiata, corpus callosum and an explicitly operational centrum semiovale channel added.
Centrum semiovale is HO cerebral white matter above the highest HO lateral-ventricle voxel centre,
excluding JHU corpus callosum. It is not an official atlas label, not equivalent to superior corona
radiata, not a validated segmentation, and does not specifically remove cortical U-fibres.
Anatomical-location reference only: https://doi.org/10.3389/fneur.2022.888511
The source study's acquisition, slice selection and segmentation are not reproduced here.

Statistics must use roi_masks_MNI_HO_JHU_2mm.nii.gz with independent_binary_channels mode.
display_MNI_HO_JHU_2mm.nii.gz is a lossy priority rendering aid only.
Original v1 masks remain unchanged. No patient images or clinical records were used to construct these atlas assets.
"""
    (out / "来源与定义.txt").write_text(notice, encoding="utf-8")
    print(json.dumps({"atlas": str(atlas_path), "sha256": sha(atlas_path), "channels": len(records),
                      "CSO_roof_z_mm": definitions["centrum_semiovale"]["MNI_z_strictly_greater_than_mm"],
                      "union_voxels": union_voxels}, ensure_ascii=False))
    return provenance


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
    parser.add_argument("--out",type=Path)
    parser.add_argument("--version", choices=("v1", "v2"), default="v1")
    args=parser.parse_args()
    output = args.out or (Path(__file__).parent / "atlas" / "v2" if args.version == "v2" else Path(__file__).parent / "atlas")
    (build_v2 if args.version == "v2" else build)(args.package, output)
