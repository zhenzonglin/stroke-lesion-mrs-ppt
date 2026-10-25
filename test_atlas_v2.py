"""V2 atlas construction: unions, index conventions and preservation of legacy regions."""
from pathlib import Path
import json

import nibabel as nib
import numpy as np
import pytest

from build_display_atlas import V2_REGIONS, label_values, make_v2_masks, sha


ROOT = Path(__file__).parent


def test_xml_index_conventions_are_atlas_type_specific(tmp_path: Path) -> None:
    label = tmp_path / "label.xml"
    label.write_text('<atlas><header><type>Label</type></header><data><label index="17">IC</label></data></atlas>')
    prob = tmp_path / "prob.xml"
    prob.write_text('<atlas><header><type>Probabilistic</type></header><data><label index="6">Pallidum</label></data></atlas>')
    assert label_values(label, "Label")["IC"] == 17
    assert label_values(prob, "Probabilistic")["Pallidum"] == 7
    with pytest.raises(ValueError, match="Unexpected atlas type"):
        label_values(label, "Probabilistic")


def test_frozen_v2_channels_keep_union_and_overlap_contract() -> None:
    folder = ROOT / "atlas/v2"
    image = nib.load(folder / "roi_masks_MNI_HO_JHU_2mm.nii.gz")
    masks = np.asanyarray(image.dataobj)
    old = np.asanyarray(nib.load(ROOT / "atlas/display_MNI9_HObrainstem_2mm.nii.gz").dataobj)
    ho = np.asanyarray(nib.load(folder / "sources/HOsub.nii.gz").dataobj)
    jhu = np.asanyarray(nib.load(folder / "sources/JHUlabels.nii.gz").dataobj)
    assert masks.shape == old.shape + (13,)
    assert np.isin(masks, (0, 1)).all()
    for key in (2, 3, 4, 5, 6, 8, 9, 10):
        assert np.array_equal(masks[..., key-1] > 0, old == key)
    assert np.array_equal(masks[..., 0] > 0, np.isin(old, (1, 7)) | np.isin(ho, (7, 18)))
    assert np.array_equal(masks[..., 6] > 0, np.isin(jhu, (17, 18, 19, 20, 21, 22)))
    assert np.array_equal(masks[..., 10] > 0, np.isin(jhu, (23, 24, 25, 26, 27, 28)))
    assert np.array_equal(masks[..., 12] > 0, np.isin(jhu, (3, 4, 5)))
    # Cross-atlas overlaps are retained, rather than erased by display priorities.
    assert np.any(masks[..., 0] & masks[..., 6])
    assert int(masks.sum()) > int(np.any(masks, axis=-1).sum())


def test_centrum_semiovale_is_explicit_white_matter_region_not_superior_cr() -> None:
    folder = ROOT / "atlas/v2"
    image = nib.load(folder / "roi_masks_MNI_HO_JHU_2mm.nii.gz")
    masks = np.asanyarray(image.dataobj)
    ho = np.asanyarray(nib.load(folder / "sources/HOsub.nii.gz").dataobj)
    jhu = np.asanyarray(nib.load(folder / "sources/JHUlabels.nii.gz").dataobj)
    cso = masks[..., 11] > 0
    coords = np.indices(image.shape[:3])
    world_z = sum(image.affine[2, a] * coords[a] for a in range(3)) + image.affine[2, 3]
    roof = float(world_z[np.isin(ho, (3, 14))].max())
    expected = np.isin(ho, (1, 12)) & (world_z > roof) & ~np.isin(jhu, (3, 4, 5))
    assert roof == 28.0
    assert np.array_equal(cso, expected)
    assert not np.any(cso & ~np.isin(ho, (1, 12)))
    assert not np.any(cso & np.isin(jhu, (3, 4, 5)))
    assert not np.array_equal(cso, np.isin(jhu, (25, 26)))
    assert np.any(cso & (masks[..., 10] > 0))


def test_manifest_hashes_counts_and_display_only_separation() -> None:
    folder = ROOT / "atlas/v2"
    manifest = json.loads((folder / "roi_manifest.json").read_text(encoding="utf-8"))
    provenance = json.loads((folder / "provenance.json").read_text(encoding="utf-8"))
    masks = np.asanyarray(nib.load(folder / manifest["statistical_atlas"]).dataobj)
    display = np.asanyarray(nib.load(folder / manifest["display_only_atlas"]).dataobj)
    assert manifest["atlas_format"] == "independent_binary_channels"
    assert provenance["atlas_sha256"] == sha(folder / manifest["statistical_atlas"])
    for i, row in enumerate(manifest["channels"]):
        assert row["label_id"] == i + 1
        assert row["channel_index"] == i
        assert row["roi_name"] == V2_REGIONS[i][1]
        assert row["mask_sha256"] == sha(folder / row["mask_path"])
        individual = np.asanyarray(nib.load(folder / row["mask_path"]).dataobj)
        assert np.array_equal(individual, masks[..., i])
        assert row["voxel_count"] == int(individual.sum())
    assert np.array_equal(display > 0, np.any(masks, axis=-1))
    assert not np.array_equal(display == 1, masks[..., 0] > 0)
    overlap = np.array(provenance["pairwise_overlap_voxels"])
    assert np.array_equal(overlap, overlap.T)
    assert int(np.any(masks, axis=-1).sum()) == provenance["union_voxels"]


def test_missing_ventricles_cannot_silently_create_cso() -> None:
    folder = ROOT / "atlas/v2/sources"
    ho_names = label_values(folder / "HOsub.xml", "Probabilistic")
    jhu_names = label_values(folder / "JHUlabels.xml", "Label")
    with pytest.raises(ValueError, match="No lateral ventricle"):
        make_v2_masks(np.zeros((2, 2, 2)), np.zeros((2, 2, 2)), np.zeros((2, 2, 2)),
                      np.eye(4), ho_names, jhu_names)
