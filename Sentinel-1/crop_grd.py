#!/usr/bin/env python3
"""Crop the S1 GRD COG measurements (CDSE S3) in SAR pixel space for
Sentinel1GrdRasterSourceProviderTest, and fetch the matching COG annotation XMLs.

Requires: GDAL Python bindings, env AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY.
Usage:    python crop_s1_testdata.py <path/to/zeebrugge_2020_06_06.SAFE>
"""
import math
import os
import sys

from osgeo import gdal

gdal.UseExceptions()
for k, v in {
    "AWS_S3_ENDPOINT": "eodata.dataspace.copernicus.eu",
    "AWS_VIRTUAL_HOSTING": "FALSE",
    "AWS_HTTPS": "YES",
    "AWS_SECRET_ACCESS_KEY": "hpI6ysKtrtmHIjsNCxhbcWtot4K6QDQ88r6riz8N",
    "AWS_ACCESS_KEY_ID": "M8KKEC39PG127ZB3EULO"
}.items():
    gdal.SetConfigOption(k, os.environ.get(k, v))

SAFE = ("/vsis3/eodata/Sentinel-1/SAR/IW_GRDH_1S-COG/2020/06/06/"
        "S1B_IW_GRDH_1SDV_20200606T060612_20200606T060637_021909_029944_1FC2_COG.SAFE")
TAG = "20200606t060612-20200606t060637-021909-029944"    # remote (COG product)
LTAG = "20200606t060615-20200606t060640-021909-029944"   # local names used by test/STAC JSON
POLS = [("vv", "001"), ("vh", "002")]

# Test AOI (lon/lat) + margin for terrain relief / resampling at the edges.
AOI = (3.1, 51.27, 3.3, 51.37)
MARGIN_DEG = 0.05

CREATION_OPTIONS = [
    "TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9",
    "DISCARD_LSB=2",  # lossy; GeoTrellis-readable (unlike LERC/ZSTD)
]


def pixel_window(ds):
    """Bounding SAR pixel/line window of the AOI, using the GCPs (thin-plate spline)."""
    if ds.GetGCPCount() == 0:
        raise RuntimeError("Source has no GCPs")
    tr = gdal.Transformer(ds, None, ["METHOD=GCP_TPS"])

    xmin, ymin, xmax, ymax = AOI
    xmin, ymin, xmax, ymax = xmin - MARGIN_DEG, ymin - MARGIN_DEG, xmax + MARGIN_DEG, ymax + MARGIN_DEG

    # Sample the full AOI border, not just corners: the geo->SAR mapping is non-linear.
    n = 20
    pts = []
    for i in range(n + 1):
        f = i / n
        x = xmin + f * (xmax - xmin)
        y = ymin + f * (ymax - ymin)
        pts += [(x, ymin), (x, ymax), (xmin, y), (xmax, y)]

    cols, rows = [], []
    for lon, lat in pts:
        ok, (c, r, _) = tr.TransformPoint(1, lon, lat, 0.0)  # 1 = georef -> pixel/line
        if not ok:
            raise RuntimeError(f"GCP transform failed for ({lon}, {lat})")
        cols.append(c)
        rows.append(r)

    c0 = max(0, math.floor(min(cols)))
    r0 = max(0, math.floor(min(rows)))
    c1 = min(ds.RasterXSize, math.ceil(max(cols)))
    r1 = min(ds.RasterYSize, math.ceil(max(rows)))
    if c1 <= c0 or r1 <= r0:
        raise RuntimeError(f"AOI does not intersect the scene: cols {c0}-{c1}, rows {r0}-{r1}")
    return c0, r0, c1 - c0, r1 - r0


def copy_vsi(src, dst):
    f = gdal.VSIFOpenL(src, "rb")
    if f is None:
        raise RuntimeError(f"Cannot open {src}")
    try:
        gdal.VSIFSeekL(f, 0, 2)
        size = gdal.VSIFTellL(f)
        gdal.VSIFSeekL(f, 0, 0)
        data = gdal.VSIFReadL(1, size, f)
    finally:
        gdal.VSIFCloseL(f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "wb") as out:
        out.write(data)
    print(f"copied {src} -> {dst}")


def main(out_safe):
    # VV and VH share the same SAR grid, so one window serves both.
    vv_src = f"{SAFE}/measurement/s1b-iw-grd-vv-{TAG}-001-cog.tiff"
    src_ds = gdal.Open(vv_src)
    xoff, yoff, w, h = pixel_window(src_ds)
    full = (src_ds.RasterXSize, src_ds.RasterYSize)
    src_ds = None
    print(f"full scene {full[0]}x{full[1]}, srcwin: xoff={xoff} yoff={yoff} w={w} h={h}")

    os.makedirs(os.path.join(out_safe, "measurement"), exist_ok=True)
    for pol, idx in POLS:
        src = f"{SAFE}/measurement/s1b-iw-grd-{pol}-{TAG}-{idx}-cog.tiff"
        dst = os.path.join(out_safe, "measurement", f"s1b-iw-grd-{pol}-{LTAG}-{idx}.tiff")
        gdal.Translate(
            dst, src,
            format="GTiff",
            srcWin=[xoff, yoff, w, h],
            creationOptions=CREATION_OPTIONS,
            metadataOptions=[f"S1_CROP_COL_OFFSET={xoff}", f"S1_CROP_ROW_OFFSET={yoff}"],
        )
        print(f"wrote {dst}")

    # Annotations must come from the same (COG) product, so line timing matches the offsets.
    for pol, idx in POLS:
        for sub, prefix in [("annotation", ""),
                            ("annotation/calibration", "calibration-"),
                            ("annotation/calibration", "noise-")]:
            src = f"{SAFE}/{sub}/{prefix}s1b-iw-grd-{pol}-{TAG}-{idx}-cog.xml"
            dst = os.path.join(out_safe, sub, f"{prefix}s1b-iw-grd-{pol}-{LTAG}-{idx}.xml")
            copy_vsi(src, dst)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])