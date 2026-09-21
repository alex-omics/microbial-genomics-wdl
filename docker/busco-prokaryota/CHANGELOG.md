# Changelog

## 5.8.0 — 28 July 2026
- Fixed: lineages consolidated under `/busco_downloads`. 5.7.2 split them across two trees, so
  only the spirochaete datasets resolved and `--auto-lineage-prok` could not run.
- Added: the complete OrthoDB v10 prokaryote lineage set (99 datasets, 1.5 GB compressed).
- Added: 18 odb10 placement files, enabling fully offline `--auto-lineage-prok`.
- Changed: lineages are fetched at build time by `fetch_busco_data.sh` with md5 verification,
  rather than shipped through the build context.
- Note: the fetch pins **odb10**. The upstream manifest also lists ~740 OrthoDB v12 datasets,
  which BUSCO 5.7.x cannot use; `busco --download prokaryota` resolves to odb12 and fails.
- Note: `busco-data.ezlab.org` redirects to `busco-data.s3.amazonaws.com`. The script fetches
  from the S3 origin with `wget`, which follows redirects (the base image has no `curl`).
- Size: ~1.75 GB compressed, ~10 GB unpacked. Tasks using this image need `bootDiskSizeGb`
  above Cromwell's 10 GB default.
- Digest: sha256:dd2c2cf6b9e225ff9b1ab5fcbefc98837fef9865549d2313b2be8cd855ecd894
  - Pin with `aarvani1/busco-prokaryota@sha256:dd2c2cf6b9e225ff9b1ab5fcbefc98837fef9865549d2313b2be8cd855ecd894`

## 5.7.2 — 23 June 2026
- Fixed: download path corrected from `/busco_downloads` to `/data/busco_downloads`.
- Changed: pre-downloaded tarballs are copied in rather than fetched with `busco --download`.
- Added: `spirochaetales_odb10`, `spirochaetia_odb10`.
- Digest: sha256:7f6bc636fe2ddeded8e80ae41ef85241f502362f762c426e9792ec5c4a4748b6

## 5.7.1 — 16 June 2026
- Initial build.
- Base: `staphb/busco:5.7.1-prok-bacteria_odb10_2024_01_08`
- Added: all prokaryote lineages via `busco --download prokaryota`.
- Digest: sha256:299d5d6d7969f0444bd4fafa8c1a87e3e08c9a45cb0eeee243166680ad8c9326
