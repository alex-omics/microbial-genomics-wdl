# Changelog

## 5.7.2 — 23 June 2026
- Fixed: corrected download path from /busco_downloads to /data/busco_downloads
- Changed: switched from busco --download to COPY of pre-downloaded tarballs (S3 redirect issue)
- Added: spirochaetales_odb10, spirochaetia_odb10
- Digest: sha256:7f6bc636fe2ddeded8e80ae41ef85241f502362f762c426e9792ec5c4a4748b6

## 5.7.1 - 16 June 2026
- Initial build
- Base: staphb/busco:5.7.1-prok-bacteria_odb10_2024_01_08
- Added: all prokaryote lineages via `busco --download prokaryota`
- Digest: sha256:299d5d6d7969f0444bd4fafa8c1a87e3e08c9a45cb0eeee243166680ad8c9326