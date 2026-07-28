# Changelog

## 5.8.0 — 28 July 2026
- Fixed: consolidated lineages back under `/busco_downloads`. 5.7.2 put them in
  `/data/busco_downloads`, leaving two disjoint trees — `bacteria_odb10` and
  `file_versions.tsv` lived at one path while the spirochaete sets lived at the other.
  The task pointed at the latter, so only spirochaetes resolved and `--auto-lineage-prok`
  could not work at all (no manifest at that path).
- Added: the complete OrthoDB v10 prokaryote lineage set — 99 datasets, 1.5 GB compressed
- Added: odb10 placement files (18), enabling fully offline `--auto-lineage-prok`
- Changed: lineages are fetched at build time by `fetch_busco_data.sh` with md5
  verification, rather than shipped through the build context. Keeps the repo small and
  the build reproducible from git.
- Note: the fetch pins **odb10**. Upstream's live manifest now also advertises ~740
  OrthoDB v12 datasets, which BUSCO 5.7.x cannot consume; `busco --download prokaryota`
  resolves to odb12 and fails on md5. The script resolves odb10 entries explicitly.
- Note: the "S3 redirect issue" from 5.7.2 is `busco-data.ezlab.org` 301-ing to
  `busco-data.s3.amazonaws.com`. Fetchers that do not follow redirects receive an HTML
  error page. The script targets the S3 origin directly and uses `wget`, which follows
  redirects by default (the base image has no `curl`).
- Size: ~1.75 GB compressed, ~10 GB unpacked (the lineage layer alone is 1.51 GB
  compressed). Tasks using this image need `bootDiskSizeGb` above Cromwell's 10 GB
  default or the pull fails before the task starts.
- Digest: sha256:dd2c2cf6b9e225ff9b1ab5fcbefc98837fef9865549d2313b2be8cd855ecd894
  - Verified against the registry, not just locally. Pin with:
    `aarvani1/busco-prokaryota@sha256:dd2c2cf6b9e225ff9b1ab5fcbefc98837fef9865549d2313b2be8cd855ecd894`

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