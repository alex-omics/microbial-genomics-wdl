#!/usr/bin/env bash
#
# Populate an offline BUSCO download directory with the complete set of
# OrthoDB v10 prokaryote lineages plus the placement files required for
# --auto-lineage-prok.
#
# Why this exists, and why it pins odb10:
#   busco --download prokaryota consults the LIVE manifest, which since
#   2026-05 also advertises ~740 OrthoDB v12 datasets. BUSCO 5.7.x cannot
#   consume odb12, and pulling that tree would add tens of GB to the image.
#   We therefore resolve the 99 odb10 Prokaryota lineages explicitly and
#   fetch them by URL with md5 verification.
#
# Note on hosts: busco-data.ezlab.org issues a 301 to busco-data.s3.amazonaws.com.
#   Any fetch that does not follow redirects silently receives an HTML error
#   page instead of data. We hit the S3 origin directly to sidestep this.

set -euo pipefail

DOWNLOAD_PATH="${1:-/busco_downloads}"
S3_BASE="https://busco-data.s3.amazonaws.com"
JOBS="${FETCH_JOBS:-8}"

mkdir -p "${DOWNLOAD_PATH}/lineages" "${DOWNLOAD_PATH}/placement_files"
workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT

echo "==> Fetching BUSCO file manifest"
# wget, not curl: the BUSCO base image ships wget only, and wget follows the
# ezlab -> S3 redirect by default.
wget -q --tries=3 -O "${DOWNLOAD_PATH}/file_versions.tsv" "${S3_BASE}/file_versions.tsv"

# Column layout: name <tab> date <tab> md5 <tab> domain <tab> category
awk -F'\t' '$4=="Prokaryota" && $5=="lineages" && $1 ~ /odb10$/ {print $1"\t"$2"\t"$3}' \
    "${DOWNLOAD_PATH}/file_versions.tsv" > "${workdir}/lineages.tsv"
awk -F'\t' '$5=="placement_files" && $1 ~ /_odb10\./ {print $1"\t"$2"\t"$3}' \
    "${DOWNLOAD_PATH}/file_versions.tsv" > "${workdir}/placement.tsv"

n_lineages=$(wc -l < "${workdir}/lineages.tsv")
n_placement=$(wc -l < "${workdir}/placement.tsv")
echo "==> Resolved ${n_lineages} odb10 prokaryote lineages, ${n_placement} placement files"

if [ "${n_lineages}" -lt 90 ]; then
    echo "ERROR: expected ~99 odb10 prokaryote lineages, manifest yielded ${n_lineages}." >&2
    echo "The upstream manifest may have dropped odb10. Refusing to build a partial image." >&2
    exit 1
fi

# Fetch one file, verify md5, extract, discard the tarball. Retries on
# hash mismatch, which the upstream CDN occasionally serves under load.
fetch_one() {
    local url="$1" tarball="$2" expected_md5="$3" dest="$4"
    local got

    for attempt in 1 2 3; do
        if wget -q --tries=3 --waitretry=2 -O "${tarball}" "${url}"; then
            got="$(md5sum "${tarball}" | cut -d' ' -f1)"
            if [ "${got}" = "${expected_md5}" ]; then
                tar -xzf "${tarball}" -C "${dest}"
                rm -f "${tarball}"
                return 0
            fi
            echo "    md5 mismatch on $(basename "${url}") (attempt ${attempt})" >&2
        fi
        sleep $(( attempt * 3 ))
    done

    echo "ERROR: failed to fetch $(basename "${url}") after 3 attempts" >&2
    return 1
}
export -f fetch_one
export S3_BASE DOWNLOAD_PATH workdir

echo "==> Downloading lineages (${JOBS} parallel)"
# Lineage tarballs expand to a bare <name>/ directory, which is exactly what
# BUSCO's offline resolver expects under lineages/.
< "${workdir}/lineages.tsv" xargs -P "${JOBS}" -n 3 bash -c '
    fetch_one "${S3_BASE}/lineages/$0.$1.tar.gz" \
              "${workdir}/$0.tar.gz" "$2" "${DOWNLOAD_PATH}/lineages"
' || { echo "ERROR: one or more lineages failed to download" >&2; exit 1; }

echo "==> Downloading placement files"
# Placement archives are keyed as <basename>.<date>.<ext>.tar.gz and BUSCO's
# offline resolver globs for "<basename>.*<ext>" — so the extracted files must
# KEEP their embedded date. Do not rename them.
while IFS=$'\t' read -r name date md5; do
    base="${name%.*}"
    ext="${name##*.}"
    fetch_one "${S3_BASE}/placement_files/${base}.${date}.${ext}.tar.gz" \
              "${workdir}/${base}.${date}.${ext}.tar.gz" "${md5}" \
              "${DOWNLOAD_PATH}/placement_files"
done < "${workdir}/placement.tsv"

got_lineages=$(find "${DOWNLOAD_PATH}/lineages" -maxdepth 1 -mindepth 1 -type d | wc -l)
got_placement=$(find "${DOWNLOAD_PATH}/placement_files" -maxdepth 1 -type f | wc -l)

echo "==> Installed ${got_lineages} lineages, ${got_placement} placement files"
if [ "${got_lineages}" -lt "${n_lineages}" ] || [ "${got_placement}" -lt "${n_placement}" ]; then
    echo "ERROR: incomplete download (wanted ${n_lineages}/${n_placement})." >&2
    exit 1
fi

# Sanity-check the datasets BUSCO falls back on during auto-lineage.
for required in bacteria_odb10 archaea_odb10; do
    if [ ! -d "${DOWNLOAD_PATH}/lineages/${required}" ]; then
        echo "ERROR: required root dataset ${required} missing." >&2
        exit 1
    fi
done

du -sh "${DOWNLOAD_PATH}"
