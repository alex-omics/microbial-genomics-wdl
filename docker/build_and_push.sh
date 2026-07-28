#!/usr/bin/env bash
#
# build_and_push.sh — build microbial-genomics-wdl Docker images and push them
# to Docker Hub.
#
# Usage:
#   ./docker/build_and_push.sh                      # build all, don't push
#   ./docker/build_and_push.sh --push               # build all, push each
#   ./docker/build_and_push.sh microgwas-mash       # build just one
#   ./docker/build_and_push.sh --push microgwas-mash microgwas-mlst
#   ./docker/build_and_push.sh --list               # show what would be built
#
# Notes for running this:
#   - Several images install from conda.anaconda.org, which is blocked on the
#     MGH network. Those are marked NEEDS_CONDA below; run them from a network
#     that can reach anaconda (hotspot / home wifi).
#   - --platform linux/amd64 is not optional. Terra/GCP runs amd64, and on
#     Apple Silicon a native build produces arm64 images that will fail there
#     with "exec format error". Building amd64 on a Mac uses emulation, so it
#     is slow but correct.
#   - `docker login` once before using --push.

set -euo pipefail

NAMESPACE="aarvani1"
PLATFORM="linux/amd64"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# image_directory : tag : needs_conda
IMAGES=(
  "microgwas-base:0.9.1:yes"
  "microgwas-pyseer:1.4.2:yes"
  "microgwas-enrich:0.9.1:yes"
  "microgwas-limix:3.0.4:yes"
  "microgwas-mash:2.1:no"
  "microgwas-mlst:2.16.2:no"
  "microgwas-bcftools:1.13:no"
  "microgwas-snpsites:2.5.1:no"
  "microgwas-snippy:4.6.0:no"
  "microgwas-panfeed:1.6.1:no"
  "microgwas-eggnog:2.1.13:no"
)

PUSH=false
SELECTED=()

for arg in "$@"; do
  case "$arg" in
    --push) PUSH=true ;;
    --list)
      printf '%-24s %-10s %s\n' "IMAGE" "TAG" "NEEDS CONDA"
      for entry in "${IMAGES[@]}"; do
        IFS=':' read -r name tag conda <<< "$entry"
        printf '%-24s %-10s %s\n' "$name" "$tag" "$conda"
      done
      exit 0
      ;;
    -*) echo "Unknown flag: $arg" >&2; exit 1 ;;
    *) SELECTED+=("$arg") ;;
  esac
done

if $PUSH; then
  # Advisory pre-flight only: warn before spending 20 minutes building, but do
  # not block on it. There is no reliable, stable way to ask Docker "am I
  # logged in" — `docker system info` stopped emitting a Username: line in
  # recent Docker Desktop, and credentials may live in a credential helper
  # rather than in config.json. The real authority is the push itself, so a
  # guess about internal state must never be able to veto a valid run.
  logged_in=false

  creds_store="$(python3 -c "
import json, os
p = os.path.expanduser('~/.docker/config.json')
try:
    print(json.load(open(p)).get('credsStore', ''))
except Exception:
    print('')
" 2>/dev/null || true)"

  if [[ -n "$creds_store" ]] && command -v "docker-credential-${creds_store}" > /dev/null 2>&1; then
    if "docker-credential-${creds_store}" list 2>/dev/null | grep -q 'index.docker.io'; then
      logged_in=true
    fi
  fi

  if ! $logged_in && grep -q 'index.docker.io' "${HOME}/.docker/config.json" 2>/dev/null; then
    logged_in=true
  fi

  if ! $logged_in; then
    echo "WARNING: could not confirm a Docker Hub login from this machine's config." >&2
    echo "         If the push fails with 'denied', run: docker login" >&2
    echo "         Continuing anyway — this check is advisory, not authoritative." >&2
    echo >&2
  fi
fi

built=()
pushed=()
skipped=()
failed=()

for entry in "${IMAGES[@]}"; do
  IFS=':' read -r name tag needs_conda <<< "$entry"

  # Honour an explicit selection, if one was given
  if [[ ${#SELECTED[@]} -gt 0 ]]; then
    match=false
    for s in "${SELECTED[@]}"; do [[ "$s" == "$name" ]] && match=true; done
    $match || continue
  fi

  image="${NAMESPACE}/${name}:${tag}"
  echo
  echo "=============================================================="
  echo "Building ${image}"
  [[ "$needs_conda" == "yes" ]] && echo "  (needs conda.anaconda.org — blocked on the MGH network)"
  echo "=============================================================="

  if docker buildx build \
       --platform "$PLATFORM" \
       -t "$image" \
       --load \
       "${SCRIPT_DIR}/${name}/"; then
    built+=("$image")
  else
    failed+=("$image")
    # Keep going: one unreachable package index should not stop the rest
    continue
  fi

  if $PUSH; then
    if docker push "$image"; then
      pushed+=("$image")
    else
      failed+=("${image} (push)")
    fi
  else
    skipped+=("$image")
  fi
done

echo
echo "=============================================================="
echo "SUMMARY"
echo "=============================================================="
[[ ${#built[@]}   -gt 0 ]] && { echo "Built:";        printf '  %s\n' "${built[@]}"; }
[[ ${#pushed[@]}  -gt 0 ]] && { echo "Pushed:";       printf '  %s\n' "${pushed[@]}"; }
[[ ${#skipped[@]} -gt 0 ]] && { echo "Built, not pushed (re-run with --push):"; printf '  %s\n' "${skipped[@]}"; }
[[ ${#failed[@]}  -gt 0 ]] && { echo "FAILED:";       printf '  %s\n' "${failed[@]}"; }

[[ ${#failed[@]} -gt 0 ]] && exit 1
exit 0
