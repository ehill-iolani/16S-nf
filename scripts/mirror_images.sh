#!/usr/bin/env bash
#
# Copies every container image this pipeline uses into our own Artifact Registry
# repo, so a Batch run doesn't depend on quay.io / Docker Hub answering from every
# VM at the moment its task starts (a single timed-out pull there once ended a
# 2-hour, 8,000-task run).
#
# The image list is read from the pipeline itself -- the `container` directives in
# modules/ and nextflow.config -- so it can't drift from what actually runs. Each
# image keeps its upstream path under the mirror
# (quay.io/biocontainers/blast:X -> $MIRROR/biocontainers/blast:X), which is what
# lets the google_batch profile swap only the registry host (params
# container_registry / dockerhub_registry).
#
#   scripts/mirror_images.sh --list       what would be mirrored, and where
#   scripts/mirror_images.sh --check      read-only: is each image in the mirror, and is
#                                         upstream's exact image what it holds? exits 1 if not
#   scripts/mirror_images.sh --dry-run    print the copy commands without running them
#   scripts/mirror_images.sh              copy what's missing (tags are pinned, so an
#                                         image already present is left alone)
#   scripts/mirror_images.sh --force      re-copy everything
#
# Needs docker with buildx, and `gcloud auth configure-docker us-central1-docker.pkg.dev`
# once for pushing. Override the destination with MIRROR_REGISTRY or --dest=HOST/PROJECT/REPO.
set -euo pipefail
cd "$(dirname "$0")/.."

MIRROR="${MIRROR_REGISTRY:-us-central1-docker.pkg.dev/iscc-400300/containers}"
mode=copy
force=0

usage() { sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; }

for arg in "$@"; do
    case "$arg" in
        --list)    mode=list ;;
        --check)   mode=check ;;
        --dry-run) mode=dry ;;
        --force)   force=1 ;;
        --dest=*)  MIRROR="${arg#--dest=}" ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $arg" >&2; usage >&2; exit 2 ;;
    esac
done

files=(modules/*.nf nextflow.config)
[ -d conf ] && files+=(conf/*.config)

# A `container` that doesn't go through the registry params would silently bypass
# the mirror, so refuse to carry on rather than copy an incomplete list.
hardcoded=$(grep -nE "^[^/]*container[[:space:]]*=?[[:space:]]*(\{[[:space:]]*)?['\"]" "${files[@]}" \
            | grep -v '\${params\.\(container\|dockerhub\)_registry}' || true)
if [ -n "$hardcoded" ]; then
    echo "container directives that bypass --container_registry / --dockerhub_registry:" >&2
    echo "$hardcoded" >&2
    exit 1
fi

# "${params.container_registry}/biocontainers/blast:X" -> quay.io/biocontainers/blast:X
sources=$(grep -hoE '"\$\{params\.(container|dockerhub)_registry\}/[^"]+"' "${files[@]}" \
          | tr -d '"' \
          | sed -e 's#\${params\.container_registry}#quay.io#' -e 's#\${params\.dockerhub_registry}#docker.io#' \
          | sort -u)
[ -n "$sources" ] || { echo "no container images found" >&2; exit 1; }

# These read buildx's plain-text output: its --format templates aren't reliable.
# digest: the top-level manifest digest (of the multi-arch index, when there is one).
digest() { docker buildx imagetools inspect "$1" 2>/dev/null | awk '/^Digest:/ { print $2; exit }'; }
# digests: that one plus the digest of every manifest inside it. Copying a
# single-platform image wraps it in a new manifest list, so the mirror's top-level
# digest differs from upstream's even though the image underneath is byte-identical:
# it counts as the same image when upstream's digest is among these.
digests() { docker buildx imagetools inspect "$1" 2>/dev/null | grep -oE 'sha256:[0-9a-f]{64}' | sort -u; }

missing=0
while read -r src; do
    dest="$MIRROR/${src#*/}"   # keep the upstream path, drop only its host
    case "$mode" in
        list)
            echo "$src -> $dest"
            ;;
        dry)
            echo "docker buildx imagetools create --tag $dest $src"
            ;;
        check)
            have=$(digests "$dest" || true)
            if [ -z "$have" ]; then
                echo "MISSING  $dest"; missing=$((missing + 1))
            elif ! grep -qx "$(digest "$src")" <<< "$have"; then
                echo "DIFFERS  $dest (upstream is $(digest "$src"), not found in it)"; missing=$((missing + 1))
            else
                echo "ok       $dest"
            fi
            ;;
        copy)
            if [ "$force" -eq 0 ] && [ -n "$(digest "$dest" || true)" ]; then
                echo "present  $dest"
            else
                echo "copying  $src -> $dest"
                docker buildx imagetools create --tag "$dest" "$src"
            fi
            ;;
    esac
done <<< "$sources"

if [ "$mode" = check ] && [ "$missing" -gt 0 ]; then
    echo "$missing image(s) missing or different -- run scripts/mirror_images.sh" >&2
    exit 1
fi
