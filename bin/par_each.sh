# Sourced helper (`. par_each.sh`), not an executable: run a bash function over
# a list of files, N at a time.
#
#   par_each <jobs> <function> <file>...
#
# The per-sample consensus processes loop over every cluster in one task rather
# than spawning a task (a VM, an image pull, a GCS stage-in) per cluster, and
# use this to keep all of the task's cpus busy while they do. Files are dealt
# round-robin to <jobs> workers; if any call fails the whole thing returns
# non-zero. Plain `xargs -P` isn't used because busybox builds (which some
# biocontainers images are) don't all support it.
#
# Call the function as a bare command: the shell ignores `set -e` inside a
# function that is tested with `||`/`if`, so a failed step in it would be
# silently skipped.
par_each() {
    local n=$1 fn=$2
    shift 2
    local w pids=() rc=0
    for ((w = 0; w < n; w++)); do
        (
            i=0
            for f in "$@"; do
                if (( i % n == w )); then
                    "$fn" "$f"
                fi
                i=$((i + 1))
            done
        ) &
        pids+=($!)
    done
    for p in "${pids[@]}"; do
        wait "$p" || rc=1
    done
    return $rc
}
