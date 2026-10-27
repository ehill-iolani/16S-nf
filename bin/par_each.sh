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
#
# The function runs inside this one's scope, so a plain assignment in it (say
# `n=...`) overwrites a variable of the same name here. Declare the function's
# own variables `local`; the _pe_ names below are a second guard. In 16S-nf
# v0.3.4 a clobbered worker count left clusters drafted twice or not at all.
par_each() {
    local _pe_n=$1 _pe_fn=$2
    shift 2
    local _pe_w _pe_p _pe_pids=() _pe_rc=0
    for ((_pe_w = 0; _pe_w < _pe_n; _pe_w++)); do
        (
            _pe_i=0
            for _pe_f in "$@"; do
                if (( _pe_i % _pe_n == _pe_w )); then
                    "$_pe_fn" "$_pe_f"
                fi
                _pe_i=$((_pe_i + 1))
            done
        ) &
        _pe_pids+=($!)
    done
    for _pe_p in "${_pe_pids[@]}"; do
        wait "$_pe_p" || _pe_rc=1
    done
    return $_pe_rc
}
