#!/usr/bin/env bash

if (( $# < 1 || $# > 2 )); then
  echo "Usage: $0 <work_dir> [output_csv]" >&2
  exit 1
fi

work_dir=$1
out_csv=${2:-mavedb_logs.csv}

[[ -d "$work_dir" ]] || { echo "No such directory: $work_dir" >&2; exit 2; }
out_dir=$(dirname "$out_csv")
mkdir -p "$out_dir"

# Gather all .command.log files
mapfile -d '' LOGS < <(find "$work_dir" -type f -name '.command.log' -print0 2>/dev/null || true)

if ((${#LOGS[@]}==0)); then
  echo "No .command.log files found under: $work_dir" >&2
  echo '"file","time","urn","step","reason","subid"' > "$out_csv"
  exit 0
fi

tmp_keys="$(mktemp)"
trap 'rm -f "$tmp_keys"' EXIT

# discover all dynamic key names
gawk '
  function get_keys(rest,    m2,key) {
    while (match(rest, /[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)=("([^"]*)"|[^[:space:]]+)/, m2)) {
      key = m2[1]
      if (!(key in seen)) { print key; seen[key]=1 }
      rest = substr(rest, RSTART + RLENGTH)
    }
  }
  {
    if (match($0, /^\[([^]]+)\]\[MaveDB\]\[URN=([^]]+)\]\[STEP=([^]]+)\]\[REASON=([^]]+)\]\[SUBID=([^]]+)\][[:space:]]*(.*)$/, m)) {
      get_keys(m[6])
    }
  }
' "${LOGS[@]}" | LC_ALL=C sort -u > "$tmp_keys"

# generate CSV
gawk -v KEYFILE="$tmp_keys" -v OFS=',' '
  function esc(s){ gsub(/"/, "\"\"", s); return "\"" s "\"" }

  BEGIN {
    nKeys=0
    while ((getline k < KEYFILE) > 0) if (k!="") keys[++nKeys]=k
    close(KEYFILE)

    printf "%s", "\"file\",\"time\",\"urn\",\"step\",\"reason\",\"subid\""
    for (i=1; i<=nKeys; i++) printf ",\"%s\"", keys[i]
    print ""
  }

  function parse_extras(rest,    m2,k,v) {
    delete extra
    while (match(rest, /[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)=("([^"]*)"|[^[:space:]]+)/, m2)) {
      k = m2[1]
      if (substr(m2[2],1,1)=="\"") v=m2[3]; else v=m2[2]
      extra[k]=v
      rest = substr(rest, RSTART + RLENGTH)
    }
  }

  {
    if (!match($0, /^\[([^]]+)\]\[MaveDB\]\[URN=([^]]+)\]\[STEP=([^]]+)\]\[REASON=([^]]+)\]\[SUBID=([^]]+)\][[:space:]]*(.*)$/, m))
      next

    fpath = FILENAME
    t     = m[1]
    urn   = m[2]
    step  = m[3]
    reason= m[4]
    subid = m[5]
    rest  = m[6]

    parse_extras(rest)

    printf "%s,%s,%s,%s,%s,%s", esc(fpath), esc(t), esc(urn), esc(step), esc(reason), esc(subid)
    for (i=1; i<=nKeys; i++) {
      val = (keys[i] in extra) ? extra[keys[i]] : ""
      printf ",%s", esc(val)
    }
    print ""
  }
' "${LOGS[@]}" > "$out_csv"

echo "Wrote: $out_csv"
