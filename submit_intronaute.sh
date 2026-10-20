#!/bin/bash
#
# Submit intronaute.sh. CPUs, memory, time, and workers are set in that file.
#
# Usage (from login node):
#   ./submit_intronaute.sh --config /path/to/config.yaml
#   ./submit_intronaute.sh --config /path/to/config.yaml analyse
#   ./submit_intronaute.sh --config /path/to/config.yaml prepare
#   ./submit_intronaute.sh --config /path/to/config.yaml quantify
#   ./submit_intronaute.sh --local --config /path/to/config.yaml all
#   ./submit_intronaute.sh --force-prepare --config /path/to/config.yaml all
#   ./submit_intronaute.sh --hit-mode coverage --config /path/to/config.yaml all
#
# --config is required. The paths in that file are passed to intronaute.sh.
# Slurm logs are written to <outdir>/Logs.
# Slurm starts in the directory that contains these scripts.
#
# The block above is the help text. usage() prints lines 3-16, so leave
# that block where it is.
#
# -e  exit if a command fails
# -u  exit if a variable is used before it is set
# -o pipefail  a pipeline fails when any command in it fails, not only the last
set -euo pipefail

#########################################---Step 1: Define Functions---###########################################################
# Directory that holds intronaute.sh and load_deps.sh.
# The first candidate is the directory of this file. That is correct for
# ./submit_intronaute.sh on a login node.
# If those two files are not next to this script, use INTRONAUTE_SCRIPT_DIR
# when it is already set. Step 4 exports that variable for the job.
resolve_script_dir() {
  local cand="${BASH_SOURCE[0]}"
  cand="$(cd "$(dirname "${cand}")" && pwd)"
  if [[ -f "${cand}/intronaute.sh" && -f "${cand}/load_deps.sh" ]]; then
    printf '%s\n' "${cand}"
    return 0
  fi
  cand="${INTRONAUTE_SCRIPT_DIR:-}"
  if [[ -n "${cand}" && -f "${cand}/intronaute.sh" && -f "${cand}/load_deps.sh" ]]; then
    printf '%s\n' "$(cd "${cand}" && pwd)"
    return 0
  fi
  return 1
}

# First top-level "key:" line. Exit if the key is missing.
# Quotes around the value are removed. A value is required, so the
# capture is one or more characters.
yaml_value() {
  local key="$1" file="$2"
  local line
  line="$(grep -E "^${key}:" "$file" | head -1)" || true
  if [[ -z "${line}" ]]; then
    echo "Missing key '${key}' in ${file}" >&2
    exit 1
  fi
  sed -E 's/^[^:]+:[[:space:]]*"?([^"]+)"?[[:space:]]*$/\1/' <<<"${line}"
}

# Same lookup as yaml_value, but a missing key is empty and not an error.
# samplesheet and hit_mode are optional. The capture allows an empty value.
yaml_optional() {
  local key="$1" file="$2"
  local line
  line="$(grep -E "^${key}:" "$file" | head -1)" || true
  if [[ -z "${line}" ]]; then
    return 0
  fi
  sed -E 's/^[^:]+:[[:space:]]*"?([^"]*)"?[[:space:]]*$/\1/' <<<"${line}"
}

# Print one path per line for a YAML key. Accepts all of these:
#   bam_dir: /path
#   bam_dirs: [/path/a, /path/b]
#   bam_dirs:
#     - /path/a
#     - /path/b
# A block list ends at the first line that is not "  - ...".
# Quotes and surrounding whitespace are stripped. Blank entries are dropped.
yaml_paths() {
  local key="$1" file="$2"
  awk -v key="${key}" '
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      gsub(/^["'\'']|["'\'']$/, "", s)
      return s
    }
    function emit_flow(rest,   n, i, part, parts) {
      sub(/^\[/, "", rest)
      sub(/\]$/, "", rest)
      n = split(rest, parts, /,/)
      for (i = 1; i <= n; i++) {
        part = trim(parts[i])
        if (length(part) > 0) print part
      }
    }
    $0 ~ "^" key ":" {
      sub("^" key ":[[:space:]]*", "")
      if ($0 ~ /^\[/) {
        emit_flow($0)
      } else {
        $0 = trim($0)
        if (length($0) > 0) print $0
      }
      found = 1
      next
    }
    found && /^[[:space:]]*-[[:space:]]*/ {
      sub(/^[[:space:]]*-[[:space:]]*/, "")
      $0 = trim($0)
      if (length($0) > 0) print $0
      next
    }
    found { exit }
  ' "${file}"
}

# A path that starts with / is already absolute. Anything else, including
# an empty samplesheet, is resolved from the script directory. It is not
# resolved from the config file's directory or from the current directory.
abs_under_script() {
  local p="$1"
  if [[ -z "${p}" || "${p}" == /* ]]; then
    printf '%s\n' "${p}"
  else
    printf '%s\n' "${SCRIPT_DIR}/${p}"
  fi
}

# Print the header at the top of this file and exit.
# Lines 3-16 are the usage examples. Keep that header in place.
usage() {
  sed -n '3,16p' "$0" | sed 's/^# \?//'
  exit "${1:-0}"
}

#############################################---Step 2: Set Up Parameters---######################################################
SCRIPT_DIR="$(resolve_script_dir)" || {
  echo "Could not find intronaute.sh and load_deps.sh next to submit_intronaute.sh" >&2
  exit 1
}
CONFIG=""
LOCAL=0
# No stage name means the whole pipeline: prepare, discover, quantify, analyse.
STEP="all"
FORCE_QUANTIFY=0
FORCE_PREPARE=0
HIT_MODE=""
# 0 means "the user did not pass --hit-mode", so the config can still decide.
HIT_MODE_EXPLICIT=0

# Flags and the stage name can be mixed. The stage is one of
# prepare, quantify, analyse, or all.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage 0 ;;
    --local) LOCAL=1; shift ;;
    --config)
      if [[ $# -lt 2 ]]; then
        echo "--config needs a path" >&2
        usage 1
      fi
      CONFIG="$2"
      shift 2
      ;;
    --force-quantify) FORCE_QUANTIFY=1; shift ;;
    --force-prepare) FORCE_PREPARE=1; shift ;;
    --hit-mode)
      if [[ $# -lt 2 ]]; then
        echo "--hit-mode needs a value" >&2
        usage 1
      fi
      HIT_MODE="$2"
      HIT_MODE_EXPLICIT=1
      shift 2
      ;;
    prepare|quantify|analyse|all) STEP="$1"; shift ;;
    *) echo "Unknown argument: $1" >&2; usage 1 ;;
  esac
done

# --config is required. A relative config path is from the current
# directory, then canonicalized. YAML paths below use a different rule:
# they are resolved from SCRIPT_DIR by abs_under_script.
if [[ -z "${CONFIG}" ]]; then
  echo "--config is required" >&2
  usage 1
fi
if [[ "${CONFIG}" != /* ]]; then
  CONFIG="${PWD}/${CONFIG}"
fi
if [[ ! -f "${CONFIG}" ]]; then
  echo "Config not found: ${CONFIG}" >&2
  exit 1
fi
CONFIG="$(readlink -f "${CONFIG}")"

##########################################---Step 3: Read the Config---############################################################
# outdir, introns, and gtf are required. samplesheet is optional.
# targets.tsv.gz is not a config key. It always lives next to outdir:
#   <parent of outdir>/annotation/targets.tsv.gz
OUTDIR="$(abs_under_script "$(yaml_value outdir "${CONFIG}")")"
INTRONS="$(abs_under_script "$(yaml_value introns "${CONFIG}")")"
GTF="$(abs_under_script "$(yaml_value gtf "${CONFIG}")")"
SAMPLESHEET="$(abs_under_script "$(yaml_optional samplesheet "${CONFIG}")")"
TARGETS="$(dirname "${OUTDIR}")/annotation/targets.tsv.gz"

# Read bam_dir and bam_dirs. Both may be set. Drop duplicate absolute paths.
# BAM_DIR_SEEN is a newline-delimited set. Each path is stored with a
# newline on both sides so a shorter path cannot match inside a longer one.
mapfile -t BAM_DIR_LIST < <(yaml_paths bam_dir "${CONFIG}"; yaml_paths bam_dirs "${CONFIG}")
BAM_DIRS=""
BAM_DIR_SEEN=$'\n'
for d in "${BAM_DIR_LIST[@]+"${BAM_DIR_LIST[@]}"}"; do
  abs="$(abs_under_script "${d}")"
  case "${BAM_DIR_SEEN}" in
    *$'\n'"${abs}"$'\n'*) continue ;;
  esac
  BAM_DIR_SEEN+="${abs}"$'\n'
  BAM_DIRS+="${abs}"$'\n'
done

# A hit_mode line in the config is applied by the tool with the other
# analyse settings. When the line is absent, the job uses combined. The
# tool's own default is coverage, so that fallback is passed explicitly.
# --hit-mode on this command wins over the file.
# When the config does have hit_mode, leave HIT_MODE unset so intronaute.sh
# does not pass --hit-mode and the tool reads the config value itself.
if [[ "${HIT_MODE_EXPLICIT}" != "1" ]]; then
  if [[ -z "$(yaml_optional hit_mode "${CONFIG}")" ]]; then
    HIT_MODE="combined"
    HIT_MODE_EXPLICIT=1
  fi
fi
mkdir -p "${OUTDIR}/Logs"

#######################################---Step 4: Export the Job Environment---####################################################
# sbatch --export=ALL copies this environment into the job. intronaute.sh
# reads these names and refuses to start when the required ones are missing.
# BAM_DIRS is a newline-separated list. An empty samplesheet is still exported;
# intronaute.sh treats an empty value as "no samplesheet, run discover".
export INTRONAUTE_SCRIPT_DIR="${SCRIPT_DIR}"
export INTRONAUTE_CONFIG="${CONFIG}"
export INTRONAUTE_SAMPLESHEET="${SAMPLESHEET}"
export INTRONAUTE_BAM_DIRS="${BAM_DIRS}"
export INTRONAUTE_OUTDIR="${OUTDIR}"
export INTRONAUTE_INTRONS="${INTRONS}"
export INTRONAUTE_GTF="${GTF}"
export INTRONAUTE_TARGETS="${TARGETS}"
export INTRONAUTE_STEP="${STEP}"
export INTRONAUTE_FORCE_QUANTIFY="${FORCE_QUANTIFY}"
export INTRONAUTE_FORCE_PREPARE="${FORCE_PREPARE}"
if [[ "${HIT_MODE_EXPLICIT}" == "1" ]]; then
  export INTRONAUTE_HIT_MODE="${HIT_MODE}"
else
  unset INTRONAUTE_HIT_MODE
fi

###########################################---Step 5: Run or Submit---#############################################################
# --local runs in this shell. SLURM_JOB_ID is set when this script is
# already inside a job, and submitting again would nest a second job.
# exec replaces this process, so a local run does not continue past the pipeline.
# Otherwise sbatch submits intronaute.sh. --chdir is the script directory.
# Stdout and stderr share one log: <outdir>/Logs/intronaute_<jobid>.log.
if [[ "${LOCAL}" == "1" || -n "${SLURM_JOB_ID:-}" ]]; then
  exec "${SCRIPT_DIR}/intronaute.sh"
fi

sbatch \
  --export=ALL \
  --job-name="intronaute_${STEP}" \
  --chdir="${SCRIPT_DIR}" \
  --output="${OUTDIR}/Logs/intronaute_%j.log" \
  --error="${OUTDIR}/Logs/intronaute_%j.log" \
  "${SCRIPT_DIR}/intronaute.sh"
