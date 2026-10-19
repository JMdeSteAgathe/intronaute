#!/bin/bash
#
# Runs intronaute. submit_intronaute.sh submits this file and passes the
# config values in the environment. Edit the Slurm settings and WORKERS below.
#
#########################################---Step 1: Set Up Parameters---###########################################################
# Slurm reads the #SBATCH lines from this file. submit_intronaute.sh does
# not pass cpus, memory, or time; change them here.
# Account and partition are the Sherlock smontgom batch allocation.
#SBATCH --job-name=intronaute
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=smontgom

# Quantify processes. Passed as --workers, so this wins over a workers key
# in the config file. Keep it the same as --cpus-per-task unless you want fewer.
WORKERS=8

# -e  exit if a command fails
# -u  exit if a variable is used before it is set
# -o pipefail  a pipeline fails when any command in it fails, not only the last
set -euo pipefail

# submit_intronaute.sh exports these. The : command does nothing on success.
# ${VAR:?message} exits when VAR is unset or empty, and prints the message.
# A direct ./intronaute.sh therefore stops here instead of running with blank paths.
: "${INTRONAUTE_CONFIG:?Run this through submit_intronaute.sh}"
: "${INTRONAUTE_OUTDIR:?Run this through submit_intronaute.sh}"
: "${INTRONAUTE_INTRONS:?Run this through submit_intronaute.sh}"
: "${INTRONAUTE_GTF:?Run this through submit_intronaute.sh}"
: "${INTRONAUTE_TARGETS:?Run this through submit_intronaute.sh}"
: "${INTRONAUTE_STEP:?Run this through submit_intronaute.sh}"

#########################################---Step 2: Define Functions---###########################################################
# Slurm runs a copy of this file from its job spool. submit_intronaute.sh
# passes the real directory in INTRONAUTE_SCRIPT_DIR.
# Prefer the directory of this file when load_deps.sh is sitting next to it
# (a --local run). Otherwise use the directory the submit script recorded.
resolve_script_dir() {
  local cand="${BASH_SOURCE[0]}"
  cand="$(cd "$(dirname "${cand}")" && pwd)"
  if [[ -f "${cand}/load_deps.sh" ]]; then
    printf '%s\n' "${cand}"
    return 0
  fi
  cand="${INTRONAUTE_SCRIPT_DIR:-}"
  if [[ -n "${cand}" && -f "${cand}/load_deps.sh" ]]; then
    printf '%s\n' "$(cd "${cand}" && pwd)"
    return 0
  fi
  return 1
}

SCRIPT_DIR="$(resolve_script_dir)" || {
  echo "Could not find load_deps.sh next to intronaute.sh" >&2
  exit 1
}

# Copy the exported values into local names. samplesheet, the force flags,
# and hit mode may be empty. BAM directories arrive as one newline-separated
# string and are split back into an array. Blank lines are dropped.
CONFIG="${INTRONAUTE_CONFIG}"
SAMPLESHEET="${INTRONAUTE_SAMPLESHEET:-}"
OUTDIR="${INTRONAUTE_OUTDIR}"
INTRONS="${INTRONAUTE_INTRONS}"
GTF="${INTRONAUTE_GTF}"
TARGETS="${INTRONAUTE_TARGETS}"
STEP="${INTRONAUTE_STEP}"
FORCE_QUANTIFY="${INTRONAUTE_FORCE_QUANTIFY:-0}"
FORCE_PREPARE="${INTRONAUTE_FORCE_PREPARE:-0}"
HIT_MODE="${INTRONAUTE_HIT_MODE:-}"
BAM_DIRS=()
if [[ -n "${INTRONAUTE_BAM_DIRS:-}" ]]; then
  while IFS= read -r d; do
    [[ -n "${d}" ]] && BAM_DIRS+=("${d}")
  done <<< "${INTRONAUTE_BAM_DIRS}"
fi

# True only when this process already has the project venv.
# VIRTUAL_ENV must be this .venv, not some other environment that happens to be active.
project_deps_ready() {
  local venv
  venv="${SCRIPT_DIR}/.venv"

  [[ -f "${venv}/bin/activate" && -x "${venv}/bin/intronaute" ]] || return 1
  [[ -n "${VIRTUAL_ENV:-}" ]] || return 1
  [[ "$(readlink -f "${VIRTUAL_ENV}")" == "$(readlink -f "${venv}")" ]] || return 1
}

##########################################---Step 3: Load Dependencies---##########################################################
# Later commands write samplesheet.tsv into this directory, and intronaute
# is the console script from .venv. Start there so those paths are stable.
cd "${SCRIPT_DIR}"

# A batch job has no terminal, so load_deps.sh is sourced with --yes.
# That replaces a conflicting python instead of waiting for a prompt.
# The || exit 1 is here because a failing source under set -e does not
# always abort the caller.
if project_deps_ready; then
  echo "Dependencies already loaded; skipping load_deps.sh"
else
  echo "### load_deps.sh ###"
  # shellcheck source=/dev/null
  source "${SCRIPT_DIR}/load_deps.sh" --yes || exit 1
fi

###########################################---Step 4: Prepare Targets---###########################################################
# prepare builds targets.tsv.gz from the minor-intron list and the GTF.
# Run it when:
#   the stage is prepare, or
#   --force-prepare was passed, or
#   a later stage was asked for and the target table is not there yet.
# quantify or analyse alone skips prepare when the table already exists.
need_prepare=0
if [[ "${STEP}" == "prepare" ]]; then
  need_prepare=1
elif [[ "${FORCE_PREPARE}" == "1" ]]; then
  need_prepare=1
elif [[ "${STEP}" == "quantify" || "${STEP}" == "analyse" || "${STEP}" == "all" ]]; then
  if [[ ! -f "${TARGETS}" ]]; then
    need_prepare=1
  fi
fi

# introns and gtf are checked only when prepare will actually run.
# --force-prepare is forwarded so the tool rebuilds an existing table.
if [[ "${need_prepare}" == "1" ]]; then
  for f in "${INTRONS}" "${GTF}"; do
    if [[ ! -f "$f" ]]; then
      echo "Required for prepare, file not found: $f" >&2
      exit 1
    fi
  done
  mkdir -p "$(dirname "${TARGETS}")"
  echo "### intronaute prepare -> ${TARGETS} ###"
  prep_args=(--config "${CONFIG}" --introns "${INTRONS}" --gtf "${GTF}" --targets "${TARGETS}")
  if [[ "${FORCE_PREPARE}" == "1" ]]; then
    prep_args+=(--force-prepare)
  fi
  intronaute prepare "${prep_args[@]}"
fi

# prepare on its own stops here. The samplesheet and the counts are not touched.
if [[ "${STEP}" == "prepare" ]]; then
  echo "Done. Targets: ${TARGETS}"
  exit 0
fi

# quantify and analyse both read the target table. Fail before discover
# if prepare was supposed to write it and the file is still missing.
if [[ ! -f "${TARGETS}" ]]; then
  echo "Targets not found after prepare: ${TARGETS}" >&2
  exit 1
fi

###########################################---Step 5: Discover Samples---##########################################################
# A samplesheet in the config wins. discover is skipped even if bam_dir
# or bam_dirs is also set. The file has to exist; this script does not
# create a samplesheet that the config already named.
if [[ -n "${SAMPLESHEET}" ]]; then
  if [[ ${#BAM_DIRS[@]} -gt 0 ]]; then
    echo "samplesheet is set; skipping discover (${SAMPLESHEET})"
  fi
  if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "Samplesheet not found: ${SAMPLESHEET}" >&2
    exit 1
  fi
else
  # No samplesheet: discover walks every BAM directory and writes
  # samplesheet.tsv next to these scripts. That path is gitignored.
  if [[ ${#BAM_DIRS[@]} -eq 0 ]]; then
    echo "Set samplesheet, bam_dir, or bam_dirs in ${CONFIG}" >&2
    exit 1
  fi
  for d in "${BAM_DIRS[@]}"; do
    if [[ ! -d "${d}" ]]; then
      echo "BAM directory not found: ${d}" >&2
      exit 1
    fi
  done

  SAMPLESHEET="${SCRIPT_DIR}/samplesheet.tsv"
  echo "### intronaute discover -> ${SAMPLESHEET} ###"
  intronaute discover --config "${CONFIG}" --bam-dir "${BAM_DIRS[@]}" --samplesheet "${SAMPLESHEET}"

  if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "discover did not write ${SAMPLESHEET}" >&2
    exit 1
  fi
fi

# Shared by quantify and analyse. WORKERS and HIT_MODE are added below
# because they apply to only one of the two stages.
common=(--config "${CONFIG}" --targets "${TARGETS}" --samplesheet "${SAMPLESHEET}" --outdir "${OUTDIR}")

##############################################---Step 6: Quantify---###############################################################
# Count each target in each sample. --workers is WORKERS from the top of
# this file, so a workers key in the config does not change the job size.
# --force re-counts. Without it, finished samples are reused from the cache.
# The stage name "all" includes this step. "analyse" and "prepare" do not.
if [[ "${STEP}" == "quantify" || "${STEP}" == "all" ]]; then
  echo "### intronaute quantify (workers=${WORKERS}) ###"
  qargs=("${common[@]}" --workers "${WORKERS}")
  if [[ "${FORCE_QUANTIFY}" == "1" ]]; then
    qargs+=(--force)
  fi
  intronaute quantify "${qargs[@]}"
fi

##############################################---Step 7: Analyse---################################################################
# Score the cohort and write the HTML report. --hit-mode is passed only
# when submit_intronaute.sh exported one: an explicit --hit-mode, or the
# combined fallback used when the config has no hit_mode line.
# An empty HIT_MODE means the config has hit_mode, and the tool reads it.
# analyse writes ${OUTDIR}/report.html. A quantify-only run still prints
# that path at the end; the file is produced by this step.
if [[ "${STEP}" == "analyse" || "${STEP}" == "all" ]]; then
  aargs=("${common[@]}")
  if [[ -n "${HIT_MODE}" ]]; then
    aargs+=(--hit-mode "${HIT_MODE}")
    echo "### intronaute analyse (hit-mode=${HIT_MODE}) ###"
  else
    echo "### intronaute analyse (hit-mode from config) ###"
  fi
  intronaute analyse "${aargs[@]}"
fi

echo "Done. Report: ${OUTDIR}/report.html"
