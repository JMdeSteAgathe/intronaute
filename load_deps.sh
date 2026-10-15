#!/usr/bin/env bash
# Load this project's Python virtualenv and R.
#
# Source it so the current shell keeps the change:
#   source ./load_deps.sh            # ask before replacing commands already on PATH
#   source ./load_deps.sh --yes      # replace them
#   source ./load_deps.sh --no       # keep the commands already on PATH
#
# Python comes from .venv (intronaute, numpy, pandas, scipy, matplotlib).
# An Rscript already on PATH is kept when it starts. Otherwise R is taken
# from ~/.pixi/envs/r-base when that install exists. hits.rds needs only base R.
#
# --help prints lines 2-9 of this file. Leave that block where it is.

#########################################---Step 1: Confirm This File Was Sourced---###############################################
# Sourcing keeps VIRTUAL_ENV and PATH in the current shell. Running the file
# as a program would change a subshell and then throw the change away.
#
# Bash parses the whole file before it runs. The zsh-only expansion
# ${(%):-%x} would be a syntax error in bash, so it stays inside eval and
# runs only when this file is sourced from zsh.
_LOAD_DEPS_SELF=""
_load_deps_sourced=0
if [[ -n "${ZSH_VERSION:-}" ]]; then
  eval '_LOAD_DEPS_SELF=${(%):-%x}'
  # ZSH_EVAL_CONTEXT contains "file" when the file was sourced.
  case ":${ZSH_EVAL_CONTEXT:-}:" in
    *:file:*) _load_deps_sourced=1 ;;
  esac
else
  _LOAD_DEPS_SELF="${BASH_SOURCE[0]}"
  # In bash, $0 is the shell while BASH_SOURCE[0] is this file when sourced.
  if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    _load_deps_sourced=1
  fi
fi
if [[ "${_load_deps_sourced}" -eq 0 ]]; then
  echo "Source this file so the current shell keeps the environment:" >&2
  echo "  source ${_LOAD_DEPS_SELF} [--yes|--no]" >&2
  exit 1
fi

#########################################---Step 2: Define Functions---###########################################################
# One function so a failure can return to the sourcing shell instead of
# exiting it. The helpers are nested because they are removed before return.
_load_deps() {
  local root venv r_env choice arg
  local -a conflicts=()
  local i name current target reply front
  local python_conflict=0 r_conflict=0

  # This file lives next to .venv. R is the user's pixi env, not a copy in the repo.
  root="$(cd "$(dirname "${_LOAD_DEPS_SELF}")" && pwd)"
  venv="${root}/.venv"
  r_env="${HOME}/.pixi/envs/r-base"
  choice=""

#############################################---Step 3: Parse Arguments---#########################################################
  # --yes and --no both answer the PATH question. Only one of them is allowed.
  # --help reprints the header and returns, so the caller's shell stays open.
  for arg in "$@"; do
    case "${arg}" in
      -h|--help)
        sed -n '2,9p' "${_LOAD_DEPS_SELF}" | sed 's/^# \{0,1\}//'
        return 0
        ;;
      -y|--yes)
        [[ -z "${choice}" ]] || { echo "Pass only one of --yes or --no." >&2; return 1; }
        choice="yes"
        ;;
      -n|--no)
        [[ -z "${choice}" ]] || { echo "Pass only one of --yes or --no." >&2; return 1; }
        choice="no"
        ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        echo "Usage: source ${_LOAD_DEPS_SELF} [--yes|--no]" >&2
        return 1
        ;;
    esac
  done

  # The job cannot install packages. Build .venv before submitting:
  # wheels from requirements.txt, then pip install -e . --no-deps.
  if [[ ! -f "${venv}/bin/activate" ]]; then
    echo "Missing ${venv}." >&2
    echo "Create it from requirements.txt (wheel install, then pip install -e . --no-deps)." >&2
    return 1
  fi

  # type -P is bash-only. This file is sourced, so it may be running under zsh.
  # Both forms print the path of an external command and ignore functions and aliases.
  _load_deps_on_path() {
    if [[ -n "${ZSH_VERSION:-}" ]]; then
      whence -p "$1" 2>/dev/null || true
    else
      type -P "$1" 2>/dev/null || true
    fi
  }

  # Same command if the strings match or if both paths are the same file
  # after resolving symlinks. Either path being empty is not a match.
  _load_deps_same() {
    local a="$1" b="$2"
    [[ -n "${a}" && -n "${b}" ]] || return 1
    [[ "${a}" == "${b}" ]] && return 0
    [[ "$(readlink -f "${a}")" == "$(readlink -f "${b}")" ]]
  }

  # One conflict is name, current path, and replacement path, joined by
  # ASCII unit separator so a path that contains spaces still splits cleanly.
  _load_deps_add() {
    conflicts+=("$1"$'\x1f'"$2"$'\x1f'"$3")
  }

  # hits.rds is written with base R (read.delim and saveRDS). A binary that
  # exists but cannot start does not count as a usable Rscript.
  _load_deps_r_ready() {
    local bin="$1"
    [[ -n "${bin}" && -x "${bin}" ]] || return 1
    "${bin}" -e 'quit(status = 0)' >/dev/null 2>&1
  }

######################################---Step 4: Find Commands That Would Change---################################################
  # Skip the scan when this venv is already active. Otherwise every executable
  # in .venv/bin is a possible replacement. activate* scripts are not commands
  # the user would call, so they are not listed.
  if ! _load_deps_same "${VIRTUAL_ENV:-}" "${venv}"; then
    local cmd
    for cmd in "${venv}/bin"/*; do
      [[ -x "${cmd}" && ! -d "${cmd}" ]] || continue
      name="$(basename "${cmd}")"
      case "${name}" in
        activate|activate.csh|activate.fish|Activate.ps1) continue ;;
      esac
      current="$(_load_deps_on_path "${name}")"
      if [[ -n "${current}" ]] && ! _load_deps_same "${current}" "${cmd}"; then
        _load_deps_add "${name}" "${current}" "${cmd}"
        python_conflict=1
      fi
    done
  fi

  # An Rscript that already starts is kept, including when --yes was passed.
  # That is not a conflict. --yes only replaces commands recorded below.
  # No Rscript on PATH is also not a conflict: the pixi one is added later.
  # ~/.pixi/bin/Rscript and ~/.pixi/bin/R are pixi shims. They already point
  # at this env, so they are not treated as a foreign R.
  local want_r="" want_rscript="" have_r=0 r_ready=0
  current="$(_load_deps_on_path Rscript)"
  if _load_deps_r_ready "${current}"; then
    r_ready=1
  fi
  if [[ -x "${r_env}/bin/R" && -x "${r_env}/bin/Rscript" ]]; then
    have_r=1
    want_r="${r_env}/bin/R"
    want_rscript="${r_env}/bin/Rscript"
    if [[ "${r_ready}" -eq 1 ]]; then
      :
    elif [[ -z "${current}" ]]; then
      :
    elif [[ "${current}" == "${HOME}/.pixi/bin/Rscript" ]] || _load_deps_same "${current}" "${want_rscript}"; then
      :
    else
      _load_deps_add "Rscript" "${current}" "${want_rscript}"
      r_conflict=1
      current="$(_load_deps_on_path R)"
      if [[ -n "${current}" && "${current}" != "${HOME}/.pixi/bin/R" ]] && ! _load_deps_same "${current}" "${want_r}"; then
        _load_deps_add "R" "${current}" "${want_r}"
      fi
    fi
  fi

#####################################---Step 5: Ask Before Replacing Commands---###################################################
  # A login shell with no --yes or --no gets a prompt. The default is No.
  # A Slurm job has no terminal. intronaute.sh sources this file with --yes
  # so the job does not stop on the prompt. A source with no flag and no
  # terminal keeps the commands already on PATH.
  if [[ ${#conflicts[@]} -gt 0 ]]; then
    local row rest
    echo "Loading dependencies would replace these commands:" >&2
    for row in "${conflicts[@]}"; do
      name="${row%%$'\x1f'*}"
      rest="${row#*$'\x1f'}"
      current="${rest%%$'\x1f'*}"
      target="${rest#*$'\x1f'}"
      echo "  ${name}" >&2
      echo "    now:  ${current}" >&2
      echo "    load: ${target}" >&2
    done
    if [[ -z "${choice}" ]]; then
      if [[ -t 0 ]]; then
        printf 'Replace them? [y/N] ' >&2
        read -r reply || reply=""
        case "${reply}" in
          y|Y|yes|YES) choice="yes" ;;
          *) choice="no" ;;
        esac
      else
        echo "No terminal to ask. Existing commands were kept. Pass --yes or --no." >&2
        choice="no"
      fi
    elif [[ "${choice}" == "yes" ]]; then
      echo "Replacing them (--yes)." >&2
    else
      echo "Keeping them (--no)." >&2
    fi
  fi

##############################################---Step 6: Load Python---############################################################
  # Already active: say so and do not source activate again.
  # Conflict and the answer was not yes: leave PATH alone and print the
  # project interpreter so it can be called by its full path.
  # Otherwise source activate, which prepends .venv/bin. hash -r drops the
  # shell's remembered command locations so the new python is the one that runs.
  if _load_deps_same "${VIRTUAL_ENV:-}" "${venv}"; then
    echo "Python: ${venv} (already active)" >&2
  elif [[ "${python_conflict}" -eq 1 && "${choice}" != "yes" ]]; then
    echo "Python left unchanged. Project interpreter:" >&2
    echo "  ${venv}/bin/python" >&2
    echo "  ${venv}/bin/intronaute" >&2
  else
    # shellcheck source=/dev/null
    source "${venv}/bin/activate"
    hash -r 2>/dev/null || true
    echo "Python: ${venv}" >&2
  fi

################################################---Step 7: Load R---###############################################################
  # Neither a working Rscript nor a pixi r-base install: say so and continue.
  # The function still returns success. intronaute.sh is what requires R.
  # A working Rscript (r_ready) is kept even if pixi R also exists.
  # The pixi shim, or a path that is already the r-base Rscript, is already
  # the right R, so PATH is left as it is.
  # Otherwise, when pixi R exists and the user allowed the replacement (or
  # there was nothing to replace), symlink R and Rscript into .deps_bin and
  # put that directory first on PATH. The symlinks avoid activating pixi
  # itself. .deps_bin is gitignored.
  if [[ "${have_r}" -eq 0 && "${r_ready}" -eq 0 ]]; then
    echo "Rscript was not found." >&2
  else
    current="$(_load_deps_on_path Rscript)"
    local r_ok=0
    if [[ "${r_ready}" -eq 1 ]]; then
      r_ok=1
    elif [[ "${current}" == "${HOME}/.pixi/bin/Rscript" ]] || _load_deps_same "${current}" "${want_rscript}"; then
      r_ok=1
    fi
    if [[ "${r_ok}" -eq 0 && "${have_r}" -eq 1 && ( "${r_conflict}" -eq 0 || "${choice}" == "yes" ) ]]; then
      front="${root}/.deps_bin"
      command mkdir -p "${front}"
      ln -sfn "${want_r}" "${front}/R"
      ln -sfn "${want_rscript}" "${front}/Rscript"
      PATH="${front}:${PATH}"
      export PATH
      hash -r 2>/dev/null || true
      r_ok=1
      r_ready=0
    fi
    if [[ "${r_ok}" -eq 1 ]]; then
      echo "R: $(_load_deps_on_path Rscript)" >&2
    else
      echo "R left unchanged: ${current:-not on PATH}" >&2
    fi
  fi

  # Drop the helpers so a later source does not keep the previous copies.
  unset -f _load_deps_on_path _load_deps_same _load_deps_add _load_deps_r_ready
  return 0
}

###############################################---Step 8: Run the Loader---########################################################
# return, not exit. A failed source should set $? in the caller's shell
# and leave that shell running. The status is saved first because unset
# would replace it.
_load_deps "$@"
_load_deps_status=$?
unset -f _load_deps
unset _LOAD_DEPS_SELF _load_deps_sourced
return "${_load_deps_status}"
