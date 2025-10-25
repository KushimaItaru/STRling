#!/bin/bash
# 35_strling_outliers_v7_1.sh
#
# STRling outliers：染色体ごとに実行し、サンプル別 *.STRs.tsv を
# ノードローカルの作業ディレクトリで生成 → その場で集約して
# str-results/outliers/STRs_chrN.tsv（1染色体=1ファイル）を出力する。
#
# v7.1 のポイント
# - race condition解消：parallel で共有ファイルに同時追記しない
# - ripgrep 出力から行番号を除去（--no-line-number）
# - 変数名修正：DST→DEST_PS
# - BLAS/NumPyのスレッド固定（オーバーサブスク回避）
# - tmp出力→成功時mv（ゼロ長TSVを残さない）
#
#SBATCH -J strling_outliers_v7_1
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-24%12
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 8:00:00
#SBATCH --mem=512G
#SBATCH -o /home/kushima-pg/strling_10142025/logfiles/outliers_v7_1_chr_%A_%a.out
#SBATCH -e /home/kushima-pg/strling_10142025/logfiles/outliers_v7_1_chr_%A_%a.err

set -euo pipefail

# ---- CPUスレッドの暴走防止（Python側の内部並列を制限）----
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export LC_ALL=C LANG=C

SCRIPT_START=$(date +%s)

echo "========================================"
echo "STRling Outlier Analysis v7.1"
echo "========================================"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "JobID: ${SLURM_JOB_ID:-N/A}"
echo "ArrayTaskID: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "CPUs: ${SLURM_CPUS_ON_NODE:-N/A}"
echo "Mem:  ${SLURM_MEM_PER_NODE:-N/A}"
echo ""

# -------- 共通設定 --------
CFG="/home/kushima-pg/strling_10142025/00_config_strling_v2.sh"
[ -f "$CFG" ] || { echo "[ERROR] Config not found: $CFG" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CFG"

mkdir -p "${STR_RES_DIR}" "${LOG_DIR}"
OUTLIER_DIR="${STR_RES_DIR}/outliers"
PER_SAMPLE_DIR="${OUTLIER_DIR}/per-sample"
mkdir -p "${OUTLIER_DIR}" "${PER_SAMPLE_DIR}"

# -------- 染色体ID --------
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
if   [ "$TASK_ID" -le 22 ]; then CHR="chr${TASK_ID}"
elif [ "$TASK_ID" -eq 23 ]; then CHR="chrX"
elif [ "$TASK_ID" -eq 24 ]; then CHR="chrY"
else echo "[ERROR] Invalid array id: ${TASK_ID}" >&2; exit 1
fi
echo "[INFO] CHR=${CHR}"

OUTPUT_TSV="${OUTLIER_DIR}/STRs_${CHR}.tsv"
OUTPUT_LOG="${LOG_DIR}/outliers_${CHR}_v7_1.log"

# 既存スキップ
if [ -s "${OUTPUT_TSV}" ]; then
  FSZ=$(stat -c%s "${OUTPUT_TSV}" 2>/dev/null || echo 0)
  if [ "${FSZ}" -ge $((500*1024)) ]; then
    echo "[INFO] Already present: ${OUTPUT_TSV} ($(printf '%.1f' "$(echo "${FSZ}/1048576" | bc -l)") MB)"
    exit 0
  else
    echo "[WARN] Existing tiny file. Rebuilding: ${OUTPUT_TSV}"
    rm -f "${OUTPUT_TSV}"
  fi
fi

# -------- 作業ディレクトリ（ノードローカル） --------
LOCAL_TMP="${SLURM_TMPDIR:-/tmp}/strling_${CHR}_$$"
WORKDIR="${LOCAL_TMP}/work_${CHR}"
mkdir -p "${WORKDIR}"
export LOCAL_TMP CHR WORKDIR

echo "[INFO] WORKDIR=${WORKDIR}"

# -------- 入力列挙 --------
mapfile -d '' GENO_FILES < <(find "${STR_RES_DIR}" -maxdepth 1 -type f -name "*-genotype.txt" -print0)
[ ${#GENO_FILES[@]} -gt 0 ] || { echo "[ERROR] No genotype files in ${STR_RES_DIR}" >&2; exit 1; }
echo "[INFO] Genotypes: ${#GENO_FILES[@]} files"

# -------- 抽出ツール & 並列度 --------
EXTRACT_J=${SLURM_CPUS_ON_NODE:-16}
EXTRACT_J=$(( EXTRACT_J>2 ? EXTRACT_J-2 : 1 ))
USE_RG="false"; USE_MAWK="false"
if command -v rg >/dev/null 2>&1; then USE_RG="true"
elif command -v mawk >/dev/null 2>&1; then USE_MAWK="true"; fi
echo "[INFO] Parallel extract: ${EXTRACT_J}  (rg=${USE_RG}, mawk=${USE_MAWK})"

# -------- 染色体抽出（並列・race free） --------
echo "[INFO] Extracting ${CHR}..."

# 個別リストの集約先
GENO_LIST="${LOCAL_TMP}/geno_list.txt"
rm -f "${GENO_LIST}"

if command -v parallel >/dev/null 2>&1; then
  if [ "${USE_RG}" = "true" ]; then
    printf '%s\0' "${GENO_FILES[@]}" | parallel -0 -j "${EXTRACT_J}" '
      f="{}"; id=$(basename "$f" -genotype.txt)
      out="'"${LOCAL_TMP}"'/${id}-genotype.txt"
      list_tmp="'"${LOCAL_TMP}"'/geno_${id}.list"
      { head -n 1 "$f"; rg --no-heading --no-line-number -e "^'"${CHR}"'\t" "$f"; } > "$out"
      if [ "$(wc -l < "$out")" -gt 1 ]; then
        printf "%s\n" "$out" > "$list_tmp"
      else
        rm -f "$out" "$list_tmp"
      fi
    '
  elif [ "${USE_MAWK}" = "true" ]; then
    printf '%s\0' "${GENO_FILES[@]}" | parallel -0 -j "${EXTRACT_J}" '
      f="{}"; id=$(basename "$f" -genotype.txt)
      out="'"${LOCAL_TMP}"'/${id}-genotype.txt"
      list_tmp="'"${LOCAL_TMP}"'/geno_${id}.list"
      mawk -v chr="'"${CHR}"'" "NR==1 || \$1==chr" "$f" > "$out"
      if [ "$(wc -l < "$out")" -gt 1 ]; then
        printf "%s\n" "$out" > "$list_tmp"
      else
        rm -f "$out" "$list_tmp"
      fi
    '
  else
    printf '%s\0' "${GENO_FILES[@]}" | parallel -0 -j "${EXTRACT_J}" '
      f="{}"; id=$(basename "$f" -genotype.txt)
      out="'"${LOCAL_TMP}"'/${id}-genotype.txt"
      list_tmp="'"${LOCAL_TMP}"'/geno_${id}.list"
      awk -v chr="'"${CHR}"'" "NR==1 || \$1==chr" "$f" > "$out"
      if [ "$(wc -l < "$out")" -gt 1 ]; then
        printf "%s\n" "$out" > "$list_tmp"
      else
        rm -f "$out" "$list_tmp"
      fi
    '
  fi
  # 結合（重複除去）
  cat "${LOCAL_TMP}"/geno_*.list 2>/dev/null | sort -u > "${GENO_LIST}" || : 
  rm -f "${LOCAL_TMP}"/geno_*.list
else
  : > "${GENO_LIST}"
  for f in "${GENO_FILES[@]}"; do
    id=$(basename "$f" -genotype.txt)
    out="${LOCAL_TMP}/${id}-genotype.txt"
    if [ "${USE_RG}" = "true" ]; then
      { head -n 1 "$f"; rg --no-heading --no-line-number -e "^${CHR}\t" "$f"; } > "$out"
    elif [ "${USE_MAWK}" = "true" ]; then
      mawk -v chr="${CHR}" 'NR==1 || $1==chr' "$f" > "$out"
    else
      awk -v chr="${CHR}" 'NR==1 || $1==chr' "$f" > "$out"
    fi
    if [ "$(wc -l < "$out")" -gt 1 ]; then printf "%s\n" "$out" >> "${GENO_LIST}"; else rm -f "$out"; fi
  done
fi

FILTERED_N=$(wc -l < "${GENO_LIST}" 2>/dev/null || echo 0)
[ "${FILTERED_N}" -gt 0 ] || { echo "[ERROR] No ${CHR} data found in any genotype" >&2; rm -rf "${LOCAL_TMP}"; exit 1; }
echo "[INFO] Samples with ${CHR} data: ${FILTERED_N}"

# -------- unplaced（該当サンプルのみ） --------
: > "${LOCAL_TMP}/unplaced_list.txt"
while IFS= read -r gfile; do
  sid=$(basename "${gfile}" -genotype.txt)
  u="${STR_RES_DIR}/${sid}-unplaced.txt"
  [ -s "$u" ] && printf "%s\n" "$u" >> "${LOCAL_TMP}/unplaced_list.txt"
done < "${GENO_LIST}"
UNPL_COUNT=$(wc -l < "${LOCAL_TMP}/unplaced_list.txt" 2>/dev/null || echo 0)
echo "[INFO] Unplaced files: ${UNPL_COUNT}"

# -------- outliers 実行体の解決 --------
PYBIN="python"; command -v pypy3 >/dev/null 2>&1 && PYBIN="pypy3"
resolve_outliers_cmd() {
  local target="${STRLING_OUTLIERS:-strling-outliers.py}"
  if command -v "$target" >/dev/null 2>&1; then echo "exe::$target"; return 0
  elif [ -f "$target" ]; then case "$target" in *.py) echo "py::$target";; *) echo "exe::$target";; esac; return 0; fi
  echo "mod::strling_outliers"
}
MODE_AND_PATH=$(resolve_outliers_cmd); MODE=${MODE_AND_PATH%%::*}; OPATH=${MODE_AND_PATH##*::}
echo "[INFO] Outliers runner: ${MODE}(${OPATH})" | tee -a "${OUTPUT_LOG}"

# -------- outliers 実行（WORKDIRでサンプル別ファイルを生成） --------
mapfile -t GENO_ARGS < "${GENO_LIST}"
UNPL_ARGS=(); [ -s "${LOCAL_TMP}/unplaced_list.txt" ] && mapfile -t UNPL_ARGS < "${LOCAL_TMP}/unplaced_list.txt"

build_cmd_repeat() {
  local -n _GA=$1; local -n _UA=$2; local _cmd=()
  case "$MODE" in
    exe) _cmd=( "$OPATH" );;
    py)  _cmd=( "$PYBIN" "$OPATH" );;
    mod) _cmd=( "$PYBIN" -m "$OPATH" );;
  esac
  local g; for g in "${_GA[@]}"; do _cmd+=( "--genotypes" "$g" ); done
  if [ ${#_UA[@]} -gt 0 ]; then _cmd+=( "--unplaced" ); local u; for u in "${_UA[@]}"; do _cmd+=( "$u" ); done; fi
  printf '%s\0' "${_cmd[@]}"
}

pushd "${WORKDIR}" >/dev/null
readarray -d '' CMD_ARR < <(build_cmd_repeat GENO_ARGS UNPL_ARGS)
echo "[INFO] Running outliers for ${CHR} ..." | tee -a "${OUTPUT_LOG}"
set +e
"${CMD_ARR[@]}" >> "${OUTPUT_LOG}" 2>&1
RET=$?
set -e
popd >/dev/null
[ $RET -eq 0 ] || { echo "[ERROR] outliers exit=${RET}" | tee -a "${OUTPUT_LOG}"; rm -rf "${LOCAL_TMP}"; exit 1; }

# -------- サンプル別 *.STRs.tsv 集約 --------
mapfile -t SFILES < <(find "${WORKDIR}" -maxdepth 1 -type f -name "*.STRs.tsv" -size +0 | sort)
[ ${#SFILES[@]} -gt 0 ] || { echo "[ERROR] No per-sample *.STRs.tsv in ${WORKDIR}" | tee -a "${OUTPUT_LOG}"; rm -rf "${LOCAL_TMP}"; exit 1; }

TMP_OUT="${OUTPUT_TSV}.tmp"; : > "${TMP_OUT}"
HEADER_WRITTEN=false
for f in "${SFILES[@]}"; do
  if [ "${HEADER_WRITTEN}" = false ]; then head -n 1 "${f}" > "${TMP_OUT}"; HEADER_WRITTEN=true; fi
  awk -v chr="${CHR}" 'NR>1 && $1==chr' "${f}" >> "${TMP_OUT}"
done

LINES=$(wc -l < "${TMP_OUT}")
[ "${LINES}" -gt 1 ] || { echo "[ERROR] Aggregated header-only: ${TMP_OUT}" | tee -a "${OUTPUT_LOG}"; rm -rf "${LOCAL_TMP}"; exit 1; }
mv "${TMP_OUT}" "${OUTPUT_TSV}"
echo "[INFO] Wrote ${OUTPUT_TSV} (rows: $((LINES-1)))"

# -------- サンプル別の保存（修正版：DEST_PSで統一） --------
DEST_PS="${PER_SAMPLE_DIR}/${CHR}"
mkdir -p "${DEST_PS}"
SAVED=0
for f in "${SFILES[@]}"; do
  base=$(basename "$f"); dst="${DEST_PS}/${base}"
  if [ -e "${dst}" ]; then mv -f "${dst}" "${dst}.bak.$(date +%s)"; fi
  mv -f "$f" "${DEST_PS}/" && SAVED=$((SAVED+1)) || echo "[WARN] move failed: $f" >&2
done
echo "[INFO] Saved ${SAVED} per-sample files into ${DEST_PS}"

# -------- お片付け・サマリー --------
rm -rf "${LOCAL_TMP}"

ELAPSED=$(( $(date +%s) - SCRIPT_START ))
echo "========================================"
echo "Summary"
echo "========================================"
echo "CHR: ${CHR}"
echo "Samples: ${FILTERED_N}"
echo "Out: ${OUTPUT_TSV}"
echo "Rows: $((LINES-1))"
echo "Time: ${ELAPSED}s ($((ELAPSED/60))m)"
echo "========================================"

