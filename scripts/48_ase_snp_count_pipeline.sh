#!/bin/bash
#
# GridEngine resources (tweak as needed)
#$ -pe parallel 32
#$ -l h_vmem=16G
#$ -l h_rt=06:00:00
#$ -j y
#$ -o ./results/ase_mf/logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -S /bin/bash
#$ -cwd
#$ -m n
#$ -t 1-33
#$ -tc 12

set -euo pipefail
THREADS=${NSLOTS:-16}
PROJECT_ROOT="${PROJECT_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}"
CFG_ASE="${CFG_ASE:-${PROJECT_ROOT}/config/ase.yml}"
CFG_UMI="${CFG_UMI:-${PROJECT_ROOT}/config/umi.yml}"
CFG_PATHS="${CFG_PATHS:-${PROJECT_ROOT}/config/paths.yml}"
DESIGN_TSV="${DESIGN_TSV:-${PROJECT_ROOT}/config/design.tsv}"

cfg_get() {
  local cfg="$1"
  local key="$2"
  local default="${3:-}"
  if [[ ! -f "$cfg" ]]; then
    echo "$default"
    return 0
  fi
  local val
  val="$(awk -F': *' -v k="$key" '
    $0 !~ /^[[:space:]]*#/ && $1 == k {
      gsub(/^"/, "", $2); gsub(/"$/, "", $2);
      print $2; exit
    }' "$cfg")"
  if [[ -z "$val" ]]; then
    echo "$default"
  else
    echo "$val"
  fi
}

count_manifest_rows() {
  local f="$1"
  [[ -f "$f" ]] || { echo 0; return 0; }
  awk 'BEGIN{FS="\t"} /^[[:space:]]*($|#)/{next} NF>=3 {c++} END{print c+0}' "$f"
}

###############################
# Allele-specific expression (ASE) with SNPsplit
#
# This script builds an N-masked genome from *fixed differences* between
# your Male (M) and Female (F) DNA samples (joint VCF from nf-core Sarek),
# creates a STAR index, aligns multiple RNA-seq libraries, splits alignments
# by allele with SNPsplit, and produces per-allele gene counts with htseq-count.
#
# NOTE:
#  - Column order in the SNP file is the SNPsplit format (ID CHR POS STRAND M/F)
#  - The allele order is **M/F** (genome1 = M, genome2 = F)
#  - Only positions where M and F are *both homozygous* and *different* are used
#    (i.e., fixed differences), which are best for ASE assignment.
###############################

########### User configuration ###########
# Working directory for this run (dedicated ASE output root inside this project)
WORKDIR=${WORKDIR:-"${PROJECT_ROOT}/results/ase_mf"}
mkdir -p "${WORKDIR}" && cd "${WORKDIR}"
mkdir -p logs

# Reference genome (same strain as M)
genome=${genome:-data/annotation/genome.fa}

# Gene annotation source: default to the same annotation configured for the main pipeline.
# If this is GFF/GFF3, we auto-convert to a GTF compatible with STAR+htseq-count.
gtf=${gtf:-$(cfg_get "$CFG_PATHS" gff_gtf "data/annotation/genes.gff")}

# Joint multi-sample VCF (Sarek output) containing BOTH DNA samples (M and F)
# Must be indexed (.tbi). Update this path to your run.
joint_vcf=${joint_vcf:-$(cfg_get "$CFG_ASE" joint_vcf "data/upstream_sarek/joint_germline.vcf.gz")}

# Sample IDs in the VCF (exact names in the VCF header)
M_SAMPLE=${M_SAMPLE:-$(cfg_get "$CFG_ASE" m_sample "P1_Bm3235M_thl")}
F_SAMPLE=${F_SAMPLE:-$(cfg_get "$CFG_ASE" f_sample "P1_Bm3235F_thl")}

# Optional: exclude contigs (regex). Example filters UV/sex scaffolds if desired.
# Leave empty to keep all contigs. Example: '(^|_)U[Vv]|(^|_)V$'
EXCLUDE_CONTIGS_REGEX="^(ChrVm|ChrUf)$"

# RNA-seq samples manifest (TSV with 3 columns, **no header**):
# sample_id\tR1.fastq.gz\tR2.fastq.gz
# You can reuse your existing paths; one row per library.
SAMPLES_TSV=${SAMPLES_TSV:-$(cfg_get "$CFG_ASE" samples_tsv "${PROJECT_ROOT}/config/ase_samples.tsv")}
SAMPLESHEET=${SAMPLESHEET:-$(cfg_get "$CFG_UMI" samplesheet "data/upstream_rnaseq/pipeline_info/samplesheet.valid.csv")}

# Tool paths. Override these environment variables if tools are not on PATH.
STAR_BIN=${STAR_BIN:-STAR}
SNPSPLIT_BIN=${SNPSPLIT_BIN:-SNPsplit}
TRIMGALORE_IMG=${TRIMGALORE_IMG:-containers/trim-galore.sif}
HTSEQ_IMG=${HTSEQ_IMG:-containers/htseq.sif}

# Directories to bind inside the container (comma-separated)
BIND_PATHS=${BIND_PATHS:-$(pwd)}

# Whether to augment the GTF for htseq-count so UTRs are counted as exons
INCLUDE_UTR_IN_HTSEQ=${INCLUDE_UTR_IN_HTSEQ:-1}

# -------- Variant QC thresholds (for fixed-difference SNP set) --------
MIN_DP=${MIN_DP:-12}           # per-parent min depth
MIN_GQ=${MIN_GQ:-30}           # per-parent min genotype quality
MAX_MINOR_FRAC=${MAX_MINOR_FRAC:-0.05}  # max minor-allele fraction allowed in a homozygote
INDEL_BUFFER=${INDEL_BUFFER:-5}  # drop SNPs within this many bp of any indel
REPEAT_BED=${REPEAT_BED:-}     # optional BED of repeats/low-complexity to exclude

# -------- Mapping bias mitigation (optional) --------
USE_WASP=${USE_WASP:-0}        # 1 = add STAR --varVCFfile + --waspOutputMode SAMtag
WASP_FILTER=${WASP_FILTER:-0}  # 1 = post-STAR, drop reads with vW:i:1 (failed WASP)

########### (Optional) Conda environments ###########
# Optional conda initialization for non-interactive shells.
if [[ -n "${CONDA_SH:-}" && -f "$CONDA_SH" ]]; then
  source "$CONDA_SH"
fi

# bcftools/bedtools/samtools env (edit to your cluster env if you prefer modules)
# conda activate /tmp/global2/alipinska/bcftools/

# Trimming env (as in your original script)

# If no explicit ASE manifest with real rows exists, auto-build from project design + samplesheet.
if [[ "$(count_manifest_rows "$SAMPLES_TSV")" -eq 0 ]]; then
  SAMPLES_TSV="${WORKDIR}/tmp/ase_samples.from_project.tsv"
  mkdir -p "$(dirname "$SAMPLES_TSV")"
  SAMPLES_TSV_TMP="${SAMPLES_TSV}.tmp.$$"
  python - "$SAMPLESHEET" "$DESIGN_TSV" "$SAMPLES_TSV_TMP" <<'PY'
import csv, sys

sheet, design, out = sys.argv[1], sys.argv[2], sys.argv[3]

include = set()
with open(design, newline="", encoding="utf-8") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for r in rdr:
        sid = (r.get("sample_id","") or "").strip()
        inc = (r.get("include","TRUE") or "").strip().upper()
        if sid and inc in {"TRUE","T","1","YES","Y"}:
            include.add(sid)

with open(sheet, newline="", encoding="utf-8-sig") as fh, open(out, "w", encoding="utf-8") as oh:
    rdr = csv.DictReader(fh)
    for r in rdr:
        s = (r.get("sample","") or "").strip()
        fq1 = (r.get("fastq_1","") or "").strip()
        fq2 = (r.get("fastq_2","") or "").strip()
        if not s or not fq1 or not fq2:
            continue
        s_base = s[:-3] if s.endswith("_T1") else s
        if s_base in include:
            oh.write(f"{s}\t{fq1}\t{fq2}\n")
PY
  mv "$SAMPLES_TSV_TMP" "$SAMPLES_TSV"
  echo "[INFO] Auto-generated ASE manifest from config/design.tsv + samplesheet: ${SAMPLES_TSV}" >&2
fi

# Auto-build ASE metadata for the R script if missing.
ASE_META="${ASE_META:-${PROJECT_ROOT}/config/ase_sample_metadata.tsv}"
if [[ ! -s "$ASE_META" || "$(awk 'BEGIN{FS="\t"} NR==1{next} /^[[:space:]]*($|#)/{next} {c++} END{print c+0}' "$ASE_META" 2>/dev/null)" -eq 0 ]]; then
  ASE_META_TMP="${ASE_META}.tmp.$$"
  python - "$DESIGN_TSV" "$ASE_META_TMP" <<'PY'
import csv, sys
design, out = sys.argv[1], sys.argv[2]
with open(design, newline="", encoding="utf-8") as fh, open(out, "w", encoding="utf-8") as oh:
    rdr = csv.DictReader(fh, delimiter="\t")
    oh.write("sample\tgroup\ttissue\treplicate\tploidy\n")
    for r in rdr:
        sid = (r.get("sample_id","") or "").strip()
        inc = (r.get("include","TRUE") or "").strip().upper()
        if not sid or inc not in {"TRUE","T","1","YES","Y"}:
            continue
        sample = sid if sid.endswith("_T1") else f"{sid}_T1"
        group = (r.get("stage","") or "").strip() or (r.get("tissue","") or "").strip() or "group1"
        tissue = (r.get("tissue","") or "").strip() or group
        rep = (r.get("replicate","") or "").strip() or "1"
        gen = (r.get("generation","") or "").strip().lower()
        if gen == "1n":
            ploidy = "1"
        elif gen == "2n":
            ploidy = "2"
        else:
            ploidy = ""
        oh.write(f"{sample}\t{group}\t{tissue}\t{rep}\t{ploidy}\n")
PY
  mv "$ASE_META_TMP" "$ASE_META"
  echo "[INFO] Auto-generated ASE metadata from config/design.tsv: ${ASE_META}" >&2
fi

##########################################
# 0) Folders
##########################################
mkdir -p logs vcf snps masked_genome star_index trimmed align snpsplit counts tmp

# Normalize annotation: if input is GFF/GFF3, convert to GTF with gene_id/transcript_id.
ANNOT_GTF_PATH="${gtf}"
case "${gtf,,}" in
  *.gff|*.gff3)
    gtf_base="$(basename "${gtf}")"
    gtf_base="${gtf_base%.gff}"
    gtf_base="${gtf_base%.gff3}"
    ANNOT_GTF_PATH="tmp/${gtf_base}.from_main_gff.gtf"
    if [[ ! -s "${ANNOT_GTF_PATH}" || "${gtf}" -nt "${ANNOT_GTF_PATH}" ]]; then
      echo "[INFO] Converting GFF -> GTF for STAR/htseq: ${gtf} -> ${ANNOT_GTF_PATH}" >&2
      awk -F'\t' '
        function trim(s){gsub(/^[[:space:]]+|[[:space:]]+$/,"",s); return s}
        function parse_attr(attr, key,   n,i,p,kv,k,v){
          n=split(attr,p,";");
          for(i=1;i<=n;i++){
            kv=trim(p[i]); if(kv=="") continue;
            split(kv, v, "="); k=trim(v[1]);
            if(k==key){sub(/^[^=]*=/,"",kv); return trim(kv)}
          }
          return ""
        }
        BEGIN{OFS="\t"}
        NR==FNR{
          if($0 ~ /^#/) next;
          id=parse_attr($9,"ID");
          parent=parse_attr($9,"Parent");
          biotype=parse_attr($9,"biotype");
          f=tolower($3);
          if(f=="gene" && id!=""){ gene_biotype[id]=biotype; next }
          if((f=="mrna" || f=="transcript") && id!=""){
            tx2gene[id]=parent;
            if(parent in gene_biotype) tx_biotype[id]=gene_biotype[parent];
            else tx_biotype[id]=biotype;
          }
          next;
        }
        {
          if($0 ~ /^#/) next;
          id=parse_attr($9,"ID");
          parent=parse_attr($9,"Parent");
          split(parent, pp, ","); parent=pp[1];
          f=tolower($3);
          gene_id=""; tx_id="";
          if(f=="gene"){
            gene_id=id; tx_id=id".t1";
          } else if(f=="mrna" || f=="transcript"){
            tx_id=id; gene_id=parent;
          } else {
            tx_id=parent;
            gene_id=tx2gene[tx_id];
            if(gene_id==""){
              gene_id=tx_id;
              sub(/\.t[0-9]+$/, "", gene_id);
            }
          }
          if(gene_id=="" || tx_id=="") next;
          attrs="gene_id \"" gene_id "\"; transcript_id \"" tx_id "\";";
          print $1,$2,$3,$4,$5,$6,$7,$8,attrs;
        }
      ' "${gtf}" "${gtf}" > "${ANNOT_GTF_PATH}.tmp.$$"
      mv "${ANNOT_GTF_PATH}.tmp.$$" "${ANNOT_GTF_PATH}"
    fi
    ;;
esac

# Prepare a GTF for htseq-count that treats UTRs as exons (helps retain informative SNPs in UTRs)
HTSEQ_GTF_PATH="${ANNOT_GTF_PATH}"
if [[ "${INCLUDE_UTR_IN_HTSEQ}" -eq 1 ]]; then
  gtf_base="$(basename "${ANNOT_GTF_PATH}")"
  gtf_base="${gtf_base%.gtf}"
  gtf_base="${gtf_base%.gff}"
  gtf_base="${gtf_base%.gff3}"
  HTSEQ_GTF_PATH="tmp/${gtf_base}.htseq_utr_exons.gtf"
  HTSEQ_GTF_LOCK="${HTSEQ_GTF_PATH}.lock"
  # Build once under lock to avoid array-task races on shared tmp path.
  while ! mkdir "${HTSEQ_GTF_LOCK}" 2>/dev/null; do
    if [[ -s "${HTSEQ_GTF_PATH}" && "${HTSEQ_GTF_PATH}" -nt "${ANNOT_GTF_PATH}" ]]; then
      break
    fi
    sleep 1
  done
  if [[ -d "${HTSEQ_GTF_LOCK}" ]]; then
    if [[ ! -s "${HTSEQ_GTF_PATH}" || "${ANNOT_GTF_PATH}" -nt "${HTSEQ_GTF_PATH}" ]]; then
      echo "[INFO] Building htseq-count GTF with duplicated UTR -> exon entries: ${HTSEQ_GTF_PATH}" >&2
      tmp_build="${HTSEQ_GTF_PATH}.tmp.$$"
      awk 'BEGIN{OFS="\t"}
           /^#/ {print; next}
           {
             print $0;
             f=tolower($3);
             if (f=="five_prime_utr" || f=="three_prime_utr" || f=="utr") {
               $3="exon";
               print;
             }
           }' "${ANNOT_GTF_PATH}" > "${tmp_build}"
      mv "${tmp_build}" "${HTSEQ_GTF_PATH}"
    fi
    rmdir "${HTSEQ_GTF_LOCK}" || true
  fi
else
  echo "[INFO] INCLUDE_UTR_IN_HTSEQ=0; using provided GTF for htseq-count without modification." >&2
fi

if [[ -z "${SGE_TASK_ID:-}" || "${SGE_TASK_ID}" -eq 1 ]]; then
  ##########################################
  # Validate that M_SAMPLE and F_SAMPLE exist in the VCF header
  ##########################################
  mapfile -t __VCF_SAMPLES < <(bcftools query -l "${joint_vcf}")
  if ! printf '%s\n' "${__VCF_SAMPLES[@]}" | grep -Fxq "${M_SAMPLE}"; then
    echo "ERROR: M_SAMPLE ('${M_SAMPLE}') not found in VCF: ${joint_vcf}" >&2
    echo "Available samples:" >&2
    printf '  - %s\n' "${__VCF_SAMPLES[@]}" >&2
    exit 2
  fi
  if ! printf '%s\n' "${__VCF_SAMPLES[@]}" | grep -Fxq "${F_SAMPLE}"; then
    echo "ERROR: F_SAMPLE ('${F_SAMPLE}') not found in VCF: ${joint_vcf}" >&2
    echo "Available samples:" >&2
    printf '  - %s\n' "${__VCF_SAMPLES[@]}" >&2
    exit 2
  fi

  ##########################################
  # 1) Derive fixed M vs F SNP set (biallelic, hom/hom, different)
  ##########################################
  # Keep biallelic SNPs for the two samples of interest, drop SNPs near indels, left-normalize
  if [[ ! -s vcf/MF.biallelic.snps.vcf.gz ]]; then
    # NOTE: do NOT require -f PASS here; Sarek joint VCF often uses '.' in FILTER.
    bcftools view -m2 -M2 -v snps -s "${M_SAMPLE},${F_SAMPLE}" -Ou "${joint_vcf}" \
    | bcftools filter --SnpGap ${INDEL_BUFFER} -Ou \
    | bcftools norm -f "${genome}" -m -both -Oz -o vcf/MF.biallelic.snps.vcf.gz
    bcftools index -ft vcf/MF.biallelic.snps.vcf.gz
  else
    echo "[SKIP] vcf/MF.biallelic.snps.vcf.gz already exists."
  fi

  # Optionally remove SNPs overlapping repeats/low-complexity
  VCF_PREFILT=vcf/MF.biallelic.snps.vcf.gz
  if [[ -n "${REPEAT_BED}" ]]; then
    if [[ ! -s vcf/MF.biallelic.snps.norep.vcf.gz ]]; then
      bedtools intersect -header -v -a vcf/MF.biallelic.snps.vcf.gz -b "${REPEAT_BED}" \
        | bgzip -c > vcf/MF.biallelic.snps.norep.vcf.gz
      tabix -f vcf/MF.biallelic.snps.norep.vcf.gz
    fi
    VCF_PREFILT=vcf/MF.biallelic.snps.norep.vcf.gz
  fi

  # Extract GT/AD/DP/GQ and compute fixed differences (M hom != F hom) with purity checks
  if [[ ! -s snps/MF.fixed_diffs.table.tsv || ! -s snps/MF.fixed_diffs.snpfile.txt ]]; then
    bcftools annotate -x INFO,^FMT/GT,FMT/AD,FMT/DP,FMT/GQ -Ou "${VCF_PREFILT}" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%AD:%DP:%GQ]\n' -s "${M_SAMPLE},${F_SAMPLE}" \
    | awk -v exc_re="${EXCLUDE_CONTIGS_REGEX}" -v m="${M_SAMPLE}" -v f="${F_SAMPLE}" \
           -v mindp="${MIN_DP}" -v mingq="${MIN_GQ}" -v maxminor="${MAX_MINOR_FRAC}" 'BEGIN{FS=OFS="\t"}
      function ok_hom(gt, ad, dp, gq,  a1,a2,ref,alt,refad,altad,minor){
        if(dp=="."||gq==".") return 0;
        if(dp+0 < mindp || gq+0 < mingq) return 0;
        # parse AD as ref,alt
        n=split(ad, aa, ","); if(n<2) return 0; refad=aa[1]+0; altad=aa[2]+0;
        if(gt=="0/0"||gt=="0|0"){ minor = (dp>0)? altad/dp : 1; return (minor <= maxminor) }
        else if(gt=="1/1"||gt=="1|1"){ minor = (dp>0)? refad/dp : 1; return (minor <= maxminor) }
        else return 0;
      }
      {
        chr=$1; pos=$2; ref=$3; alt=$4; mgt=""; mad=""; mdp=""; mgq=""; fgt=""; fad=""; fdp=""; fgq="";
        for(i=5;i<=NF;i++){
          split($i,b, "="); if(length(b)!=2) continue;
          split(b[2],c, ":"); if(length(c)<4) continue;
          if(b[1]==m){ mgt=c[1]; mad=c[2]; mdp=c[3]; mgq=c[4]; }
          else if(b[1]==f){ fgt=c[1]; fad=c[2]; fdp=c[3]; fgq=c[4]; }
        }
        if(exc_re!="" && chr ~ exc_re) next;
        # normalize GT phasing so 1|1 == 1/1 and 0|0 == 0/0 before comparisons
        gsub(/\|/,"/",mgt); gsub(/\|/,"/",fgt);
        if(!ok_hom(mgt,mad,mdp,mgq) || !ok_hom(fgt,fad,fdp,fgq)) next;
        if(mgt==fgt) next; # must differ after normalization
        mallele = (mgt ~ /^0\/0$/) ? ref : alt;
        fallele = (fgt ~ /^0\/0$/) ? ref : alt;
        print chr, pos, ref, alt, mallele, fallele;
      }' > snps/MF.fixed_diffs.table.tsv

    awk 'BEGIN{FS=OFS="\t"}{print NR,$1,$2,1,$5"/"$6}' snps/MF.fixed_diffs.table.tsv > snps/MF.fixed_diffs.snpfile.txt
  else
    echo "[SKIP] snps/MF.fixed_diffs.table.tsv and snps/MF.fixed_diffs.snpfile.txt exist."
  fi

  # Create a positions list and a reduced VCF containing just those SNPs
  if [[ ! -s vcf/MF.fixed_diffs.vcf.gz ]]; then
    cut -f1,2 snps/MF.fixed_diffs.table.tsv > snps/MF.fixed_diffs.positions.tsv
    if [[ ! -s snps/MF.fixed_diffs.positions.tsv ]] || [[ $(wc -l < snps/MF.fixed_diffs.positions.tsv) -eq 0 ]]; then
      echo "ERROR: No fixed M↔F sites passed the current thresholds (positions file is empty)." >&2
      echo "Hint: relax MIN_DP/MIN_GQ or skip AD-based purity by setting MAX_MINOR_FRAC=1, or check sample IDs in ${joint_vcf}." >&2
      exit 3
    fi
    bcftools view -R snps/MF.fixed_diffs.positions.tsv -Oz -o vcf/MF.fixed_diffs.vcf.gz "${VCF_PREFILT}"
    bcftools index -ft vcf/MF.fixed_diffs.vcf.gz
  else
    echo "[SKIP] vcf/MF.fixed_diffs.vcf.gz already exists."
  fi

  ##########################################
  # 2) N-mask the reference at these SNPs
  ##########################################
  # Sort the VCF for bedtools and mask
  if [[ ! -s vcf/MF.fixed_diffs.sorted.vcf.gz ]]; then
    bcftools sort -Oz -o vcf/MF.fixed_diffs.sorted.vcf.gz vcf/MF.fixed_diffs.vcf.gz
  fi
  # Create/refresh an uncompressed varVCF for STAR WASP (avoid per-task races)
  if [[ ! -s vcf/MF.fixed_diffs.sorted.vcf || vcf/MF.fixed_diffs.sorted.vcf -ot vcf/MF.fixed_diffs.sorted.vcf.gz ]]; then
    gzip -cd vcf/MF.fixed_diffs.sorted.vcf.gz > vcf/MF.fixed_diffs.sorted.vcf
  fi
  if [[ ! -s masked_genome/Bm3235.masked_MFfixed.fa || ! -s masked_genome/Bm3235.masked_MFfixed.fa.fai ]]; then
    bedtools maskfasta -fi "${genome}" -bed vcf/MF.fixed_diffs.sorted.vcf.gz -fo masked_genome/Bm3235.masked_MFfixed.fa
    samtools faidx masked_genome/Bm3235.masked_MFfixed.fa
  else
    echo "[SKIP] masked genome already exists."
  fi

  ##########################################
  # 3) Build STAR index on N-masked genome
  ##########################################
  if [[ ! -s star_index/Bm3235_masked_MFfixed/SA ]]; then
    mkdir -p star_index/Bm3235_masked_MFfixed
    "${STAR_BIN}" \
      --runThreadN ${THREADS} \
      --runMode genomeGenerate \
      --genomeDir star_index/Bm3235_masked_MFfixed \
      --genomeFastaFiles masked_genome/Bm3235.masked_MFfixed.fa \
      --genomeSAindexNbases 12 \
      --sjdbGTFfile "${ANNOT_GTF_PATH}"
  else
    echo "[SKIP] STAR index already present."
  fi

fi

# Wait for prep artifacts (when running as an array task not equal to 1)
for f in snps/MF.fixed_diffs.snpfile.txt masked_genome/Bm3235.masked_MFfixed.fa star_index/Bm3235_masked_MFfixed/SA; do
  until [ -s "$f" ]; do
    echo "[TASK ${SGE_TASK_ID:-NA}] Waiting for $f ..." >&2
    sleep 60
  done
done

# If WASP is enabled, also wait for the uncompressed varVCF to be ready
if [[ "${USE_WASP}" -eq 1 ]]; then
  until [ -s "vcf/MF.fixed_diffs.sorted.vcf" ]; do
    echo "[TASK ${SGE_TASK_ID:-NA}] Waiting for vcf/MF.fixed_diffs.sorted.vcf ..." >&2
    sleep 30
  done
fi

# ---------------- Manifest sanitization (robust TSV reader) ----------------
# Create a cleaned manifest to avoid CRLF/whitespace/blank-lines issues that often
# affect the last line. We also ensure there is a trailing newline.
MANIFEST_CLEAN=tmp/ase_samples.clean.tsv
if [[ ! -s "$MANIFEST_CLEAN" || "$SAMPLES_TSV" -nt "$MANIFEST_CLEAN" ]]; then
  mkdir -p tmp
  MANIFEST_CLEAN_TMP="${MANIFEST_CLEAN}.tmp.$$"
  # Strip CR, drop blank/comment lines, trim fields, keep first 3 tab-separated fields only.
  # NOTE: Avoid awkward "i in $" awk construct; use NF and explicit trimming.
  tr -d '\r' < "$SAMPLES_TSV" \
  | awk 'BEGIN{FS="\t"; OFS="\t"}
          /^[[:space:]]*($|#)/{next}
          {
            # ensure at least 3 fields exist
            for(i=1;i<=3;i++){ if(i>NF) $i="" }
            # trim whitespace around each field
            for(i=1;i<=3;i++){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i) }
            if($1!="" && $2!="" && $3!="") print $1,$2,$3;
          }' > "$MANIFEST_CLEAN_TMP"
  # ensure newline at EOF
  tail -c1 "$MANIFEST_CLEAN_TMP" | read -r _ || echo >> "$MANIFEST_CLEAN_TMP"
  mv "$MANIFEST_CLEAN_TMP" "$MANIFEST_CLEAN"
fi

N_MANIFEST=$(wc -l < "$MANIFEST_CLEAN" | tr -d ' ')
echo "[INFO] Using manifest: $MANIFEST_CLEAN" >&2

##########################################
# 4) Process RNA-seq libraries (trim -> map -> split -> count)
##########################################
process_sample() {
  local sample="$1" r1="$2" r2="$3"
  [ -z "$sample" ] && return 0
  echo "[START] ${sample}"

  # Define output directories up front
  local sdir_trim=trimmed/${sample}
  local sdir_aln=align/${sample}
  local sdir_split=snpsplit/${sample}
  local sdir_counts=counts/${sample}

  # If a previous run completed fully, skip early without touching the lock
  if [[ -f "${sdir_counts}/.done" ]]; then
    echo "[SKIP] ${sample}: sample marked complete (${sdir_counts}/.done)."
    return 0
  fi

  # Sanity check: inputs exist
  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "ERROR: FASTQs not found for ${sample}:" >&2
    echo "  R1: ${r1}" >&2
    echo "  R2: ${r2}" >&2
    return 5
  fi


  # Per-sample lock to avoid duplicate work; remove stale locks
  local lock="tmp/${sample}.lock"
  if [[ -e "${lock}" ]]; then
    read -r old_pid old_job < "${lock}" 2>/dev/null || true
    if [[ -n "${old_pid:-}" ]] && ps -p "${old_pid}" >/dev/null 2>&1; then
      echo "[SKIP] ${sample}: another process holds lock (pid=${old_pid}, job=${old_job})." >&2
      return 0
    else
      echo "[WARN] ${sample}: removing stale lock ${lock}" >&2
      rm -f "${lock}"
    fi
  fi
  if ( set -o noclobber; printf '%s %s\n' "$$" "${JOB_ID}.${SGE_TASK_ID:-NA}" > "${lock}" ) 2>/dev/null; then
    trap "rm -f '${lock}'" EXIT
  else
    echo "[SKIP] ${sample}: could not acquire lock." >&2
    return 0
  fi

  # Ensure output directories are writable (helps with container outdir errors)
  for d in "${sdir_trim}" "${sdir_aln}" "${sdir_split}" "${sdir_counts}"; do
    mkdir -p "$d"
    if [[ ! -w "$d" ]]; then
      echo "ERROR: Output directory is not writable: $d" >&2
      return 4
    fi
  done

  ################ Trim adapters ################
  if [[ -f "${sdir_trim}/.ok" ]]; then
    echo "[SKIP] ${sample}: trimming already OK."
  else
    # Clean partials if present without .ok
    if compgen -G "${sdir_trim}/*val_*.fq" > /dev/null; then
      echo "[CLEAN] ${sample}: removing partial trimmed reads (no .ok)" >&2
      rm -f "${sdir_trim}"/*val_*.fq || true
    fi
    apptainer exec --cleanenv --bind "${BIND_PATHS}" "${TRIMGALORE_IMG}" \
      trim_galore --dont_gzip --stringency 1 --fastqc --length 5 --paired \
      --output_dir "${sdir_trim}" "${r1}" "${r2}" |& tee "${sdir_trim}/trim_galore.log"
    # Verify outputs
    if ! compgen -G "${sdir_trim}/*val_1.fq" > /dev/null || ! compgen -G "${sdir_trim}/*val_2.fq" > /dev/null; then
      echo "ERROR: trimming did not produce paired *_val_1.fq and *_val_2.fq for ${sample}" >&2
      return 6
    fi
    touch "${sdir_trim}/.ok"
  fi

  local fq1
  fq1=$(ls "${sdir_trim}"/*val_1.fq | paste -sd "," -)
  local fq2
  fq2=$(ls "${sdir_trim}"/*val_2.fq | paste -sd "," -)
  if [[ -z "${fq1}" || -z "${fq2}" ]]; then
    echo "ERROR: Did not find trimmed FASTQs for ${sample} after trimming (.ok present)." >&2
    ls -l "${sdir_trim}" >&2 || true
    return 6
  fi

  ################ STAR alignment ################
  local bam="${sdir_aln}/${sample}.Aligned.sortedByCoord.out.bam"
  if [[ -f "${sdir_aln}/.ok" ]]; then
    echo "[SKIP] ${sample}: STAR alignment already OK."
  else
    # Use a single, pre-generated uncompressed VCF for STAR WASP to avoid race conditions across array tasks
    local VAR_VCF_TXT="vcf/MF.fixed_diffs.sorted.vcf"
    local NVARS=0
    if [[ -s "$VAR_VCF_TXT" ]]; then
      NVARS=$(awk '!/^#/ {c++} END{print c+0}' "$VAR_VCF_TXT" 2>/dev/null || echo 0)
    fi

    STAR_ARGS=(
      --genomeDir star_index/Bm3235_masked_MFfixed
      --outFileNamePrefix "${sdir_aln}/${sample}."
      --readFilesIn "${fq1}" "${fq2}"
      --runThreadN ${THREADS}
      --quantMode GeneCounts
      --alignEndsType EndToEnd
      --outSAMattributes NH HI NM MD
      --outSAMtype BAM SortedByCoordinate
    )
    if [[ "${USE_WASP}" -eq 1 ]]; then
      if [[ "$NVARS" -gt 0 ]]; then
        STAR_ARGS+=( --varVCFfile "$VAR_VCF_TXT" --waspOutputMode SAMtag )
      else
        echo "[INFO] No variants found in $VAR_VCF_TXT; disabling WASP for ${sample}." >&2
      fi
    fi
    echo "[INFO] ${sample}: launching STAR (USE_WASP=${USE_WASP}, NVARS=${NVARS})" >&2
    "${STAR_BIN}" "${STAR_ARGS[@]}"
    samtools index "${bam}"
    if [[ "${WASP_FILTER}" -eq 1 ]]; then
      echo "[INFO] Applying WASP read filter (dropping vW:i:1)" >&2
      samtools view -h "${bam}" \
        | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {fail=0; for(i=12;i<=NF;i++){ if($i=="vW:i:1"){fail=1; break} } if(!fail) print }' \
        | samtools view -b -o "${bam%.bam}.wasp.bam" -
      mv -f "${bam%.bam}.wasp.bam" "${bam}"
      samtools index -f "${bam}"
    fi
    # Verify BAM and index, then mark OK
    if [[ ! -s "${bam}" || ! -s "${bam}.bai" ]]; then
      echo "ERROR: STAR did not produce indexed BAM for ${sample}" >&2
      return 7
    fi
    touch "${sdir_aln}/.ok"
  fi

  ################ SNPsplit ################
  if [[ -f "${sdir_split}/.ok" ]]; then
    echo "[SKIP] ${sample}: SNPsplit already OK."
    g1_bam=$(ls "${sdir_split}"/*genome1.bam 2>/dev/null | head -1 || true)
    g2_bam=$(ls "${sdir_split}"/*genome2.bam 2>/dev/null | head -1 || true)
    if [[ -z "${g1_bam}" || -z "${g2_bam}" ]] || [[ ! -s "${g1_bam}" || ! -s "${g2_bam}" ]]; then
      echo "ERROR: ${sample}: SNPsplit marked OK but genome1/genome2 BAMs are missing; removing stale marker." >&2
      rm -f "${sdir_split}/.ok"
      return 8
    fi
  else
    # Clean partials if present without .ok
    if compgen -G "${sdir_split}/*genome*.bam" > /dev/null; then
      echo "[CLEAN] ${sample}: removing partial SNPsplit BAMs (no .ok)" >&2
      rm -f "${sdir_split}"/*genome*.bam || true
    fi
    "${SNPSPLIT_BIN}" \
      --snp_file snps/MF.fixed_diffs.snpfile.txt \
      --paired \
      --no_sort \
      --output_dir "${sdir_split}" \
      "${bam}"
    # Verify outputs
    g1_bam=$(ls "${sdir_split}"/*genome1.bam 2>/dev/null | head -1 || true)
    g2_bam=$(ls "${sdir_split}"/*genome2.bam 2>/dev/null | head -1 || true)
    if [[ -z "${g1_bam}" || -z "${g2_bam}" ]] || [[ ! -s "${g1_bam}" || ! -s "${g2_bam}" ]]; then
      echo "ERROR: SNPsplit did not produce genome1/genome2 BAMs for ${sample}" >&2
      return 8
    fi
    touch "${sdir_split}/.ok"
  fi

  ################ htseq-count ################
  # Adjust strandedness (-s) if your libraries are stranded
  if [[ -f "${sdir_counts}/.ok" ]]; then
    echo "[SKIP] ${sample}: htseq-count already OK."
  else
    # Use the BAMs found/verified in SNPsplit step
    apptainer exec --cleanenv --bind "${BIND_PATHS}" "${HTSEQ_IMG}" \
      htseq-count -f bam -r pos -s no -i gene_id "${g1_bam}" "${HTSEQ_GTF_PATH}" > "${sdir_counts}/${sample}.M.counts.txt"
    apptainer exec --cleanenv --bind "${BIND_PATHS}" "${HTSEQ_IMG}" \
      htseq-count -f bam -r pos -s no -i gene_id "${g2_bam}" "${HTSEQ_GTF_PATH}" > "${sdir_counts}/${sample}.F.counts.txt"
    # Verify outputs and mark OK
    if [[ ! -s "${sdir_counts}/${sample}.M.counts.txt" || ! -s "${sdir_counts}/${sample}.F.counts.txt" ]]; then
      echo "ERROR: htseq-count did not produce both M and F counts for ${sample}" >&2
      return 9
    fi
    touch "${sdir_counts}/.ok"
  fi

  # Mark the whole sample as done only after all steps succeeded
  touch "${sdir_counts}/.done"

  echo "[DONE] ${sample}"
}

if [[ -n "${SGE_TASK_ID:-}" ]]; then
  # Array mode: pick the N-th line from the cleaned manifest
  if [[ ${SGE_TASK_ID} -gt ${N_MANIFEST} ]]; then
    echo "[SKIP] Task ${SGE_TASK_ID} > manifest lines (${N_MANIFEST}); nothing to do." >&2
    exit 0
  fi
  IFS=$'\t' read -r sample r1 r2 < <(sed -n "${SGE_TASK_ID}p" "$MANIFEST_CLEAN")
  if [[ -z "${sample:-}" ]]; then
    echo "No sample found at line ${SGE_TASK_ID} of $MANIFEST_CLEAN. Exiting." >&2
    exit 0
  fi
  process_sample "${sample}" "${r1}" "${r2}"
else
  # Fallback: process all samples sequentially from the cleaned manifest
  while IFS=$'\t' read -r sample r1 r2; do
    [[ -z "${sample}" || -z "${r1}" || -z "${r2}" ]] && continue
    process_sample "${sample}" "${r1}" "${r2}"
  done < "$MANIFEST_CLEAN"
fi

##########################################
# 5) Run summaries (only in non-array mode)
##########################################
if [[ -z "${SGE_TASK_ID:-}" ]]; then
  # Record the allele order used by SNPsplit
  printf "genome1\tM\ngenome2\tF\n" > snps/MF.genome_map.txt
  # Collect a simple inventory of outputs
  find counts -type f -name "*.counts.txt" | sort > logs/count_files.list
  echo "All done. Outputs in: ${WORKDIR}"
fi
