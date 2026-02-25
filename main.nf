nextflow.enable.dsl=2

/*
  EPI2ME refreshes/inspects workflows without user params.
  So we MUST NOT fail if params.reads / params.db_fasta are not set yet.
*/

workflow {

  if( !params.reads || !params.db_fasta ) {
    log.info "Parameters --reads and/or --db_fasta not set yet (UI refresh). Skipping execution."
    return
  }

  // Input: READS_ROOT/Sample/*fastq.gz
  Channel
    .fromPath("${params.reads}/*/*.fastq.gz", checkIfExists: true)
    .map { f -> tuple(f.parent.name, f) }
    .groupTuple()
    .map { sample, files -> tuple(sample, files.sort()) }
    .set { ch_samples }

  merged   = MERGE_FASTQ(ch_samples)
  qc_raw   = QC_SEQKIT(merged)

  trimmed  = (params.primer_fwd?.trim() && params.primer_rev?.trim()) ? TRIM_CUTADAPT(merged) : merged
  filtered = FILTER_NANOFILT(trimmed)

  fasta     = FASTQ_TO_FASTA(filtered)
  clustered = CLUSTER_VSEARCH(fasta)
  kept      = FILTER_CLUSTERS(clustered)

  // Build BLAST DB once (from params.db_fasta) and attach to every sample tuple
  ref_fa   = Channel.value( file(params.db_fasta) )
  dbdir_ch = MAKEBLASTDB(ref_fa)

  kept_with_db = kept.combine(dbdir_ch)
  blasted      = BLASTN(kept_with_db)

  tax = JOIN_COUNTS_BLAST(blasted)   // emits (sample, sample.taxonomy.tsv)

  // Aggregate all sample tax tables
  tax_tables_ch   = tax.map { sample, taxfile -> taxfile }
  tax_tables_list = tax_tables_ch.collect()   // emits ONE item: a list of all taxonomy files

  summary = AGGREGATE_RESULTS(tax_tables_list)

  // summary.out ist das Output-Tuple: (all_samples.long.tsv, abundance_matrix.tsv)
  tsv_ch = summary
    .map { long_tsv, mat_tsv -> [ long_tsv, mat_tsv ] }
    .flatten()

  // TSV -> CSV (ein Prozess-Aufruf, zwei Dateien rein)
  csv_ch = TSV_TO_CSV(tsv_ch)

  // beide CSVs sammeln und in ein Tuple für REPORT_HTML packen
  csv_list = csv_ch.collect()
  csv_tuple = csv_list.map { files ->
      def m = files.collectEntries { [(it.name): it] }
      tuple(m["all_samples.long.csv"], m["abundance_matrix.csv"])
  }

  REPORT_HTML(csv_tuple)

}


// --------------------
// Processes
// --------------------

process MERGE_FASTQ {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("${sample}.fastq.gz")

  """
  # Concatenating gz streams is valid
  cat ${reads.join(' ')} > ${sample}.fastq.gz
  """
}

process QC_SEQKIT {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/qc", mode: 'copy'
  container "quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0"

  input:
  tuple val(sample), path(fq)

  output:
  tuple val(sample), path("${sample}.seqkit.stats.tsv")

  """
  seqkit stats -a -T ${fq} > ${sample}.seqkit.stats.tsv
  """
}

process TRIM_CUTADAPT {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'
  container "python:3.11-slim"

  input:
  tuple val(sample), path(fq)

  output:
  tuple val(sample), path("${sample}.trimmed.fastq.gz")

  def discard = params.require_primers ? "--discard-untrimmed" : ""

  """
  python -m pip install --no-cache-dir cutadapt==4.9 >/dev/null

  cutadapt \
    --revcomp \
    -e ${params.trim_error_rate} \
    -g ${params.primer_fwd} \
    -a ${params.primer_rev} \
    ${discard} \
    -o ${sample}.trimmed.fastq.gz \
    ${fq}
  """
}

process FILTER_NANOFILT {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'
  container "python:3.11-slim"

  input:
  tuple val(sample), path(fq)

  output:
  tuple val(sample), path("${sample}.filtered.fastq.gz")

  """
  python - << 'PY'
  import gzip
  import sys

  infile  = "${fq}"
  outfile = "${sample}.filtered.fastq.gz"

  min_q   = int(${params.min_q})
  min_len = int(${params.min_len})
  max_len = int(${params.max_len})

  def open_in(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

  kept = 0
  total = 0

  with open_in(infile) as fin, gzip.open(outfile, "wt") as fout:
    while True:
      h = fin.readline()
      if not h:
        break
      seq  = fin.readline().rstrip("\\n")
      plus = fin.readline()
      qual = fin.readline().rstrip("\\n")
      total += 1

      L = len(seq)
      if L < min_len or L > max_len:
        continue
      if len(qual) != L:
        continue

      # mean Q from Phred+33
      qsum = 0
      for c in qual:
        qsum += (ord(c) - 33)
      mean_q = qsum / L if L else 0.0

      if mean_q < min_q:
        continue

      fout.write(h)
      fout.write(seq + "\\n")
      fout.write(plus)
      fout.write(qual + "\\n")
      kept += 1

  sys.stderr.write(f"FILTER: kept {kept}/{total} reads\\n")
  PY
  """
}

process FASTQ_TO_FASTA {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'
  container "python:3.11-slim"

  input:
  tuple val(sample), path(fq)

  output:
  tuple val(sample), path("${sample}.fasta")

  """
  python - << 'PY'
  import gzip

  infile  = "${fq}"
  outfile = "${sample}.fasta"

  def open_maybe_gzip(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

  with open_maybe_gzip(infile, "rt") as fin, open(outfile, "wt") as fout:
    while True:
      h = fin.readline()
      if not h:
        break
      seq = fin.readline().strip()
      fin.readline()  # +
      fin.readline()  # qual
      if h.startswith("@"):
        h = h[1:]
      fout.write(">" + h.strip() + "\\n")
      fout.write(seq + "\\n")
  PY
  """
}

process CLUSTER_VSEARCH {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'
  container "quay.io/biocontainers/vsearch:2.28.1--h6a68c12_0"

  input:
  tuple val(sample), path(fa)

  output:
  tuple val(sample),
        path("${sample}.centroids.fasta"),
        path("${sample}.clusters.uc")

  """
  vsearch \
    --cluster_fast ${fa} \
    --id ${params.cluster_id} \
    --strand both \
    --threads ${task.cpus} \
    --uc ${sample}.clusters.uc \
    --centroids ${sample}.centroids.fasta
  """
}

process FILTER_CLUSTERS {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'
  container "quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0"

  input:
  tuple val(sample), path(centroids), path(uc)

  output:
  tuple val(sample),
        path("${sample}.cluster_counts.tsv"),
        path("${sample}.centroids.kept.fasta")

  """
  awk -F'\\t' '
    BEGIN{OFS="\\t"}
    \$1=="S"{cl=\$2; id=\$9; centroid[cl]=id; count[id]=1; next}
    \$1=="H"{cl=\$2; id=centroid[cl]; count[id]++; next}
    END{for (id in count) print id, count[id]}
  ' ${uc} | sort -k2,2nr > ${sample}.cluster_counts.tsv

  awk -v m=${params.min_cluster_reads} '\$2>=m{print \$1}' ${sample}.cluster_counts.tsv > keep_ids.txt || true

  if [ -s keep_ids.txt ]; then
    seqkit grep -f keep_ids.txt ${centroids} > ${sample}.centroids.kept.fasta
  else
    : > ${sample}.centroids.kept.fasta
  fi
  """
}

process MAKEBLASTDB {
  publishDir "${params.out_dir}/refdb", mode: 'copy'
  container "ncbi/blast:2.16.0"

  input:
  path(db_fasta)

  output:
  path("blastdb")

  """
  mkdir -p blastdb
  cp ${db_fasta} blastdb/db.fasta
  makeblastdb -in blastdb/db.fasta -dbtype nucl -out blastdb/fusarium_tef1
  """
}

process BLASTN {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/taxonomy", mode: 'copy'
  container "ncbi/blast:2.16.0"

  input:
  tuple val(sample), path(counts), path(centroids_fa), path(dbdir)

  output:
  tuple val(sample), path(counts), path("${sample}.blast.tsv")

  """
  if [ ! -s ${centroids_fa} ]; then
    : > ${sample}.blast.tsv
    exit 0
  fi

  blastn \
    -query ${centroids_fa} \
    -db ${dbdir}/fusarium_tef1 \
    -max_target_seqs ${params.blast_topn} \
    -num_threads ${task.cpus} \
    -outfmt "6 qseqid sseqid pident length qlen bitscore evalue qcovs stitle" \
    > ${sample}.blast.tsv
  """
}

process JOIN_COUNTS_BLAST {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/taxonomy", mode: 'copy'
  container "python:3.11-slim"

  input:
  tuple val(sample), path(counts_tsv), path(blast_tsv)

  output:
  tuple val(sample), path("${sample}.taxonomy.tsv")

  """
  python - << 'PY'
  import csv
  from pathlib import Path

  sample = "${sample}"
  counts_path = Path("${counts_tsv}")
  blast_path  = Path("${blast_tsv}")
  out_path    = Path(f"{sample}.taxonomy.tsv")

  best = {}
  if blast_path.exists() and blast_path.stat().st_size > 0:
    with blast_path.open() as f:
      for row in csv.reader(f, delimiter='\\t'):
        q = row[0]
        if q in best:
          continue
        best[q] = {"sseqid": row[1], "pident": row[2], "qcovs": row[7], "stitle": row[8] if len(row)>8 else ""}

  with out_path.open("w", newline="") as out:
    w = csv.writer(out, delimiter='\\t')
    w.writerow(["cluster_id","read_count","best_sseqid","pident","qcovs","best_hit_title"])
    with counts_path.open() as f:
      for cid, n in csv.reader(f, delimiter='\\t'):
        hit = best.get(cid, {"sseqid":"NA","pident":"NA","qcovs":"NA","stitle":"NA"})
        w.writerow([cid, n, hit["sseqid"], hit["pident"], hit["qcovs"], hit["stitle"]])
  PY
  """
}

process AGGREGATE_RESULTS {
  publishDir "${params.out_dir}/summary", mode: 'copy'
  container "python:3.11-slim"

  input:
  path(tables)

  output:
  tuple path("all_samples.long.tsv"), path("abundance_matrix.tsv")

  """
  # 'tables' may expand to multiple file paths
  printf "%s\\n" ${tables} > tables.list

  python - << 'PY'
  import csv
  from pathlib import Path
  from collections import defaultdict

  tables = [Path(l.strip()) for l in open("tables.list") if l.strip()]
  long_rows = []
  taxa = set()
  samples = set()

  for t in tables:
    sample = t.name.split(".taxonomy.tsv")[0]
    samples.add(sample)
    with t.open() as f:
      r = csv.DictReader(f, delimiter='\\t')
      for row in r:
        tax = row["best_hit_title"] if row["best_hit_title"] != "NA" else "NA"
        taxa.add(tax)
        long_rows.append([sample, row["cluster_id"], int(row["read_count"]), tax, row["pident"], row["qcovs"], row["best_hit_title"]])

  with open("all_samples.long.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter='\\t')
    w.writerow(["sample","cluster_id","read_count","taxon","pident","qcovs","title"])
    w.writerows(long_rows)

  matrix = defaultdict(lambda: defaultdict(int))
  for s, _, n, tax, *_ in long_rows:
    matrix[s][tax] += n

  taxa = sorted(taxa)
  samples = sorted(samples)

  with open("abundance_matrix.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter='\\t')
    w.writerow(["sample"] + taxa)
    for s in samples:
      w.writerow([s] + [matrix[s].get(t, 0) for t in taxa])
  PY
  """
}

process TSV_TO_CSV {
  publishDir "${params.out_dir}/summary", mode: 'copy'
  container "python:3.11-slim"

  input:
  path tsv

  output:
  path "${tsv.baseName}.csv"

  """
  python - << 'PY'
  import csv
  from pathlib import Path

  tsv = Path("${tsv}")
  out = tsv.with_suffix(".csv")

  with tsv.open("r", encoding="utf-8", newline="") as fin, out.open("w", encoding="utf-8", newline="") as fout:
    r = csv.reader(fin, delimiter="\\t")
    w = csv.writer(fout)
    for row in r:
      w.writerow(row)
  PY
  """
}

process REPORT_HTML {
  publishDir "${params.out_dir}", mode: 'copy'
  container "python:3.11-slim"

  input:
  tuple path(long_csv), path(mat_csv)

  output:
  path("wf-fusarium-tef1-report.html")

  """
  python - << 'PY'
  from pathlib import Path

  long_name = Path("${long_csv}").name
  mat_name  = Path("${mat_csv}").name

  html = f'''<!doctype html>
  <html><head><meta charset="utf-8"/>
  <title>Fusarium TEF1 metabarcoding report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; }}
    .box {{ border: 1px solid #ddd; border-radius: 10px; padding: 14px; margin: 12px 0; }}
  </style>
  </head><body>
  <h1>Fusarium TEF1 metabarcoding report</h1>

  <div class="box">
    <h2>Key outputs</h2>
    <ul>
      <li><a href="summary/{long_name}">{long_name}</a></li>
      <li><a href="summary/{mat_name}">{mat_name}</a></li>
    </ul>
  </div>

  </body></html>
  '''

  Path("wf-fusarium-tef1-report.html").write_text(html, encoding="utf-8")
  PY
  """
}
