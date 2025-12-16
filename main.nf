nextflow.enable.dsl = 2

// --------------------
// Channels (input: reads are already separated by sample/barcode folders)
// --------------------
Channel
  .fromPath("${params.reads}/*/*.fastq.gz", checkIfExists: true)
  .map { f -> tuple(f.parent.name, f) }
  .groupTuple()
  .map { sample, files -> tuple(sample, files.sort()) }
  .set { ch_samples }

// --------------------
// Workflow
// --------------------
workflow {

  merged = MERGE_FASTQ(ch_samples)

  qc_raw = QC_SEQKIT(merged)

  trimmed = (params.primer_fwd?.trim() && params.primer_rev?.trim()) ? TRIM_CUTADAPT(merged) : merged

  filtered = FILTER_NANOFILT(trimmed)

  fasta = FASTQ_TO_FASTA(filtered)

  clustered = CLUSTER_VSEARCH(fasta)

  kept = FILTER_CLUSTERS(clustered)

  db = MAKEBLASTDB(params.db_fasta)

  blasted = BLASTN(kept, db)

  tax = JOIN_COUNTS_BLAST(blasted)

  agg = AGGREGATE_RESULTS(tax)

  REPORT_HTML(agg)
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
  # Concatenating gz streams is valid (keeps compression, fast, no re-gzip)
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
  tuple val(sample), path("${sample}.seqkit.stats.tsv"), emit: stats

  """
  seqkit stats -a -T ${fq} > ${sample}.seqkit.stats.tsv
  """
}

process TRIM_CUTADAPT {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/reads", mode: 'copy'

  container "quay.io/biocontainers/cutadapt:4.9--py311h1f90b4d_0"

  input:
  tuple val(sample), path(fq)

  output:
  tuple val(sample), path("${sample}.trimmed.fastq.gz")

  def discard = params.require_primers ? "--discard-untrimmed" : ""
  """
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

  container "quay.io/biocontainers/nanofilt:2.8.0--pyhdfd78af_0"

  input:
  tuple val(sample), path(fq)

  output:
  tuple val(sample), path("${sample}.filtered.fastq.gz")

  """
  zcat ${fq} \
    | NanoFilt -q ${params.min_q} -l ${params.min_len} --maxlength ${params.max_len} \
    | gzip > ${sample}.filtered.fastq.gz
  """
}

process FASTQ_TO_FASTA {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/clusters", mode: 'copy'

  container "quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0"

  input:
  tuple val(sample), path(fq)

  output:
  tuple val(sample), path("${sample}.fasta")

  """
  seqkit fq2fa ${fq} > ${sample}.fasta
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
  # Build centroid counts from UC:
  # - "S" = seed/centroid line (query label in col 9)
  # - "H" = hit assigned to a centroid (cluster number in col 2)
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
    # No clusters pass threshold -> create empty fasta so downstream steps are predictable
    : > ${sample}.centroids.kept.fasta
  fi
  """
}

process MAKEBLASTDB {
  publishDir "${params.out_dir}/refdb", mode: 'copy'

  container "quay.io/biocontainers/blast:2.16.0--pl5321h6f7f691_0"

  input:
  path(db_fasta)

  output:
  path("blastdb"), emit: dbdir

  """
  mkdir -p blastdb
  cp ${db_fasta} blastdb/db.fasta
  makeblastdb -in blastdb/db.fasta -dbtype nucl -parse_seqids -out blastdb/fusarium_tef1
  """
}

process BLASTN {
  tag "$sample"
  publishDir "${params.out_dir}/per_sample/${sample}/taxonomy", mode: 'copy'

  container "quay.io/biocontainers/blast:2.16.0--pl5321h6f7f691_0"

  input:
  tuple val(sample), path(counts), path(centroids_fa)
  path(dbdir)

  output:
  tuple val(sample),
        path(counts),
        path("${sample}.blast.tsv")

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
  path("${sample}.taxonomy.tsv")

  """
  python - << 'PY'
  import csv
  from pathlib import Path

  sample = "${sample}"
  counts_path = Path("${counts_tsv}")
  blast_path  = Path("${blast_tsv}")
  out_path    = Path(f"{sample}.taxonomy.tsv")

  # Read best hit per query (first occurrence; BLAST is typically sorted by score)
  best = {}
  if blast_path.exists() and blast_path.stat().st_size > 0:
    with blast_path.open() as f:
      for row in csv.reader(f, delimiter='\\t'):
        qseqid = row[0]
        if qseqid in best:
          continue
        best[qseqid] = {
          "sseqid": row[1],
          "pident": row[2],
          "qcovs": row[7],
          "stitle": row[8] if len(row) > 8 else ""
        }

  with out_path.open("w", newline="") as out:
    w = csv.writer(out, delimiter='\\t')
    w.writerow(["cluster_id","read_count","best_sseqid","pident","qcovs","best_hit_title"])
    with counts_path.open() as f:
      for row in csv.reader(f, delimiter='\\t'):
        cid, n = row[0], row[1]
        hit = best.get(cid, {"sseqid":"NA","pident":"NA","qcovs":"NA","stitle":"NA"})
        w.writerow([cid, n, hit["sseqid"], hit["pident"], hit["qcovs"], hit["stitle"]])
  PY
  """
}

process AGGREGATE_RESULTS {
  publishDir "${params.out_dir}/summary", mode: 'copy'

  container "python:3.11-slim"

  input:
  path(tables) from Channel.fromPath("${params.out_dir}/per_sample/*/taxonomy/*.taxonomy.tsv").collect()

  output:
  path("all_samples.long.tsv"),
  path("abundance_matrix.tsv")

  """
  python - << 'PY'
  import csv
  from pathlib import Path
  from collections import defaultdict

  tables = [Path(p) for p in "${tables}".split()]
  long_rows = []
  taxa = set()
  samples = set()

  # Long format
  for t in tables:
    sample = t.name.split(".taxonomy.tsv")[0]
    samples.add(sample)
    with t.open() as f:
      r = csv.DictReader(f, delimiter='\\t')
      for row in r:
        tax = row["best_sseqid"] if row["best_sseqid"] != "NA" else "NA"
        taxa.add(tax)
        long_rows.append({
          "sample": sample,
          "cluster_id": row["cluster_id"],
          "read_count": int(row["read_count"]),
          "taxon": tax,
          "pident": row["pident"],
          "qcovs": row["qcovs"],
          "title": row["best_hit_title"],
        })

  # Write long table
  with open("all_samples.long.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter='\\t')
    w.writerow(["sample","cluster_id","read_count","taxon","pident","qcovs","title"])
    for row in long_rows:
      w.writerow([row["sample"], row["cluster_id"], row["read_count"], row["taxon"], row["pident"], row["qcovs"], row["title"]])

  # Pivot abundance: sample x taxon
  matrix = defaultdict(lambda: defaultdict(int))
  for row in long_rows:
    matrix[row["sample"]][row["taxon"]] += row["read_count"]

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

process REPORT_HTML {
  publishDir "${params.out_dir}", mode: 'copy'

  container "python:3.11-slim"

  input:
  path("all_samples.long.tsv"),
  path("abundance_matrix.tsv")

  output:
  path("wf-fusarium-tef1-report.html")

  """
  python - << 'PY'
  from pathlib import Path

  long_tsv = Path("all_samples.long.tsv")
  mat_tsv  = Path("abundance_matrix.tsv")

  html = f\"\"\"<!doctype html>
  <html>
  <head>
    <meta charset="utf-8"/>
    <title>Fusarium TEF1 metabarcoding report</title>
    <style>
      body {{ font-family: Arial, sans-serif; margin: 24px; }}
      code, pre {{ background: #f5f5f5; padding: 2px 4px; }}
      .box {{ border: 1px solid #ddd; border-radius: 10px; padding: 14px; margin: 12px 0; }}
      a {{ text-decoration: none; }}
    </style>
  </head>
  <body>
    <h1>Fusarium TEF1 metabarcoding report</h1>

    <div class="box">
      <h2>Key outputs</h2>
      <ul>
        <li><a href="summary/all_samples.long.tsv">all_samples.long.tsv</a></li>
        <li><a href="summary/abundance_matrix.tsv">abundance_matrix.tsv</a></li>
      </ul>
    </div>

    <div class="box">
      <h2>Notes</h2>
      <ul>
        <li>Taxonomy is based on BLAST against the provided TEF1 reference FASTA.</li>
        <li>Cluster representatives are VSEARCH centroids (not a polished consensus yet).</li>
        <li>You can tighten/relax clustering with <code>cluster_id</code> and <code>min_cluster_reads</code>.</li>
      </ul>
    </div>
  </body>
  </html>
  \"\"\"

  Path("wf-fusarium-tef1-report.html").write_text(html, encoding="utf-8")
  PY
  """
}
