nextflow.enable.dsl=2

// =============================================================================
// BRB-seq pipeline: FastQC + multi-genome STARsolo + h5ad splitting by
// project + per-plate and per-project HTML QC reports.
//
// Alithea MERCURIUS BRB-seq kit PN 10813, V5A barcode set:
//   Read 1 = 14 bp cell barcode + 14 bp UMI (28 bp total), Read 2 = cDNA
//
// INPUT FILES
// -----------
// --input (required): samplesheet.csv, one row per pooled library/plate.
// Note: genome is NOT here -- a plate can mix genomes across wells, so
// genome lives in the plate map instead (see below).
//     library,fastq_1,fastq_2,plate_map
//     LibraryA,LibraryA_R1.fastq.gz,LibraryA_R2.fastq.gz,platemap_LibraryA.csv
//
// --genomes (required): genomes.csv, maps genome keys to STAR indices:
//     genome,star_index
//     mouse,/path/to/mouse/star_index
//     human,/path/to/human/star_index
//
// Each library's plate_map CSV maps every USED well/barcode to a project,
// sample name, AND genome:
//     barcode,project,sample_name,genome
//     TACGTTATTCCGAA,ProjectX,Sample01,mouse
//     AACAGGATAACTCC,ProjectX,Sample02,human
//     ACTCAGGCACCTCC,ProjectY,Sample03,mouse
//   A whole plate on one genome/project is just every row sharing the same
//   value -- this format doesn't force per-well granularity, it just allows
//   for it.
//
// MIXED-GENOME STRATEGY
// ----------------------
// STAR aligns a whole fastq file against one genome index per run -- it has
// no concept of "use this index for these barcodes, that index for those."
// So for each plate, we determine every distinct genome present (from its
// plate map) and run STARsolo once per genome, aligning the WHOLE plate's
// reads against each one (wasteful of compute if a plate mixes genomes --
// every genome's run processes 100% of the plate's reads -- but avoids a
// fragile pre-alignment demultiplexing-by-raw-barcode step). mtx_to_h5ad.py
// then keeps ONLY the wells actually assigned to that run's genome in the
// plate map and drops everything else; merge_plate_genomes.py recombines
// the (disjoint) per-genome fragments back into one correct whole-plate
// h5ad before anything else touches it.
//
// KNOWN LIMITATION: if a single project's samples are split across plates
// aligned to DIFFERENT genomes, the merged per-project h5ad will contain
// disjoint blocks of genes from each genome (outer-joined, zero-padded),
// not anything biologically reconciled across species. See
// bin/merge_project_h5ad.py for details.
// =============================================================================

params.input      = null   // samplesheet.csv (library,fastq_1,fastq_2,plate_map)
params.genomes    = null   // genomes.csv (genome,star_index)
params.whitelist  = null   // barcode whitelist, one sequence per line, no header
params.outdir     = './results'

process FASTQC {
    tag "$library"
    publishDir { "${params.outdir}/plates/${library}/fastqc" }, mode: 'copy'
    cpus 2

    input:
    tuple val(library), path(r1), path(r2), path(plate_map)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc --threads ${task.cpus} ${r1} ${r2}
    """
}

process SPLIT_FASTQ_BY_SAMPLE {
    tag "$library"
    publishDir { "${params.outdir}" }, mode: 'copy'
    cpus 8

    input:
    tuple val(library), path(r1), path(r2), path(plate_map)

    output:
    path "projects/**"

    script:
    """
    build_barcode_fasta.py --platemap ${plate_map} --output barcodes.fasta

    # -e 1 --no-indels: exactly 1 mismatch allowed, no indels -- matches
    # STARsolo's --soloCBmatchWLtype 1MM policy exactly (Hamming distance
    # only). --action=none: route reads without trimming/modifying them, so
    # output fastqs keep the full original R1 (barcode+UMI intact) and R2.
    # Reads that don't match any barcode land in unknown_R1/R2.fastq.gz.
    cutadapt \\
        -e 1 --no-indels \\
        -g ^file:barcodes.fasta \\
        --action=none \\
        -j ${task.cpus} \\
        -o "{name}_R1.fastq.gz" -p "{name}_R2.fastq.gz" \\
        ${r1} ${r2}

    route_by_project.py --platemap ${plate_map} --library ${library} --indir . --outdir projects
    """
}

process STARSOLO {
    tag "${library}/${genome}"
    publishDir { "${params.outdir}/plates/${library}/${genome}" }, mode: 'copy'
    cpus 32
    memory '120 GB'

    input:
    tuple val(library), path(r1), path(r2), path(plate_map), val(genome)
    path whitelist
    path genomes_csv

    output:
    tuple val(library), val(genome), path("${library}.${genome}.Solo.out"), path(plate_map), emit: solo_out
    tuple val(library), val(genome), path("${library}.${genome}.Log.final.out"),             emit: log_final
    tuple val(library), val(genome), path("${library}.${genome}.provenance.txt"),            emit: provenance

    script:
    """
    STAR_INDEX=\$(awk -F, -v g="${genome}" 'NR>1 && \$1==g {print \$2}' ${genomes_csv})
    SPECIES=\$(awk -F, -v g="${genome}" 'NR>1 && \$1==g {print \$4}' ${genomes_csv})
    ANNOTATION_SOURCE=\$(awk -F, -v g="${genome}" 'NR>1 && \$1==g {print \$5}' ${genomes_csv})
    if [ -z "\$STAR_INDEX" ]; then
        echo "ERROR: genome '${genome}' not found in ${genomes_csv}" >&2
        exit 1
    fi

    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir "\$STAR_INDEX" \\
        --readFilesIn ${r2} ${r1} \\
        --readFilesCommand zcat \\
        --soloType CB_UMI_Simple \\
        --soloCBstart 1 --soloCBlen 14 \\
        --soloUMIstart 15 --soloUMIlen 14 \\
        --soloCBwhitelist ${whitelist} \\
        --soloCBmatchWLtype 1MM \\
        --clipAdapterType CellRanger4 \\
        --soloUMIdedup 1MM_Directional \\
        --soloStrand Forward \\
        --soloCellFilter TopCells 96 \\
        --soloFeatures Gene \\
        --soloCellReadStats Standard \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        --outFileNamePrefix ${library}.${genome}.

    # STAR doesn't gzip Solo.out matrix/barcode/feature files by default,
    # but scanpy's read_10x_mtx() (used in mtx_to_h5ad.py) expects the
    # gzipped 10x-style naming convention -- mirrors what nf-core/scrnaseq's
    # own STAR_ALIGN module does for the same reason.
    find ${library}.${genome}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;

    # Provenance record: genome/annotation metadata + the EXACT command that
    # ran, captured from Nextflow's own .command.sh (not retyped/reconstructed
    # -- zero drift risk between what's documented and what actually executed).
    {
        echo "library: ${library}"
        echo "genome key: ${genome}"
        echo "species: \$SPECIES"
        echo "annotation source: \$ANNOTATION_SOURCE"
        echo "STAR index: \$STAR_INDEX"
        echo "STAR version: \$(STAR --version)"
        echo ""
        echo "Exact command executed (captured from Nextflow's own task script, not retyped):"
        echo "----------------------------------------------------------------------"
        cat .command.sh
    } > ${library}.${genome}.provenance.txt
    """
}

process MTX_TO_H5AD {
    tag "${library}/${genome}"
    publishDir { "${params.outdir}/plates/${library}/${genome}" }, mode: 'copy'
    cpus 2

    input:
    tuple val(library), val(genome), path(solo_out), path(plate_map)

    output:
    tuple val(library), path("${library}.${genome}.h5ad"),               emit: h5ad
    tuple val(library), path("${library}.${genome}.well_totals.csv"),    emit: well_totals

    script:
    """
    mtx_to_h5ad.py \\
        --solo-dir ${solo_out}/Gene/raw \\
        --cell-reads-stats ${solo_out}/Gene/CellReads.stats \\
        --platemap ${plate_map} \\
        --library ${library} \\
        --genome ${genome} \\
        --output ${library}.${genome}.h5ad \\
        --well-totals-output ${library}.${genome}.well_totals.csv
    """
}

process MERGE_PLATE_GENOMES {
    tag "$library"
    publishDir { "${params.outdir}/plates/${library}" }, mode: 'copy'
    cpus 2

    input:
    tuple val(library), path(h5ads)

    output:
    tuple val(library), path("${library}.h5ad"), emit: h5ad

    script:
    """
    merge_plate_genomes.py --inputs ${h5ads} --library ${library} --output ${library}.h5ad
    """
}

process SPLIT_H5AD_BY_PROJECT {
    tag "$library"
    cpus 2

    input:
    tuple val(library), path(h5ad)

    output:
    path "split/*.h5ad", emit: fragments

    script:
    """
    split_h5ad_by_project.py --input ${h5ad} --library ${library} --outdir split
    """
}

process PLATE_REPORT {
    tag "$library"
    publishDir { "${params.outdir}/plates/${library}" }, mode: 'copy'
    cpus 2

    input:
    tuple val(library), path(h5ad), path(log_finals), path(well_totals)

    output:
    path "${library}_report.html"

    script:
    """
    plate_report.py --h5ad ${h5ad} --log-final ${log_finals} --well-totals ${well_totals} --library ${library} --output ${library}_report.html
    """
}

process MERGE_PROJECT_H5AD {
    tag "$project"
    publishDir { "${params.outdir}/projects/${project}" }, mode: 'copy'
    cpus 2

    input:
    tuple val(project), path(fragments)

    output:
    tuple val(project), path("${project}.h5ad"), emit: h5ad

    script:
    """
    merge_project_h5ad.py --inputs ${fragments} --project "${project}" --output "${project}.h5ad"
    """
}

process PROJECT_REPORT {
    tag "$project"
    publishDir { "${params.outdir}/projects/${project}" }, mode: 'copy'
    cpus 2

    input:
    tuple val(project), path(h5ad)

    output:
    path "${project}_report.html"

    script:
    """
    project_report.py --h5ad ${h5ad} --project "${project}" --output "${project}_report.html"
    """
}

process EXPORT_GENE_TABLE {
    tag "$project"
    publishDir { "${params.outdir}/projects/${project}" }, mode: 'copy'
    cpus 2

    input:
    tuple val(project), path(h5ad)

    output:
    tuple path("${project}_counts.tsv"), path("${project}_metadata.tsv"), path("${project}_gene_id_map.tsv")

    script:
    """
    export_gene_table.py --h5ad ${h5ad} --project "${project}" \\
        --counts-output "${project}_counts.tsv" \\
        --metadata-output "${project}_metadata.tsv" \\
        --gene-id-map-output "${project}_gene_id_map.tsv"
    """
}

process BUILD_GENE_ANNOTATION {
    tag "$genome"
    publishDir { "${params.outdir}/reference" }, mode: 'copy'
    cpus 1

    input:
    tuple val(genome), path(gtf)

    output:
    tuple val(genome), path("${genome}_gene_annotation.tsv")

    script:
    """
    gene_annotation_table.py --gtf ${gtf} --output ${genome}_gene_annotation.tsv
    """
}

process EXPORT_PUBLICATION_TABLE {
    tag "$project"
    publishDir { "${params.outdir}/projects/${project}" }, mode: 'copy'
    cpus 2

    input:
    tuple val(project), path(h5ad), path(gene_annotations)

    output:
    path "${project}_publication_table.tsv"

    script:
    """
    export_publication_table.py --h5ad ${h5ad} \\
        --gene-annotation ${gene_annotations} \\
        --project "${project}" --output "${project}_publication_table.tsv"
    """
}

process RENDER_METHODS {
    tag "$project"
    publishDir { "${params.outdir}/projects/${project}" }, mode: 'copy'
    cpus 1

    input:
    tuple val(project), path(h5ad), path(provenance)

    output:
    path "${project}_methods.txt"

    script:
    """
    render_methods.py --h5ad ${h5ad} --project "${project}" \\
        --provenance ${provenance} --output "${project}_methods.txt"
    """
}

workflow {
    ch_genomes = file(params.genomes)
    ch_whitelist = file(params.whitelist)

    // one BUILD_GENE_ANNOTATION run per distinct genome in genomes.csv
    // (not per plate/library -- a genome's annotation is fixed regardless
    // of how many plates use it)
    ch_genome_gtfs = Channel.fromPath(params.genomes)
        .splitCsv(header: true)
        .map { row -> tuple(row.genome, file(row.gtf)) }
    BUILD_GENE_ANNOTATION(ch_genome_gtfs)
    // combined into one list usable by every project's export, since a
    // project's samples could in principle span more than one genome
    ch_all_gene_annotations = BUILD_GENE_ANNOTATION.out.map { genome, tsv -> tsv }.collect()

    ch_libraries = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.library, file(row.fastq_1), file(row.fastq_2), file(row.plate_map)) }

    FASTQC(ch_libraries)
    SPLIT_FASTQ_BY_SAMPLE(ch_libraries)

    // Determine the distinct genomes present on each plate from its plate
    // map, and fan out one whole-plate STARsolo run per (library, genome)
    // pair. Plain Groovy CSV parsing here rather than a Nextflow-specific
    // file.splitCsv() call, since that API's behavior on a bare Path
    // (outside a channel) isn't something to rely on without checking --
    // readLines()+split(',') is guaranteed to work.
    ch_libraries_by_genome = ch_libraries
        .flatMap { library, r1, r2, plate_map ->
            def lines = plate_map.readLines()
            def header = lines[0].split(',')
            def genomeIdx = header.findIndexOf { it.trim() == 'genome' }
            if (genomeIdx == -1) {
                error "plate_map ${plate_map} has no 'genome' column"
            }
            def genomes = lines.drop(1)
                .findAll { it.trim() }
                .collect { it.split(',')[genomeIdx].trim() }
                .unique()
            genomes.collect { genome -> tuple(library, r1, r2, plate_map, genome) }
        }

    STARSOLO(ch_libraries_by_genome, ch_whitelist, ch_genomes)

    // every (library, genome) provenance record, collected once and handed
    // to every project -- render_methods.py filters down to what's actually
    // relevant to each project's own h5ad
    ch_all_provenance = STARSOLO.out.provenance.map { library, genome, txt -> txt }.collect()

    MTX_TO_H5AD(STARSOLO.out.solo_out)

    // recombine each plate's per-genome fragments (disjoint wells) into one
    // correct whole-plate h5ad
    ch_plate_fragments = MTX_TO_H5AD.out.h5ad.groupTuple()
    MERGE_PLATE_GENOMES(ch_plate_fragments)

    // per-plate report needs the combined h5ad + every STARsolo log for
    // that plate (one per genome) + every per-genome well_totals CSV
    // (needed to cross-check assigned vs. off-target genome UMI counts)
    ch_plate_logs = STARSOLO.out.log_final
        .map { library, genome, log -> tuple(library, log) }
        .groupTuple()
    ch_plate_totals = MTX_TO_H5AD.out.well_totals.groupTuple()
    ch_plate_report_in = MERGE_PLATE_GENOMES.out.h5ad
        .join(ch_plate_logs)
        .join(ch_plate_totals)
    PLATE_REPORT(ch_plate_report_in)

    // split each plate's (already genome-correct) h5ad into per-project
    // fragments, then regroup ACROSS plates by project alone so a project
    // spanning multiple plates gets one merged file/report
    SPLIT_H5AD_BY_PROJECT(MERGE_PLATE_GENOMES.out.h5ad)

    ch_by_project = SPLIT_H5AD_BY_PROJECT.out.fragments
        .flatten()
        .map { f -> tuple(f.baseName.replaceFirst(/^.*__/, ''), f) }  // filename: library__project.h5ad (greedy match handles underscores in library names)
        .groupTuple()

    MERGE_PROJECT_H5AD(ch_by_project)
    PROJECT_REPORT(MERGE_PROJECT_H5AD.out.h5ad)
    EXPORT_GENE_TABLE(MERGE_PROJECT_H5AD.out.h5ad)
    EXPORT_PUBLICATION_TABLE(MERGE_PROJECT_H5AD.out.h5ad.combine(ch_all_gene_annotations))
    RENDER_METHODS(MERGE_PROJECT_H5AD.out.h5ad.combine(ch_all_provenance))
}
