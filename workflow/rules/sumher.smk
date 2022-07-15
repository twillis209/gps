from scipy.stats import chi2
import re

rule deduplicate_variants:
    input:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
    output:
        "resources/1000g/euro/qc/nodup/{chr}.bed",
        "resources/1000g/euro/qc/nodup/{chr}.bim",
        "resources/1000g/euro/qc/nodup/{chr}.fam"
    params:
        input_stem = "resources/1000g/euro/qc/{chr}_qc",
        output_stem = "resources/1000g/euro/qc/nodup/{chr}"
    threads: 1
    resources:
        mem_mb=get_mem_mb
    group: "sumher"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.input_stem} --rm-dup 'force-first' --make-bed --silent --out {params.output_stem}"

rule retain_only_snp_variants:
    input:
        "resources/1000g/euro/qc/nodup/{chr}.bed",
        "resources/1000g/euro/qc/nodup/{chr}.bim",
        "resources/1000g/euro/qc/nodup/{chr}.fam"
    output:
        "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}.fam"
    params:
        input_stem = "resources/1000g/euro/qc/nodup/{chr}",
        output_stem = "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}"
    threads: 1
    resources:
        mem_mb=get_mem_mb
    group: "sumher"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.input_stem} --snps-only --make-bed --silent --out {params.output_stem}"

rule subset_snp_variants:
    input:
        "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}.fam",
        range_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/matching_ids/{chr}.txt"
    output:
        temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.bed"),
        temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.bim"),
        temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.fam"),
    params:
        input_stem = "resources/1000g/euro/qc/nodup/snps_only/1000g/{chr}",
        output_stem = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}"
    threads: 1
    resources:
        mem_mb=get_mem_mb
    group: "sumher"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.input_stem} --snps-only --make-bed --extract {input.range_file} --silent --out {params.output_stem}"

rule thin_predictors:
    input:
        "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.bed",
        "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.bim",
        "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.fam",
    output:
        thin_file = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/thin.in"),
        weights_file = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/weights.thin")
    log:
        log_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/thin.log"
    params:
        input_stem = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}",
        output_stem = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/thin"
    group: "sumher"
    shell:
        """
        $ldakRoot/ldak --thin {params.output_stem} --bfile {params.input_stem} --window-prune .98 --window-kb 100 > {log.log_file};
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

rule calculate_ldak_thin_taggings_for_chromosome:
    input:
        "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.bed",
        "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.bim",
        "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}.fam",
        weights_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/weights.thin"
    output:
        tagging_file = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/{chr}.tagging"),
    log:
        log_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/{chr}.tagging.log"
    params:
        input_stem = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}",
        output_stem = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/{chr}/{chr}"
    group: "sumher"
    shell:
        "$ldakRoot/ldak --calc-tagging {params.output_stem} --bfile {params.input_stem} --weights {input.weights_file} --chr {wildcards.chr} --window-kb 1000 --power -.25 > {log.log_file}"

rule join_ldak_thin_taggings:
    input:
        [f"results/merged_gwas/{{trait_A}}_{{trait_B}}/{{trait_A}}_{{trait_B}}_{{snp_set}}/ldak/chr{x}/chr{x}.tagging" for x in range(1, 23)]
    output:
        wg_tagging_file = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/whole_genome.tagging"),
        chrom_taggings_file = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/taggings.txt")
    log:
        log_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/whole_genome.tagging.log"
    params:
        output_stem = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/whole_genome"
    group: "sumher"
    shell:
        """
        for x in {input}; do
            echo $x >> {output.chrom_taggings_file}
        done;

        $ldakRoot/ldak --join-tagging {params.output_stem} --taglist {output.chrom_taggings_file} > {log.log_file}
        """

rule process_sum_stats:
    input:
        gwas_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_A}_{trait_B}_{snp_set}.tsv.gz",
        metadata_file = "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        gwas_file_A = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_A}.assoc"),
        gwas_file_B = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_B}.assoc")
    params:
        trait_A = lambda wildcards: wildcards.trait_A,
        trait_B = lambda wildcards: wildcards.trait_B,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        beta_a_col = 'BETA.A',
        beta_b_col = 'BETA.B',
        se_a_col = 'SE.A',
        se_b_col = 'SE.B'
    threads: 8
    resources:
        time = 20
    group: "sumher"
    script:
        "../scripts/process_sum_stats.R"

rule estimate_rg_with_ldak_thin:
    input:
        wg_tagging_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/ldak/whole_genome.tagging",
        gwas_file_A = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_A}.assoc",
        gwas_file_B = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_B}.assoc"
    output:
        progress_file = "results/ldak/ldak-thin/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}.progress",
        cors_file = "results/ldak/ldak-thin/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}.cors",
        cors_full_file = "results/ldak/ldak-thin/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}.cors.full",
        labels_file = "results/ldak/ldak-thin/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}.labels",
        overlap_file = "results/ldak/ldak-thin/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}.overlap",
    log:
        log_file = "results/ldak/ldak-thin/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}.log"
    params:
        output_stem = "results/ldak/ldak-thin/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}"
    resources:
        time = 5
    group: "sumher"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.gwas_file_A} --summary2 {input.gwas_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """
#rule compile_overlap_files_for_top_imds:
#    input:
#    output:
#        "results/ldak/ldak-thin/top_imds.overlap"
#    shell:
#        """
#        echo -e "trait_A trait_B snp_set tag_predictors assoc_predictors overlap" >{output}
#        for x in {input}; do echo -n -e "$(basename $x) " | tr '_' ' '; cut -d' ' -f2 $x | tr '\n' ' '; echo; done | sed 's/\.overlap//' >>{output}
#        """
