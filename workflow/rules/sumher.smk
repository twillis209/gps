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
    threads: 4
    resources:
        time = 10
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
#        [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.overlap" for trait_pair in top_imd_pairs if os.path.exists(f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.overlap")],
#        [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_sans-mhc.overlap" for trait_pair in top_imd_pairs if os.path.exists(f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.overlap")]
#    output:
#        "results/ldak/ldak-thin/top_imds.overlap"
#    shell:
#        """
#        echo -e "trait_A trait_B snp_set tag_predictors assoc_predictors overlap" >{output}
#        for x in {input}; do echo -n -e "$(basename $x) " | tr '_' ' '; cut -d' ' -f2 $x | tr '\n' ' '; echo; done | sed 's/\.overlap//' >>{output}
#        """
#
#rule estimate_rg_with_ldak_thin_for_ukbb:
#    input:
#        wg_tagging_file = "results/ldak/ldak-thin/ukbb/whole_genome.tagging",
#        sum_stats_file_A = "resources/ukbb_sum_stats/{trait_A}.assoc",
#        sum_stats_file_B = "resources/ukbb_sum_stats/{trait_B}.assoc"
#    output:
#        progress_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.progress"
##        cors_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors",
##        cors_full_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors.full",
##        labels_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors.labels",
##        overlap_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors.overlap",
#    log:
#        log_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.log"
#    params:
#        output_stem = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}"
#    resources:
#        time = 5
#    group: "ukbb_sumher"
#    shell:
#        """
#        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
#        """
#
#rule ukbb_sumher:
#    input:
#       sumher_files = [f"results/ldak/ldak-thin/ukbb/rg/{trait_pair}.cors.full" for trait_pair in ukbb_trait_pairs],
#       metadata_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
#    output:
#        "results/ldak/ldak-thin/ukbb/rg/compiled_ukbb_sumher_results.tsv"
#    run:
#        meta_daf = pd.read_csv(input.metadata_file, sep = '\t')
#        with open(output[0], 'w') as outfile:
#            outfile.write("trait_A_code\ttrait_B_code\ttrait_A_name\ttrait_B_name\th2.A\th2.A.se\th2.B\th2.B.se\tgcov\tgcov.se\trg\trg.se\trg.z\trg.p\n")
#            for x in input.sumher_files:
#                head, tail = os.path.split(x)
#
#                trait_A, trait_B = re.match("(\w+)-(\w+).cors.full", tail).groups()
#
#                trait_A_name = meta_daf.loc[meta_daf.code == trait_A, 'long_abbrv'].values[0]
#                trait_B_name = meta_daf.loc[meta_daf.code == trait_B, 'long_abbrv'].values[0]
#
#                with open(x, 'r') as infile:
#                    line = infile.readline()
#                    line = infile.readline()
#
#                    # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
#                    _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()
#                    rg_z = float(rg)/float(rg_se)
#
#                    rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)
#                    outfile.write(f"{trait_A}\t{trait_B}\t{trait_A_name}\t{trait_B_name}\t{float(h2_A)}\t{float(h2_A_se)}\t{float(h2_B)}\t{float(h2_B_se)}\t{float(cov)}\t{float(cov_se)}\t{float(rg)}\t{float(rg_se)}\t{rg_z}\t{rg_p}\n")
