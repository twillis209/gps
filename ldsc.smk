rule download_ld_scores:
    output:
        "resources/ldsc/eur_w_ld_chr.tar.bz2"
    shell:
        "wget -O {output} https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2"

rule extract_ld_scores:
    input:
        "resources/ldsc/eur_w_ld_chr.tar.bz2"
    output:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        "resources/ldsc/eur_w_ld_chr/w_hm3.snplist",
        "resources/ldsc/eur_w_ld_chr/README"
    params:
        output_root = "resources/ldsc"
    shell:
        "tar -xjf {input} -C {params.output_root}"

rule munge_whole_genome_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_sum_stats.tsv.gz"
    output:
        # NB: Awkward output filename due to the way the LDSC script works
        temp("results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases,\d+}_{ncontrols,\d+}/{effect_blocks}_{tag,\w}.tsv.sumstats.gz")
    params:
         output_filename = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv",
         signed_sumstats_col = lambda wildcards: tag_signed_sumstats_dict[wildcards.tag],
         pvalue_col = lambda wildcards: tag_pvalue_dict[wildcards.tag]
    resources:
        disk_mb = 2048
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-con-col ncontrols --N-cas-col ncases --snp id --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a0 --a2 a1 --frq EUR --a1-inc;
        """
