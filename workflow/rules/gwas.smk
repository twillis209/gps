rule download_gwas:
    output:
        temp("resources/gwas/{trait}.tsv.gz")
    params:
        url = lambda w: metadata_daf.loc[metadata_daf['abbrv'] == w.trait, 'URL'].values[0]
    shell:
        """
        wget -O {output} {params.url}
        """

        # TODO rewrite pipeline to handle temporary dir, work without cd etc.
        # TODO pipeline is currently handling rm of a lot of stuff
        # TODO can only be run serially
rule process_gwas:
    input:
        ancient("resources/gwas/{trait}.tsv.gz")
    output:
        temp_input_cp = temp("workflow/scripts/GWAS_tools/01-Pipeline/{trait}.tsv.gz"),
        processed_file = "results/processed_gwas/{trait,[^\_]+}.tsv.gz"
    params:
        gwas_tools_dir = "workflow/scripts/GWAS_tools/01-Pipeline",
        pipeline_output_file = lambda w: f"workflow/scripts/GWAS_tools/01-Pipeline/{w.trait}-hg38.tsv.gz",
        is_ukbb = lambda w: "true" if w.trait in ukbb_traits else "false",
        temp_input_cp_decompressed_name = "{trait}.tsv",
        temp_input_cp_name = "{trait}.tsv.gz"
    shell:
        """
        # Need to remove the dash as this has special meaning in the pipeline

        cp {input} {output.temp_input_cp}
        cd {params.gwas_tools_dir}

        if [ "{params.is_ukbb}" = "true" ]; then
            zcat {params.temp_input_cp_name} | tr ':' '\t' >{params.temp_input_cp_decompressed_name};
            sed -i 's/variant/chrom\tbp\tref\talt/' {params.temp_input_cp_decompressed_name};
            gzip -f {params.temp_input_cp_decompressed_name}
        fi

        ./pipeline_v5.3.2_beta.sh -f {wildcards.trait}.tsv.gz -b {wildcards.trait}
        cd ../../../..
        mv {params.pipeline_output_file} {output.processed_file}
        """

rule recalculate_p_values:
    input:
        ancient("results/processed_gwas/{trait}.tsv.gz")
    output:
        "results/processed_gwas/{trait,[^\_]+}_recalculated_p.tsv.gz"
    params:
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P'
    threads: 4
    script: "../scripts/recalculate_p_values.R"

rule join_pair_gwas:
    input:
        A = "results/processed_gwas/{trait_A}_recalculated_p.tsv.gz",
        B = "results/processed_gwas/{trait_B}_recalculated_p.tsv.gz"
    output:
        AB = temp("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_A}_{trait_B}_{snp_set}.tsv.gz")
    threads: 4
    params:
        mhc = lambda wildcards: False if wildcards.snp_set == 'sans-mhc' else True,
        join = 'inner',
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P',
        beta_col = 'BETA',
        se_col = 'SE',
        recalculate_p = True
    group: "gps"
    script:
        "../scripts/join_pair_gwas_stats.R"

rule join_multiple_gwas:
    input:
        principal_trait_gwas_file = "results/processed_gwas/{prin_trait}.tsv.gz",
        auxiliary_trait_gwas_files = lambda wildcards: [f"results/processed_gwas/{x}.tsv.gz" for x in wildcards.aux_traits.split('+')],
        pruned_auxiliary_trait_gwas_files = lambda wildcards: [f"results/processed_gwas/{wildcards.prin_trait}_{x}/{wildcards.prin_trait}_{x}_{wildcards.snp_set}/pruned_{wildcards.prin_trait}_{x}_{wildcards.snp_set}.tsv.gz" for x in wildcards.aux_traits.split('+')]
    output:
        temp("results/merged_gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_with_prune_{snp_set}.tsv.gz")
    threads: 4
    params:
        mhc = lambda wildcards: False if wildcards.snp_set == 'sans-mhc' else True,
        bp_col = 'BP38',
        chr_col = 'CHR38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P'
    resources:
        time = 10
    group: "gps"
    script:
        "../scripts/join_multiple_gwas_stats.R"

rule add_ld_block_column_to_merged_gwas:
    input:
        gwas_file = "results/merged_gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_with_prune_{snp_set}.tsv.gz",
        block_file = "resources/ldetect/blocks.tsv"
    output:
        "results/merged_gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_with_prune_and_blocks_{snp_set}.tsv.gz"
    params:
        bp_col = 'BP38',
        chr_col = 'CHR38',
        ref_col = 'REF',
        alt_col = 'ALT'
    threads: 4
    resources:
        time = 180
    group: "gps"
    script:
        "../scripts/add_ld_block_column_to_merged_gwas.R"

rule make_plink_ranges:
     input:
      ("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      gwas_file = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_A}_{trait_B}_{snp_set}.tsv.gz"
     output:
      ("results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/matching_ids/chr%d.txt" % x for x in range(1,23))
     params:
        input_dir = "resources/1000g/euro/qc",
        output_dir = "results/merged_gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/matching_ids"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "Rscript workflow/scripts/make_plink_ranges.R -i {input.gwas_file} -b {params.input_dir} -r chr%d_qc.bim -o {params.output_dir} -nt {threads}"

rule subset_reference:
     input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam",
      range_file = "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/matching_ids/{chr}.txt"
     output:
      temp("results/merged_gwas/pid_{imd,[^\+]+}/pid_{imd}_{snp_set}/plink/{chr}.bed"),
      temp("results/merged_gwas/pid_{imd,[^\+]+}/pid_{imd}_{snp_set}/plink/{chr}.bim"),
      temp("results/merged_gwas/pid_{imd,[^\+]+}/pid_{imd}_{snp_set}/plink/{chr}.fam")
     params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_pruned_ranges:
     input:
      "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}.bed",
      "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}.bim",
      "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}.fam"
     output:
      "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/{chr}.prune.in",
      "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/{chr}.prune.out"
     params:
       bfile = "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}",
       prune_out = "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/{chr}"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise 1000kb 50 0.2 --out {params.prune_out}"

rule cat_pruned_ranges:
     input:
      ("results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/chr%d.prune.in" % x for x in range(1,23))
     output:
      "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/all.prune.in"
     threads: 1
     group: "gps"
     shell:
      "for x in {input}; do cat $x >>{output}; done"

rule prune_gwas:
     input:
      prune_file = ancient("results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/all.prune.in"),
      gwas_file = "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/pid_{imd}_{snp_set}.tsv.gz"
     output:
      "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv.gz"
     threads: 1
     group: "gps"
     shell:
      "Rscript workflow/scripts/prune_gwas.R -i {input.gwas_file} -p {input.prune_file} -o {output} -nt {threads}"

rule unzip_pruned_joined_gwas:
    input:
        "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv.gz"
    output:
        "results/merged_gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv"
    group: "gps"
    shell:
        "gunzip -c {input} >{output}"
