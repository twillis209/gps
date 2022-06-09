rule join_pair_gwas:
    input:
        A = "resources/gwas/{trait_A}.tsv.gz",
        B = "resources/gwas/{trait_B}.tsv.gz"
    output:
        AB = temp("resources/gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_A}_{trait_B}_{snp_set}.tsv.gz")
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
        principal_trait_gwas_file = "resources/gwas/{prin_trait}.tsv.gz",
        auxiliary_trait_gwas_files = lambda wildcards: [f"resources/gwas/{x}.tsv.gz" for x in wildcards.aux_traits.split('+')],
        pruned_auxiliary_trait_gwas_files = lambda wildcards: [f"resources/gwas/{wildcards.prin_trait}_{x}/{wildcards.prin_trait}_{x}_{wildcards.snp_set}/pruned_{wildcards.prin_trait}_{x}_{wildcards.snp_set}.tsv.gz" for x in wildcards.aux_traits.split('+')]
    output:
        temp("resources/gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_with_prune_{snp_set}.tsv.gz")
    threads: 4
    params:
        mhc = lambda wildcards: '-sans_mhc' if wildcards.snp_set == 'sans-mhc' else '',
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

rule make_plink_ranges:
     input:
      ("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      gwas_file = "resources/gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/{trait_A}_{trait_B}_{snp_set}.tsv.gz"
     output:
      ("resources/gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/matching_ids/chr%d.txt" % x for x in range(1,23))
     params:
        input_dir = "resources/1000g/euro/qc",
        output_dir = "resources/gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{snp_set}/matching_ids"
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
      range_file = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/matching_ids/{chr}.txt"
     output:
      temp("resources/gwas/pid_{imd,[^\+]+}/pid_{imd}_{snp_set}/plink/{chr}.bed"),
      temp("resources/gwas/pid_{imd,[^\+]+}/pid_{imd}_{snp_set}/plink/{chr}.bim"),
      temp("resources/gwas/pid_{imd,[^\+]+}/pid_{imd}_{snp_set}/plink/{chr}.fam")
     params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_pruned_ranges:
     input:
      "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}.bed",
      "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}.bim",
      "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}.fam"
     output:
      "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/{chr}.prune.in",
      "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/{chr}.prune.out"
     params:
       bfile = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/plink/{chr}",
       prune_out = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/{chr}"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise 1000kb 50 0.2 --out {params.prune_out}"

rule cat_pruned_ranges:
     input:
      ("resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/chr%d.prune.in" % x for x in range(1,23))
     output:
      "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/all.prune.in"
     threads: 1
     group: "gps"
     shell:
      "for x in {input}; do cat $x >>{output}; done"

rule prune_gwas:
     input:
      prune_file = ancient("resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/prune/all.prune.in"),
      gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pid_{imd}_{snp_set}.tsv.gz"
     output:
      "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv.gz"
     threads: 1
     group: "gps"
     shell:
      "Rscript workflow/scripts/prune_gwas.R -i {input.gwas_file} -p {input.prune_file} -o {output} -nt {threads}"

rule unzip_pruned_joined_gwas:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv.gz"
    output:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv"
    group: "gps"
    shell:
        "gunzip -c {input} >{output}"
