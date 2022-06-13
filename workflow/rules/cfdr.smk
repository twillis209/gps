rule run_pairwise_cfdr:
    input:
        gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pid_{imd}_{snp_set}.tsv.gz",
        pruned_gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv.gz"
    output:
        results_file = "results/cfdr/pid_{imd}/pid_{imd}_{snp_set}/pid_{imd}.tsv.gz"
    benchmark:
        "benchmarks/run_pairwise_cfdr/pid_{imd}_{snp_set}_benchmark.txt"
    threads: 6
    resources:
        time = 90
    params:
        chrom_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        prin_col = 'P.A',
        aux_col = 'P.B',
        p_threshold = 1e-2,
        v_col = 'v.B'
    script:
         "../scripts/pairwise_cfdr.R"

rule run_iterative_cfdr:
    input:
        gwas_file = "resources/gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_with_prune_{snp_set}.tsv.gz",
        ldetect_blocks_file = "resources/ldetect/blocks.txt"
    output:
        "results/cfdr/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}.tsv.gz"
    benchmark:
        "benchmarks/run_iterative_cfdr/{prin_trait}_{aux_traits}_{snp_set}_benchmark.txt"
    threads: 8
    resources:
        time = 600
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        prin_col = 'P',
        aux_cols = lambda wildcards: [f'P.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('+'))],
        prune_cols = lambda wildcards: [f'prune_in.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('+'))],
        v_cols = lambda wildcards: [f'v.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('+'))],
        p_threshold = 1e-2,
        use_ldetect_blocks = False
    script:
        "../scripts/iterative_cfdr.R"

rule run_iterative_cfdr_using_ldetect_blocks:
    input:
        "resources/gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_with_prune_{snp_set}.tsv.gz"
    output:
        "results/cfdr/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_ldetect.tsv.gz"
    benchmark:
        "benchmarks/run_iterative_cfdr/{prin_trait}_{aux_traits}_{snp_set}_ldetect_benchmark.txt"
    threads: 8
    resources:
        time = 600
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        prin_col = 'P',
        aux_cols = lambda wildcards: [f'P.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('+'))],
        prune_cols = lambda wildcards: [f'prune_in.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('+'))],
        v_cols = lambda wildcards: [f'v.{i+1}' for i,x in enumerate(wildcards.aux_traits.split('+'))],
        p_threshold = 1e-2,
        use_ldetect_blocks = True
    script:
        "../scripts/iterative_cfdr.R"

rule run_iterative_cfdr_for_top_traits_per_gps:
    input:
        "results/cfdr/pid_uc-delange+sys-scl+lada+t1d-cooper+psc+jia+addi+ms+igad/pid_uc-delange+sys-scl+lada+t1d-cooper+psc+jia+addi+ms+igad_sans-mhc/pid_uc-delange+sys-scl+lada+t1d-cooper+psc+jia+addi+ms+igad.tsv.gz"

rule draw_manhattan_plot:
    input:
        "resources/gwas/{trait}.tsv.gz"
    output:
        "results/plots/gwas/{trait}_manhattan.png"
    params:
        chrom_col = 'CHR38',
        bp_col = 'BP38',
        prin_col = 'P'
    threads: 4
    script:
        "../scripts/plot_gwas_manhattan.R"

rule draw_cfdr_v_manhattan_plot:
    input:
        results_file = "results/cfdr/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}.tsv.gz"
    output:
        "results/cfdr/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{snp_set}/{prin_trait}_{aux_traits}_it_{iteration}_manhattan.png"
    params:
        chrom_col = 'CHR38',
        bp_col = 'BP38',
        prin_col = 'P',
        aux_col = lambda wildcards: f"P.{len(wildcards.aux_traits.split('+'))}",
        v_col = lambda wildcards: f"v.{wildcards.iteration}",
        prin_label = lambda wildcards: wildcards.prin_trait,
        aux_label = lambda wildcards: wildcards.aux_traits
    threads: 4
    script:
        "../scripts/plot_back_to_back_manhattan.R"

rule run_cfdr_for_top_scoring_traits:
    input:
        [f'results/cfdr/pid_{imd}/pid_{imd}_sans-mhc/pid_{imd}_manhattan.png' for imd in ['uc-delange', 'sys-scl', 'lada', 't1d', 'psc', 'jia', 'addi', 'ms', 'igad', 'sle', 'pso', 'aster', 'pbc', 'ra', 't1d-cooper', 'uc', 'igm']]
