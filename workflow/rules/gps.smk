rule test_perturbation_parameters:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv"
    output:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/perturbation/pid_{imd}_{snp_set}_{epsilon}_{no_of_perts}.tsv"
    params:
        epsilon = lambda wildcards: float(wildcards.epsilon),
        no_of_perts = lambda wildcards: int(wildcards.no_of_perts)
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteProblem -i {input} -a P.A -b P.B -e {params.epsilon} -p {params.no_of_perts} &>{output}"

rule compute_gps_for_trait_pair:
     input:
         "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv"
     output:
      "results/pid_{imd}/pid_{imd}_{snp_set}_gps_value.tsv"
     log:
      "results/pid_{imd}/pid_{imd}_{snp_set}_gps_value.log"
     params:
         no_of_pert_iterations = lambda wildcards: 300 if wildcards.imd == 'psc' else 100,
         epsilon_multiple = lambda wildcards: 10.0 if wildcards.imd == 'psc' else 2.0
     threads: 1
     group: "gps"
     shell:
      "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c pid -d {wildcards.imd} -n {threads} -p {params.no_of_pert_iterations} -e {params.epsilon_multiple} -l -o {output} -g {log}"

rule compute_gps_for_trait_pair_with_naive_ecdf_algo:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv"
    output:
        "results/pid_{imd}/pid_{imd}_{snp_set}_gps_value_naive.tsv"
    log:
        "results/pid_{imd}/pid_{imd}_{snp_set}_gps_value_naive.log"
    params:
        no_of_pert_iterations = 0
    threads: 8
    resources:
        time = 120
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c pid -d {wildcards.imd} -n {threads} -p {params.no_of_pert_iterations} -o {output} -g {log}"

rule get_missing_rows:
    input:
        gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv",
        log_file = "results/pid_{imd}/pid_{imd}_{snp_set}_gps_value_naive.log"
    output:
        "results/pid_{imd}/pid_{imd}_{snp_set}_missing_rows.tsv"
    script:
        "../scripts/get_missing_rows.R"

rule permute_trait_pair:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv"
    output:
        "results/permutations/{draws}_draws/pid_{imd}_{snp_set}.tsv"
    threads: 10
    resources:
        mem_mb = get_mem_mb,
        time = get_permute_time,
        no_of_pert_iterations = 100
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input} -o {output} -a P.A -b P.B -c {threads} -n {wildcards.draws} -p {params.no_of_pert_iterations}"

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
      gps_file = "results/pid_{imd}/pid_{imd}_{snp_set}_gps_value.tsv",
      perm_file = "results/permutations/{draws}_draws/pid_{imd}_{snp_set}.tsv"
    output:
      "results/pid_{imd}/pid_{imd}_{snp_set}_{draws}_draws_gps_pvalue.tsv"
    group: "gps"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a pid -b {wildcards.imd} -o {output}"

rule fit_gev_and_compute_gps_pvalue_for_impermutable_trait_pair:
    input:
        gps_file = "results/pid_{imd}/pid_{imd}_{snp_set}_gps_value_naive.tsv",
        collated_file = "results/{snp_set}_{draws}_draws_permutable_gps_pvalues.tsv"
    output:
        all_pvalues_file = "results/pid_{imd}/pid_{imd}_{snp_set}_{draws}_draws_gps_imputed_pvalues.tsv",
        median_pvalue_file = "results/pid_{imd}/pid_{imd}_{snp_set}_{draws}_draws_gps_pvalue_median_imputed.tsv"
    params:
        trait_name = lambda wildcards: wildcards.imd
    group: "gps"
    script:
        "../scripts/fit_gev_and_compute_gps_pvalue_for_imperm.R"

rule collate_naive_gps_value_data:
    input:
        [f"results/{x}/{x}_{{snp_set}}_gps_value_naive.tsv" for x in trait_pairs]
    output:
        "results/{snp_set}_naive_gps_values.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tgps\n")
            for x in input:
                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                    outfile.write(f"{line}")

rule collate_permutable_gps_pvalue_data:
    input:
        [f"results/{x}/{x}_{{snp_set}}_{{draws}}_draws_gps_pvalue.tsv" for x in permutable_trait_pairs]
    output:
        "results/{snp_set}_{draws}_draws_permutable_gps_pvalues.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tgps\tn\tloc\tloc.sd\tscale\tscale.sd\tshape\tshape.sd\tpval\n")
            for x in input:
                with open(x, 'r') as infile:
                    trait_B = re.match(f'results/pid_([A-Za-z0-9-_]+)/pid_[A-Za-z0-9-_]+_{wildcards.snp_set}_{wildcards.draws}_draws_gps_pvalue.tsv', x).groups()[0]
                    line = infile.readline()
                    line = infile.readline()

                outfile.write(f"pid\t{trait_B}\t{line}")

rule collate_gps_pvalue_data:
    input:
        [f"results/{x}/{x}_{{snp_set}}_{{draws}}_draws_gps_pvalue.tsv" for x in permutable_trait_pairs]+
        [f"results/{x}/{x}_{{snp_set}}_{{draws}}_draws_gps_pvalue_median_imputed.tsv" for x in impermutable_trait_pairs]
    output:
        "results/{snp_set}_{draws}_draws_gps_pvalues.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tgps\tn\tloc\tloc.sd\tscale\tscale.sd\tshape\tshape.sd\tpval\tpermutable\n")
            for x in input:
                with open(x, 'r') as infile:
                    trait_B = re.match(f'results/pid_([A-Za-z0-9-_]+)/pid_[A-Za-z0-9-_]+_{wildcards.snp_set}_{wildcards.draws}_draws.+', x).groups()[0]
                    line = infile.readline()
                    line = infile.readline().strip('\n')

                    permutable = f"pid_{trait_B}" in permutable_trait_pairs

                    # TODO fix omission of trait names in script which generates permutable pvalue output file
                    if permutable:
                        outfile.write(f"pid\t{trait_B}\t{line}\t{permutable}\n")
                    else:
                        outfile.write(f"{line}\t{permutable}\n")

rule perturb_pvalues:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_pid_{imd}_{snp_set}.tsv"
    output:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_perturbed_pid_{imd}_{snp_set}.tsv"
    params:
        no_of_pert_iterations = lambda wildcards: 300 if wildcards.imd == 'psc' else 100,
        epsilon_multiple = lambda wildcards: 10.0 if wildcards.imd == 'psc' else 2.0,
        trait_b = lambda wildcards: wildcards.imd
    threads: 1
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/perturbPvalues -i {input} -a P.A -b P.B -c pid -d {params.trait_b} -e {params.epsilon_multiple} -p {params.no_of_pert_iterations} &>{output}"

rule compute_hoeffdings_for_trait_pair:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{snp_set}/pruned_perturbed_pid_{imd}_{snp_set}.tsv"
    output:
        "results/pid_{imd}/pid_{imd}_{snp_set}_hoeffdings.tsv"
    params:
        trait_A_col = 'pid',
        trait_B_col = lambda wildcards: wildcards.imd
    threads: 2
    group: "gps"
    script:
        "../scripts/compute_hoeffdings.R"

rule collate_hoeffdings_results:
    input:
        [f"results/{x}/{x}_{{snp_set}}_hoeffdings.tsv" for x in trait_pairs]
    output:
        "results/{snp_set}_hoeffdings.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tn\tDn\tscaled\tpval\n")
            for x in input:
                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                    outfile.write(line)
