import datetime
import re
import os

wildcard_constraints:
  chr = "chr[0-9X]{1,2}"

max_time = 240

def get_mem_mb(wildcards, threads):
    return threads * 3420

def get_permute_time(wildcards, threads):
    return min(max_time, 300+3*((int(wildcards.draws)/300)*3600)/int(threads))

trait_pairs = ["pid_acne",
"pid_fibrom",
"pid_pbc",
"pid_addi",
"pid_ges-dia",
"pid_agran",
"pid_glom",
"pid_polymyal",
"pid_all-conj",
"pid_gout",
"pid_polyp",
"pid_all-urt",
"pid_hyperpara",
"pid_pros",
"pid_ank-spond",
"pid_hyperthy",
"pid_psc",
"pid_arthr-nos",
"pid_hypothy",
"pid_pso-arthr",
"pid_arthr",
"pid_ibd",
"pid_pso",
"pid_aster",
"pid_ibs",
"pid_ra",
"pid_asthma",
"pid_igad",
"pid_rheum-fev",
"pid_auto-dis",
"pid_inf-bre",
"pid_rhin",
"pid_cardiomyo",
"pid_inf-cer-ute",
"pid_rosacea",
"pid_coeliac",
"pid_inf-gyn",
"pid_sarc",
"pid_copd",
"pid_inf-ute",
"pid_sjogren",
"pid_crohns-med",
"pid_jia",
"pid_sle",
"pid_crohns",
"pid_lada",
"pid_spondylo",
"pid_derm-ecz",
"pid_licp",
"pid_still",
"pid_desq",
"pid_licsa",
"pid_sulf-all",
"pid_divert",
"pid_misc-blood",
"pid_sys-scl",
"pid_drug-all",
"pid_misc-immune",
"pid_t1d",
"pid_dyschr-vit",
"pid_ms",
"pid_t2d",
"pid_endocard",
"pid_myal",
"pid_thy",
"pid_endomet",
"pid_myas",
"pid_uc",
"pid_ent-col",
"pid_pad",
"pid_eos-eso",
"pid_paed-all",
"pid_igg",
"pid_iga",
"pid_igm",
"pid_t1d-cooper",
"pid_uc-delange"]

impermutable_trait_pairs = ['pid_ra', 'pid_psc', 'pid_t1d-cooper', 'pid_uc-delange']

permutable_trait_pairs = [x for x in trait_pairs if x not in impermutable_trait_pairs]

#include: 'sumher.smk'

rule download_1000g_genotype_data:
     output:
      "resources/1000g/{chr}.vcf.gz"
     run:
      if wildcards.chr == 'chrX':
         shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz")
      else:
        shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz")

rule download_1000g_sample_metadata:
     output:
      temp("resources/1000g/20130606_g1k_3202_samples_ped_population.txt")
     shell:
      "wget -O resources/1000g/20130606_g1k_3202_samples_ped_population.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule vcf_to_bed:
     input:
      "resources/1000g/{chr}.vcf.gz"
     output:
      temp("resources/1000g/{chr}.bed"),
      temp("resources/1000g/{chr}.bim"),
      temp("resources/1000g/{chr}.fam")
     threads: 8
     resources:
         mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --vcf resources/1000g/{wildcards.chr}.vcf.gz --make-bed --out resources/1000g/{wildcards.chr}"

rule make_euro_fam:
     input:
      "resources/1000g/20130606_g1k_3202_samples_ped_population.txt"
     output:
      "resources/1000g/euro.fam"
     shell:
      "Rscript scripts/get_euro_fam_file.R --ped_file {input} --output_file {output}"

rule get_euro_samples:
     input:
      "resources/1000g/{chr}.bed",
      "resources/1000g/{chr}.bim",
      "resources/1000g/{chr}.fam",
      fam_file = "resources/1000g/euro.fam"
     output:
      temp("resources/1000g/euro/{chr}_euro.bed"),
      temp("resources/1000g/euro/{chr}_euro.bim"),
      temp("resources/1000g/euro/{chr}_euro.fam")
     params:
      bfile = "resources/1000g/{chr}",
      out = "resources/1000g/euro/{chr}_euro"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
         "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --keep {input.fam_file} --make-bed --silent --out {params.out}"

rule qc:
     input:
      "resources/1000g/euro/{chr}_euro.bed",
      "resources/1000g/euro/{chr}_euro.bim",
      "resources/1000g/euro/{chr}_euro.fam"
     output:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam"
     params:
      bfile = "resources/1000g/euro/{chr}_euro",
      out = "resources/1000g/euro/qc/{chr}_qc",
     threads: 8
     resources:
       mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --make-bed --silent --out {params.out}"

rule join_multiple_gwas:
    input:
        principal_trait_gwas_file = "resources/gwas/{prin_trait}.tsv.gz",
        auxiliary_trait_gwas_files = lambda wildcards: [f"resources/gwas/{x}.tsv.gz" for x in wildcards.aux_traits.split('+')],
        pruned_auxiliary_trait_gwas_files = lambda wildcards: [f"resources/gwas/{wildcards.prin_trait}_{x}/{wildcards.prin_trait}_{x}_{wildcards.join}/pruned_{wildcards.prin_trait}_{x}_{wildcards.join}.tsv.gz" for x in wildcards.aux_traits.split('+')]
    output:
        "resources/gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{join}/{prin_trait}_{aux_traits}_with_prune_{join}.tsv.gz"
    threads: 4
    params:
        mhc = lambda wildcards: '-sans_mhc' if wildcards.join == 'sans-mhc' else '',
        bp_col = 'BP38',
        chr_col = 'CHR38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P'
    resources:
        time = 10
    group: "gps"
    script:
        "scripts/join_multiple_gwas_stats.R"

rule make_plink_ranges:
     input:
      ("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      gwas_file = "resources/gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{join}/{trait_A}_{trait_B}_{join}.tsv.gz"
     output:
      ("resources/gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{join}/matching_ids/chr%d.txt" % x for x in range(1,23))
     params:
        input_dir = "resources/1000g/euro/qc",
        output_dir = "resources/gwas/{trait_A}_{trait_B}/{trait_A}_{trait_B}_{join}/matching_ids"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "Rscript scripts/make_plink_ranges.R -i {input.gwas_file} -b {params.input_dir} -r chr%d_qc.bim -o {params.output_dir} -nt {threads}"

rule subset_reference:
     input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam",
      range_file = "resources/gwas/pid_{imd}/pid_{imd}_{join}/matching_ids/{chr}.txt"
     output:
      temp("resources/gwas/pid_{imd,[^\+]+}/pid_{imd}_{join}/plink/{chr}.bed"),
      temp("resources/gwas/pid_{imd,[^\+]+}/pid_{imd}_{join}/plink/{chr}.bim"),
      temp("resources/gwas/pid_{imd,[^\+]+}/pid_{imd}_{join}/plink/{chr}.fam")
     params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_pruned_ranges:
     input:
      "resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}.bed",
      "resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}.bim",
      "resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}.fam"
     output:
      "resources/gwas/pid_{imd}/pid_{imd}_{join}/prune/{chr}.prune.in",
      "resources/gwas/pid_{imd}/pid_{imd}_{join}/prune/{chr}.prune.out"
     params:
       bfile = "resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}",
       prune_out = "resources/gwas/pid_{imd}/pid_{imd}_{join}/prune/{chr}"
     threads: 1
     resources:
        mem_mb=get_mem_mb
     group: "gps"
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise 1000kb 50 0.2 --out {params.prune_out}"

rule cat_pruned_ranges:
     input:
      ("resources/gwas/pid_{imd}/pid_{imd}_{join}/prune/chr%d.prune.in" % x for x in range(1,23))
     output:
      "resources/gwas/pid_{imd}/pid_{imd}_{join}/prune/all.prune.in"
     threads: 1
     group: "gps"
     shell:
      "for x in {input}; do cat $x >>{output}; done"

rule prune_gwas:
     input:
      prune_file = "resources/gwas/pid_{imd}/pid_{imd}_{join}/prune/all.prune.in",
      gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{join}/pid_{imd}_{join}.tsv.gz"
     output:
      "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv.gz"
     threads: 1
     group: "gps"
     shell:
      "Rscript scripts/prune_gwas.R -i {input.gwas_file} -p {input.prune_file} -o {output} -nt {threads}"

rule unzip_pruned_joined_gwas:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv.gz"
    output:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    group: "gps"
    shell:
        "gunzip -c {input} >{output}"

rule test_perturbation_parameters:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    output:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/perturbation/pid_{imd}_{join}_{epsilon}_{no_of_perts}.tsv"
    params:
        epsilon = lambda wildcards: float(wildcards.epsilon),
        no_of_perts = lambda wildcards: int(wildcards.no_of_perts)
    group: "gps"
    shell:
        "scripts/gps_cpp/build/apps/permuteProblem -i {input} -a P.A -b P.B -e {params.epsilon} -p {params.no_of_perts} &>{output}"

rule compute_gps_for_trait_pair:
     input:
         "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
     output:
      "results/pid_{imd}/pid_{imd}_{join}_gps_value.tsv"
     log:
      "results/pid_{imd}/pid_{imd}_{join}_gps_value.log"
     params:
         no_of_pert_iterations = lambda wildcards: 300 if wildcards.imd == 'psc' else 100,
         epsilon_multiple = lambda wildcards: 10.0 if wildcards.imd == 'psc' else 2.0
     threads: 1
     group: "gps"
     shell:
      "scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c pid -d {wildcards.imd} -n {threads} -p {params.no_of_pert_iterations} -e {params.epsilon_multiple} -l -o {output} -g {log}"

rule compute_gps_for_trait_pair_with_naive_ecdf_algo:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    output:
        "results/pid_{imd}/pid_{imd}_{join}_gps_value_naive.tsv"
    log:
        "results/pid_{imd}/pid_{imd}_{join}_gps_value_naive.log"
    params:
        no_of_pert_iterations = 0
    threads: 8
    resources:
        time = 120
    group: "gps"
    shell:
        "scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c pid -d {wildcards.imd} -n {threads} -p {params.no_of_pert_iterations} -o {output} -g {log}"

rule get_missing_rows:
    input:
        gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv",
        log_file = "results/pid_{imd}/pid_{imd}_{join}_gps_value_naive.log"
    output:
        "results/pid_{imd}/pid_{imd}_{join}_missing_rows.tsv"
    script:
        "scripts/get_missing_rows.R"

rule permute_trait_pair:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    output:
        "results/permutations/{draws}_draws/pid_{imd}_{join}.tsv"
    threads: 10
    resources:
        mem_mb = get_mem_mb,
        time = get_permute_time
    group: "gps"
    shell:
        "scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input} -o {output} -a P.A -b P.B -c {threads} -n {wildcards.draws}"

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
      gps_file = "results/pid_{imd}/pid_{imd}_{join}_gps_value.tsv",
      perm_file = "results/permutations/{draws}_draws/pid_{imd}_{join}.tsv"
    output:
      "results/pid_{imd}/pid_{imd}_{join}_{draws}_draws_gps_pvalue.tsv"
    group: "gps"
    shell:
      "Rscript scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a pid -b {wildcards.imd} -o {output}"

rule fit_gev_and_compute_gps_pvalue_for_impermutable_trait_pair:
    input:
        gps_file = "results/pid_{imd}/pid_{imd}_{join}_gps_value_naive.tsv",
        collated_file = "results/{join}_{draws}_draws_permutable_gps_pvalues.tsv"
    output:
        all_pvalues_file = "results/pid_{imd}/pid_{imd}_{join}_{draws}_draws_gps_imputed_pvalues.tsv",
        median_pvalue_file = "results/pid_{imd}/pid_{imd}_{join}_{draws}_draws_gps_pvalue_median_imputed.tsv"
    params:
        trait_name = lambda wildcards: wildcards.imd
    group: "gps"
    script:
        "scripts/fit_gev_and_compute_gps_pvalue_for_imperm.R"

rule collate_naive_gps_value_data:
    input:
        [f"results/{x}/{x}_{{join}}_gps_value_naive.tsv" for x in trait_pairs]
    output:
        "results/{join}_naive_gps_values.tsv"
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
        [f"results/{x}/{x}_{{join}}_{{draws}}_draws_gps_pvalue.tsv" for x in permutable_trait_pairs]
    output:
        "results/{join}_{draws}_draws_permutable_gps_pvalues.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tgps\tn\tloc\tloc.sd\tscale\tscale.sd\tshape\tshape.sd\tpval\n")
            for x in input:
                with open(x, 'r') as infile:
                    trait_B = re.match(f'results/pid_([A-Za-z0-9-_]+)/pid_[A-Za-z0-9-_]+_{wildcards.join}_{wildcards.draws}_draws_gps_pvalue.tsv', x).groups()[0]
                    line = infile.readline()
                    line = infile.readline()

                outfile.write(f"pid\t{trait_B}\t{line}")

rule collate_gps_pvalue_data:
    input:
        [f"results/{x}/{x}_{{join}}_{{draws}}_draws_gps_pvalue.tsv" for x in permutable_trait_pairs]+
        [f"results/{x}/{x}_{{join}}_{{draws}}_draws_gps_pvalue_median_imputed.tsv" for x in impermutable_trait_pairs]
    output:
        "results/{join}_{draws}_draws_gps_pvalues.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tgps\tn\tloc\tloc.sd\tscale\tscale.sd\tshape\tshape.sd\tpval\tpermutable\n")
            for x in input:
                with open(x, 'r') as infile:
                    trait_B = re.match(f'results/pid_([A-Za-z0-9-_]+)/pid_[A-Za-z0-9-_]+_{wildcards.join}_{wildcards.draws}_draws.+', x).groups()[0]
                    line = infile.readline()
                    line = infile.readline().strip('\n')

                    permutable = f"pid_{trait_B}" in permutable_trait_pairs

                    # TODO fix omission of trait names in script which generates permutable pvalue output file
                    if permutable:
                        outfile.write(f"pid\t{trait_B}\t{line}\t{permutable}\n")
                    else:
                        outfile.write(f"{line}\t{permutable}\n")

rule run_pairwise_cfdr:
    input:
        gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{join}/pid_{imd}_{join}.tsv.gz",
        pruned_gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv.gz"
    output:
        results_file = "results/cfdr/pid_{imd}/pid_{imd}_{join}/pid_{imd}.tsv.gz"
    benchmark:
        "benchmarks/run_pairwise_cfdr/pid_{imd}_{join}_benchmark.txt"
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
         "scripts/pairwise_cfdr.R"

rule run_iterative_cfdr:
    input:
        "resources/gwas/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{join}/{prin_trait}_{aux_traits}_with_prune_{join}.tsv.gz"
    output:
        "results/cfdr/{prin_trait}_{aux_traits}/{prin_trait}_{aux_traits}_{join}/{prin_trait}_{aux_traits}.tsv.gz"
    benchmark:
        "benchmarks/run_iterative_cfdr/{prin_trait}_{aux_traits}_{join}_benchmark.txt"
    threads: 6
             # TODO write function to handle times
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
        p_threshold = 1e-2
    script:
        "scripts/iterative_cfdr.R"

rule run_iterative_cfdr_for_top_traits_per_gps:
    input:
        "results/cfdr/pid_uc-delange+sys-scl+lada+t1d-cooper+psc+jia+addi+ms+igad/pid_uc-delange+sys-scl+lada+t1d-cooper+psc+jia+addi+ms+igad_sans-mhc/pid_uc-delange+sys-scl+lada+t1d-cooper+psc+jia+addi+ms+igad.tsv.gz"


rule draw_pairwise_cfdr_v_manhattan_plot:
    input:
        results_file = "results/cfdr/pid_{imd}/pid_{imd}_{join}/pid_{imd}.tsv.gz"
    output:
        "results/cfdr/pid_{imd}/pid_{imd}_{join}/pid_{imd}_manhattan.png"
    params:
        chrom_col = 'CHR38',
        bp_col = 'BP38',
        prin_col = 'P.A',
        aux_col = 'P.B',
        v_col = 'v.B',
        prin_label = 'PID',
        aux_label = lambda wildcards: wildcards.imd
    threads: 4
    script:
        "scripts/plot_back_to_back_manhattan.R"

rule run_cfdr_for_top_scoring_traits:
    input:
        [f'results/cfdr/pid_{imd}/pid_{imd}_sans-mhc/pid_{imd}_manhattan.png' for imd in ['uc-delange', 'sys-scl', 'lada', 't1d', 'psc', 'jia', 'addi', 'ms', 'igad', 'sle', 'pso', 'aster', 'pbc', 'ra', 't1d-cooper', 'uc', 'igm']]

rule perturb_pvalues:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    output:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_perturbed_pid_{imd}_{join}.tsv"
    params:
        no_of_pert_iterations = lambda wildcards: 300 if wildcards.imd == 'psc' else 100,
        epsilon_multiple = lambda wildcards: 10.0 if wildcards.imd == 'psc' else 2.0,
        trait_b = lambda wildcards: wildcards.imd
    threads: 1
    group: "gps"
    shell:
        "scripts/gps_cpp/build/apps/perturbPvalues -i {input} -a P.A -b P.B -c pid -d {params.trait_b} -e {params.epsilon_multiple} -p {params.no_of_pert_iterations} &>{output}"

rule compute_hoeffdings_for_trait_pair:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_perturbed_pid_{imd}_{join}.tsv"
    output:
        "results/pid_{imd}/pid_{imd}_{join}_hoeffdings.tsv"
    params:
        trait_A_col = 'pid',
        trait_B_col = lambda wildcards: wildcards.imd
    threads: 2
    group: "gps"
    script:
        "scripts/compute_hoeffdings.R"

rule collate_hoeffdings_results:
    input:
        [f"results/{x}/{x}_{{join}}_hoeffdings.tsv" for x in trait_pairs]
    output:
        "results/{join}_hoeffdings.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tn\tDn\tscaled\tpval\n")
            for x in input:
                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                    outfile.write(line)
