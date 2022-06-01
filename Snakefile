import datetime
import re
import os

wildcard_constraints:
  chr = "chr[0-9X]{1,2}"

max_time = datetime.timedelta(seconds=6*60*60)

def get_mem_mb(wildcards, threads):
    return threads * 3420

def get_permute_time(wildcards, threads):
    return str(min(max_time, datetime.timedelta(seconds=300+3*((int(wildcards.draws)/300)*3600)/int(threads))))

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
"pid_paed-all"]

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

rule join_gwas:
     input:
      A = "resources/gwas/pid.tsv.gz",
      B = "resources/gwas/{imd}.tsv.gz"
     output:
      AB = "resources/gwas/pid_{imd}/pid_{imd}_{join}/pid_{imd}_{join}.tsv.gz"
     threads: 2
     params:
         mhc = lambda wildcards: '-sans_mhc' if wildcards.join == 'sans_mhc' else ''
     group: "gps"
     shell:
      "Rscript scripts/join_gwas_stats.R -a {input.A} -b {input.B} {params.mhc} -o {output.AB} -nt {threads}"

rule make_plink_ranges:
     input:
      ("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      gwas_file = "resources/gwas/pid_{imd}/pid_{imd}_{join}/pid_{imd}_{join}.tsv.gz"
     output:
      ("resources/gwas/pid_{imd}/pid_{imd}_{join}/matching_ids/chr%d.txt" % x for x in range(1,23))
     params:
      input_dir = "resources/1000g/euro/qc",
      output_dir = "resources/gwas/pid_{imd}/pid_{imd}_{join}/matching_ids"
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
      temp("resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}.bed"),
      temp("resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}.bim"),
      temp("resources/gwas/pid_{imd}/pid_{imd}_{join}/plink/{chr}.fam")
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
      "Rscript scripts/prune_gwas.R -i {input.gwas_file}  -p {input.prune_file} -o {output} -nt {threads}"

rule unzip_pruned_joined_gwas:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv.gz"
    output:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    group: "gps"
    shell:
        "gunzip -c {input} >{output}"

rule compute_gps_for_trait_pair:
     input:
         "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
     output:
      "results/pid_{imd}/pid_{imd}_{join}_gps_value.tsv"
     threads: 1
     group: "gps"
     shell:
      "scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c pid -d {wildcards.imd} -n {threads} -p -f addEpsilon -l -o {output}"

rule compute_gps_for_trait_pair_with_naive_ecdf_algo:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    output:
        "results/pid_{imd}/pid_{imd}_{join}_gps_value_naive.tsv"
    threads: 4
    resources:
        time = 40
    group: "gps"
    shell:
        "scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c pid -d {wildcards.imd} -n {threads} -o {output}"

rule permute_trait_pair:
    input:
        "resources/gwas/pid_{imd}/pid_{imd}_{join}/pruned_pid_{imd}_{join}.tsv"
    output:
        "results/permutations/{draws}_draws/pid_{imd}_{join}.tsv"
    threads: 10
    resources:
        mem_mb = get_mem_mb,
        time = get_permute_time
    shell:
        "scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input} -o {output} -a P.A -b P.B -c {threads} -n {wildcards.draws}"

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
      gps_file = "results/pid_{imd}/pid_{imd}_{join}_gps_value.tsv",
      perm_file = "results/permutations/{draws}_draws/pid_{imd}_{join}.tsv"
    output:
      "results/pid_{imd}/pid_{imd}_{join}_{draws}_draws_gps_pvalue.tsv"
    shell:
      "Rscript scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a pid -b {wildcards.imd} -o {output}"

rule collate_gps_pvalue_data:
    input:
        [y for y in [f"results/{x}/{x}_all_3000_draws_gps_pvalue.tsv" for x in trait_pairs] if os.path.exists(y)]
    output:
        "results/all_3000_draws_gps_pvalues.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\tgps\tn\tloc\tloc.sd\tscale\tscale.sd\tshape\tshape.sd\tpval\n")
            for x in input:
                with open(x, 'r') as infile:
                    print(x)
                    trait_B = re.match('results/pid_([A-Za-z0-9-_]+)/pid_[A-Za-z0-9-_]+_all_3000_draws_gps_pvalue.tsv', x).groups()[0]
                    line = infile.readline()
                    line = infile.readline()

                outfile.write(f"pid\t{trait_B}\t{line}")

rule gps_for_all_selected_imds:
     input:
         [f"results/{x}/{x}_all_gps_value.tsv" for x in trait_pairs]+
         [f"results/{x}/{x}_all_gps_value_naive.tsv" for x in trait_pairs]

rule gps_for_all_selected_imds_sans_mhc:
    input:
        [f"results/{x}/{x}_sans_mhc_gps_value.tsv" for x in trait_pairs]+
        [f"results/{x}/{x}_sans_mhc_gps_value_naive.tsv" for x in trait_pairs]

rule gps_pvalues_for_all_selected_imds:
    input:
        [f"results/{x}/{x}_all_3000_draws_gps_pvalue.tsv" for x in trait_pairs]
