wildcard_constraints:
  chr = "chr[0-9X]{1,2}"

def get_mem_mb(wildcards, threads):
    return threads * 3420

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
      "resources/1000g/20130606_g1k_3202_samples_ped_population.txt"
     shell:
      "wget -O resources/1000g/20130606_g1k_3202_samples_ped_population.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule vcf_to_bed:
     input:
      "resources/1000g/{chr}.vcf.gz"
     output:
      "resources/1000g/{chr}.bed",
      "resources/1000g/{chr}.bim",
      "resources/1000g/{chr}.fam"
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
      "resources/1000g/euro/{chr}_euro.bed",
      "resources/1000g/euro/{chr}_euro.bim",
      "resources/1000g/euro/{chr}_euro.fam"
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
      AB = "resources/gwas/pid_{imd}/pid_{imd}.tsv.gz"
     threads: 8
     shell:
      "Rscript scripts/join_gwas_stats.R -a {input.A} -b {input.B} -o {output.AB} -nt {threads}"

rule make_plink_ranges:
     input:
      ("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      gwas_file = "resources/gwas/pid_{imd}/pid_{imd}.tsv.gz"
     output:
      ("resources/gwas/pid_{imd}/matching_ids/chr%d.txt" % x for x in range(1,23))
     params:
      input_dir = "resources/1000g/euro/qc",
      output_dir = "resources/gwas/pid_{imd}/matching_ids"
     threads: 2
     shell:
      "Rscript scripts/make_plink_ranges.R -i {input.gwas_file} -b {params.input_dir} -r chr%d_qc.bim -o {params.output_dir} -nt {threads}"

rule subset_reference:
     input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam",
      range_file = "resources/gwas/pid_{imd}/matching_ids/{chr}.txt"
     output:
      temp("resources/gwas/pid_{imd}/plink/{chr}.bed"),
      temp("resources/gwas/pid_{imd}/plink/{chr}.bim"),
      temp("resources/gwas/pid_{imd}/plink/{chr}.fam")
     params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "resources/gwas/pid_{imd}/plink/{chr}"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_prune_ranges:
     input:
      "resources/gwas/pid_{imd}/plink/{chr}.bed",
      "resources/gwas/pid_{imd}/plink/{chr}.bim",
      "resources/gwas/pid_{imd}/plink/{chr}.fam"
     output:
      "resources/gwas/pid_{imd}/prune/{chr}.prune.in",
      "resources/gwas/pid_{imd}/prune/{chr}.prune.out"
     params:
       bfile = "resources/gwas/pid_{imd}/plink/{chr}",
       prune_out = "resources/gwas/pid_{imd}/prune/{chr}"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise 1000kb 50 0.2 --out {params.prune_out}"

rule cat_prune_ranges:
     input:
      ("resources/gwas/pid_{imd}/prune/chr%d.prune.in" % x for x in range(1,23))
     output:
      "resources/gwas/pid_{imd}/prune/all.prune.in"
     shell:
      "for x in {input}; do cat $x >>{output}; done"

rule prune_gwas:
     input:
      prune_file = "resources/gwas/pid_{imd}/prune/all.prune.in",
      gwas_file = "resources/gwas/pid_{imd}/pid_{imd}.tsv.gz"
     output:
      "resources/gwas/pid_{imd}/pruned_pid_{imd}.tsv.gz"
     threads: 2
     shell:
      "Rscript scripts/prune_gwas.R -i {input.gwas_file}  -p {input.prune_file} -o {output} -nt {threads}"

rule unzip_pruned_joined_gwas:
    input:
        "resources/gwas/pid_{imd}/pruned_pid_{imd}.tsv.gz"
    output:
        "resources/gwas/pid_{imd}/pruned_pid_{imd}.tsv"
    shell:
        "gunzip -c {input} >{output}"

rule compute_gps:
     input:
      "resources/gwas/pid_{imd}/pruned_pid_{imd}.tsv"
     output:
      "results/pid_{imd}_gps.tsv"
     shell:
      "scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c pid -d {wildcards.imd} -o {output}"

rule permute_trait_pair:
    input:
        "resources/gwas/pid_{imd}/pruned_pid_{imd}.tsv"
    output:
        "results/permutations/pid_{imd}.tsv"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input} -o {output} -a P.A -b P.B -c {threads} -n 375"

rule compute_gps_p_value:
    input:
        gps_file = "results/pid_{imd}_gps.tsv",
        perm_file = "results/permutations/pid_{imd}.tsv"
    output:
        "results/pid_{imd}_gps_pvalue.tsv"
    shell:
        "Rscript scripts/compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a pid -b {wildcards.imd} -o {output}"

rule gps_for_all_selected_imds:
     input:
        "results/pid_acne_gps_pvalue.tsv",
        "results/pid_fibrom_gps_pvalue.tsv",
        "results/pid_pbc_gps_pvalue.tsv",
        "results/pid_addi_gps_pvalue.tsv",
        "results/pid_ges-dia_gps_pvalue.tsv",
        "results/pid_agran_gps_pvalue.tsv",
        "results/pid_glom_gps_pvalue.tsv",
        "results/pid_polymyal_gps_pvalue.tsv",
        "results/pid_all-conj_gps_pvalue.tsv",
        "results/pid_gout_gps_pvalue.tsv",
        "results/pid_polyp_gps_pvalue.tsv",
        "results/pid_all-urt_gps_pvalue.tsv",
        "results/pid_hyperpara_gps_pvalue.tsv",
        "results/pid_pros_gps_pvalue.tsv",
        "results/pid_ank-spond_gps_pvalue.tsv",
        "results/pid_hyperthy_gps_pvalue.tsv",
        "results/pid_psc_gps_pvalue.tsv",
        "results/pid_arthr-nos_gps_pvalue.tsv",
        "results/pid_hypothy_gps_pvalue.tsv",
        "results/pid_pso-arthr_gps_pvalue.tsv",
        "results/pid_arthr_gps_pvalue.tsv",
        "results/pid_ibd_gps_pvalue.tsv",
        "results/pid_pso_gps_pvalue.tsv",
        "results/pid_aster_gps_pvalue.tsv",
        "results/pid_ibs_gps_pvalue.tsv",
        "results/pid_ra_gps_pvalue.tsv",
        "results/pid_asthma_gps_pvalue.tsv",
        "results/pid_igad_gps_pvalue.tsv",
        "results/pid_rheum-fev_gps_pvalue.tsv",
        "results/pid_auto-dis_gps_pvalue.tsv",
        "results/pid_inf-bre_gps_pvalue.tsv",
        "results/pid_rhin_gps_pvalue.tsv",
        "results/pid_cardiomyo_gps_pvalue.tsv",
        "results/pid_inf-cer-ute_gps_pvalue.tsv",
        "results/pid_rosacea_gps_pvalue.tsv",
        "results/pid_coeliac_gps_pvalue.tsv",
        "results/pid_inf-gyn_gps_pvalue.tsv",
        "results/pid_sarc_gps_pvalue.tsv",
        "results/pid_copd_gps_pvalue.tsv",
        "results/pid_inf-ute_gps_pvalue.tsv",
        "results/pid_sjogren_gps_pvalue.tsv",
        "results/pid_crohns-med_gps_pvalue.tsv",
        "results/pid_jia_gps_pvalue.tsv",
        "results/pid_sle_gps_pvalue.tsv",
        "results/pid_crohns_gps_pvalue.tsv",
        "results/pid_lada_gps_pvalue.tsv",
        "results/pid_spondylo_gps_pvalue.tsv",
        "results/pid_derm-ecz_gps_pvalue.tsv",
        "results/pid_licp_gps_pvalue.tsv",
        "results/pid_still_gps_pvalue.tsv",
        "results/pid_desq_gps_pvalue.tsv",
        "results/pid_licsa_gps_pvalue.tsv",
        "results/pid_sulf-all_gps_pvalue.tsv",
        "results/pid_divert_gps_pvalue.tsv",
        "results/pid_misc-blood_gps_pvalue.tsv",
        "results/pid_sys-scl_gps_pvalue.tsv",
        "results/pid_drug-all_gps_pvalue.tsv",
        "results/pid_misc-immune_gps_pvalue.tsv",
        "results/pid_t1d_gps_pvalue.tsv",
        "results/pid_dyschr-vit_gps_pvalue.tsv",
        "results/pid_ms_gps_pvalue.tsv",
        "results/pid_t2d_gps_pvalue.tsv",
        "results/pid_endocard_gps_pvalue.tsv",
        "results/pid_myal_gps_pvalue.tsv",
        "results/pid_thy_gps_pvalue.tsv",
        "results/pid_endomet_gps_pvalue.tsv",
        "results/pid_myas_gps_pvalue.tsv",
        "results/pid_uc_gps_pvalue.tsv",
        "results/pid_ent-col_gps_pvalue.tsv",
        "results/pid_pad_gps_pvalue.tsv",
        "results/pid_eos-eso_gps_pvalue.tsv",
        "results/pid_paed-all_gps_pvalue.tsv"

"""
rule maf:
     input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam"
     output:
      "resources/1000g/euro/qc/maf/{chr}.frq"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/qc/{wildcards.chr}_qc --freq --out resources/1000g/euro/qc/maf/{wildcards.chr}"

rule ld:
     input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam"
     output:
      "resources/1000g/euro/qc/ld/{chr}.ld"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/qc/{wildcards.chr}_qc --r2 --ld-window-r2 0.2 --ld-window-kb 1000 --out resources/1000g/euro/qc/ld/{wildcards.chr}"

rule prune_whole_set:
     input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam"
     output:
      "resources/1000g/euro/qc/prune_whole/{chr}.prune.in",
      "resources/1000g/euro/qc/prune_whole/{chr}.prune.out",
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/qc/{wildcards.chr}_qc --indep-pairwise 1000kb 50 0.2 --out resources/1000g/euro/qc/prune_whole/{wildcards.chr}"

rule all:
     input:
      ("resources/gwas/pid_{imd}/plink/pruned/chr%d.bed" % x for x in range(1,23)),
      ("resources/gwas/pid_{imd}/plink/pruned/chr%d.bim" % x for x in range(1,23)),
      ("resources/gwas/pid_{imd}/plink/pruned/chr%d.fam" % x for x in range(1,23))
     output:
      "resources/gwas/pid_{imd}/done.txt"
     shell:
      "touch {output}"

rule prune_subset:
     input:
      "resources/gwas/pid_{imd}/plink/{chr}.bed",
      "resources/gwas/pid_{imd}/plink/{chr}.bim",
      "resources/gwas/pid_{imd}/plink/{chr}.fam",
      range_file = "resources/gwas/pid_{imd}/prune/{chr}.prune.in"
     output:
      "resources/gwas/pid_{imd}/plink/pruned/{chr}.bed",
      "resources/gwas/pid_{imd}/plink/pruned/{chr}.bim",
      "resources/gwas/pid_{imd}/plink/pruned/{chr}.fam",
     params:
      bfile = "resources/gwas/pid_{imd}/plink/{chr}",
      out = "resources/gwas/pid_{imd}/plink/pruned/{chr}"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule subset_ld:
     input:
      ("resources/1000g/euro/qc/ld/chr%d.ld" % x for x in range(1,23)),
      "resources/1000g/euro/qc/ld/chrX.ld",
      "resources/gwas/pid_{imd}/pid_{imd}.tsv.gz"
     output:
      ("resources/gwas/pid_{imd}/ld/chr%d.ld" % x for x in range(1,23)),
      "resources/gwas/pid_{imd}/ld/chrX.ld"
     shell:
      "Rscript scripts/subset_ld_with_merged_gwas.R -i resources/gwas/pid_{wildcards.imd}  -ld resources/1000g/euro/qc/ld -o resources/gwas/pid_{wildcards.imd}/ld"

rule prune_all:
     input:
      ("resources/1000g/euro/qc/prune_whole/chr%d.prune.in" % x for x in range(1,23)),
      "resources/1000g/euro/qc/prune_whole/chrX.prune.in"

rule ld_all:
     input:
      ("resources/1000g/euro/qc/ld/chr%d.ld" % x for x in range(1,23)),
      "resources/1000g/euro/qc/ld/chrX.ld"
"""
