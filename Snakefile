wildcard_constraints:
  chr = "chr[0-9X]{1,2}"

rule download_1000g_genotype_data:
     output:
      "1000g/{chr}.vcf.gz"
     run:
      if wildcards.chr == 'chrX':
         shell("wget -O 1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz")
      else:
        shell("wget -O 1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz")

rule download_1000g_sample_metadata:
     output:
      "1000g/20130606_g1k_3202_samples_ped_population.txt"
     shell:
      "wget -O 1000g/20130606_g1k_3202_samples_ped_population.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule vcf_to_bed:
     input:
      "1000g/{chr}.vcf.gz"
     output:
      "1000g/{chr}.bed",
      "1000g/{chr}.bim",
      "1000g/{chr}.fam"
     shell:
      "plink --memory 16000 --threads 8 --vcf 1000g/{wildcards.chr}.vcf.gz --make-bed --out 1000g/{wildcards.chr}"

rule make_euro_fam:
     input:
      "1000g/20130606_g1k_3202_samples_ped_population.txt"
     output:
      "1000g/euro.fam"
     shell:
      "Rscript scripts/get_euro_fam_file.R"

rule get_euro_samples:
     input:
      "1000g/{chr}.bed",
      "1000g/{chr}.bim",
      "1000g/{chr}.fam",
      "1000g/euro.fam"
     output:
      "1000g/euro/{chr}_euro.bed",
      "1000g/euro/{chr}_euro.bim",
      "1000g/euro/{chr}_euro.fam"
     shell:
      "plink --memory 16000 --threads 8 --bfile 1000g/{wildcards.chr} --keep 1000g/euro.fam --make-bed --silent --out 1000g/euro/{wildcards.chr}_euro"

rule qc:
     input:
      "1000g/euro/{chr}_euro.bed",
      "1000g/euro/{chr}_euro.bim",
      "1000g/euro/{chr}_euro.fam"
     output:
      "1000g/euro/qc/{chr}_qc.bed",
      "1000g/euro/qc/{chr}_qc.bim",
      "1000g/euro/qc/{chr}_qc.fam"
     shell:
      "plink --memory 16000 --threads 8 --bfile 1000g/euro/{wildcards.chr}_euro --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --make-bed --silent --out 1000g/euro/qc/{wildcards.chr}_qc"

rule join_gwas:
     input:
      A = "gwas/pid.tsv.gz",
      B = "gwas/{imd}.tsv.gz"
     output:
      AB = "gwas/pid_{imd}/pid_{imd}.tsv.gz"
     threads: 8
     shell:
      "Rscript scripts/join_gwas_stats.R -a {input.A} -b {input.B} -o {output.AB} -nt 8"

rule make_plink_ranges:
     input:
      ("1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      gwas_file = "gwas/pid_{imd}/pid_{imd}.tsv.gz"
     output:
      ("gwas/pid_{imd}/matching_ids/chr%d.txt" % x for x in range(1,23))
     params:
      input_dir = "1000g/euro/qc",
      output_dir = "gwas/pid_{imd}/matching_ids",
     threads: 2
     shell:
      "Rscript scripts/make_plink_ranges.R -i {input.gwas_file} -b {params.input_dir} -r chr%d_qc.bim -o {params.output_dir} -nt {threads}"

rule subset_reference:
     input:
      "1000g/euro/qc/{chr}_qc.bed",
      "1000g/euro/qc/{chr}_qc.bim",
      "1000g/euro/qc/{chr}_qc.fam",
      range_file = "gwas/pid_{imd}/matching_ids/{chr}.txt"
     output:
      "gwas/pid_{imd}/plink/{chr}.bed",
      "gwas/pid_{imd}/plink/{chr}.bim",
      "gwas/pid_{imd}/plink/{chr}.fam"
     params:
      bfile = "1000g/euro/qc/{chr}_qc",
      out = "gwas/pid_{imd}/plink/{chr}"
     threads: 8
     shell:
      "plink --memory 16000 --threads 8 --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_prune_ranges:
     input:
      "gwas/pid_{imd}/plink/{chr}.bed",
      "gwas/pid_{imd}/plink/{chr}.bim",
      "gwas/pid_{imd}/plink/{chr}.fam"
     output:
      "gwas/pid_{imd}/prune/{chr}.prune.in",
      "gwas/pid_{imd}/prune/{chr}.prune.out"
     params:
      bfile = "gwas/pid_{imd}/plink/{chr}",
      prune_out = "gwas/pid_{imd}/prune/{chr}"
     threads: 8
     shell:
      "plink --memory 16000 --threads 8 --bfile {params.bfile} --indep-pairwise 1000kb 50 0.2 --out {params.prune_out}"

rule cat_prune_ranges:
     input:
      ("gwas/pid_{imd}/prune/chr%d.prune.in" % x for x in range(1,23))
     output:
      "gwas/pid_{imd}/prune/all.prune.in"
     shell:
      "for x in {input}; do cat $x >>{output}; done"

rule prune_gwas:
     input:
      prune_file = "gwas/pid_{imd}/prune/all.prune.in",
      gwas_file = "gwas/pid_{imd}/pid_{imd}.tsv.gz"
     output:
      "gwas/pid_{imd}/pruned_pid_{imd}.tsv.gz"
     threads: 2
     shell:
      "Rscript scripts/prune_gwas.R -i {input.gwas_file}  -p {input.prune_file} -o {output} -nt {threads}"

rule compute_gps:
     input:
      pruned_gwas_file = "gwas/pid_{imd}/pruned_pid_{imd}.tsv.gz"
     output:
      "gwas/pid_{imd}/gps.tsv"
     threads: 8
     shell:
      "Rscript scripts/run_gps.R -i {input.pruned_gwas_file} -o {output} -nt {threads}"

rule gps_for_all_selected_imds:
     input:
      ["gwas/pid_acne/gps.tsv", "gwas/pid_fibrom/gps.tsv", "gwas/pid_pbc/gps.tsv",
"gwas/pid_addi/gps.tsv", "gwas/pid_ges-dia/gps.tsv",
"gwas/pid_agran/gps.tsv", "gwas/pid_glom/gps.tsv", "gwas/pid_polymyal/gps.tsv",
"gwas/pid_all-conj/gps.tsv", "gwas/pid_gout/gps.tsv", "gwas/pid_polyp/gps.tsv",
"gwas/pid_all-urt/gps.tsv", "gwas/pid_hyperpara/gps.tsv", "gwas/pid_pros/gps.tsv",
"gwas/pid_ank-spond/gps.tsv", "gwas/pid_hyperthy/gps.tsv", "gwas/pid_psc/gps.tsv",
"gwas/pid_arthr-nos/gps.tsv", "gwas/pid_hypothy/gps.tsv", "gwas/pid_pso-arthr/gps.tsv",
"gwas/pid_arthr/gps.tsv", "gwas/pid_ibd/gps.tsv", 	 "gwas/pid_pso/gps.tsv",
"gwas/pid_aster/gps.tsv", "gwas/pid_ibs/gps.tsv", 	 "gwas/pid_ra/gps.tsv",
"gwas/pid_asthma/gps.tsv", "gwas/pid_igad/gps.tsv", "gwas/pid_rheum-fev/gps.tsv",
"gwas/pid_auto-dis/gps.tsv", "gwas/pid_inf-bre/gps.tsv", "gwas/pid_rhin/gps.tsv",
"gwas/pid_cardiomyo/gps.tsv", "gwas/pid_inf-cer-ute/gps.tsv", "gwas/pid_rosacea/gps.tsv",
"gwas/pid_coeliac/gps.tsv", "gwas/pid_inf-gyn/gps.tsv", "gwas/pid_sarc/gps.tsv",
"gwas/pid_copd/gps.tsv", "gwas/pid_inf-ute/gps.tsv", "gwas/pid_sjogren/gps.tsv",
"gwas/pid_crohns-med/gps.tsv", "gwas/pid_jia/gps.tsv", 	 "gwas/pid_sle/gps.tsv",
"gwas/pid_crohns/gps.tsv", "gwas/pid_lada/gps.tsv", "gwas/pid_spondylo/gps.tsv",
"gwas/pid_derm-ecz/gps.tsv", "gwas/pid_licp/gps.tsv", "gwas/pid_still/gps.tsv",
"gwas/pid_desq/gps.tsv", "gwas/pid_licsa/gps.tsv", "gwas/pid_sulf-all/gps.tsv",
"gwas/pid_divert/gps.tsv", "gwas/pid_misc-blood/gps.tsv", "gwas/pid_sys-scl/gps.tsv",
"gwas/pid_drug-all/gps.tsv", "gwas/pid_misc-immune/gps.tsv", "gwas/pid_t1d/gps.tsv",
"gwas/pid_dyschr-vit/gps.tsv", "gwas/pid_ms/gps.tsv", 	 "gwas/pid_t2d/gps.tsv",
"gwas/pid_endocard/gps.tsv", "gwas/pid_myal/gps.tsv", "gwas/pid_thy/gps.tsv",
"gwas/pid_endomet/gps.tsv", "gwas/pid_myas/gps.tsv", "gwas/pid_uc/gps.tsv",
"gwas/pid_ent-col/gps.tsv", "gwas/pid_pad/gps.tsv",
"gwas/pid_eos-eso/gps.tsv", "gwas/pid_paed-all/gps.tsv"]

"""
rule maf:
     input:
      "1000g/euro/qc/{chr}_qc.bed",
      "1000g/euro/qc/{chr}_qc.bim",
      "1000g/euro/qc/{chr}_qc.fam"
     output:
      "1000g/euro/qc/maf/{chr}.frq"
     shell:
      "plink --memory 16000 --threads 8 --bfile 1000g/euro/qc/{wildcards.chr}_qc --freq --out 1000g/euro/qc/maf/{wildcards.chr}"

rule ld:
     input:
      "1000g/euro/qc/{chr}_qc.bed",
      "1000g/euro/qc/{chr}_qc.bim",
      "1000g/euro/qc/{chr}_qc.fam"
     output:
      "1000g/euro/qc/ld/{chr}.ld"
     threads: 8
     shell:
      "plink --memory 16000 --threads 8 --bfile 1000g/euro/qc/{wildcards.chr}_qc --r2 --ld-window-r2 0.2 --ld-window-kb 1000 --out 1000g/euro/qc/ld/{wildcards.chr}"

rule prune_whole_set:
     input:
      "1000g/euro/qc/{chr}_qc.bed",
      "1000g/euro/qc/{chr}_qc.bim",
      "1000g/euro/qc/{chr}_qc.fam"
     output:
      "1000g/euro/qc/prune_whole/{chr}.prune.in",
      "1000g/euro/qc/prune_whole/{chr}.prune.out",
     threads: 8
     shell:
      "plink --memory 16000 --threads 8 --bfile 1000g/euro/qc/{wildcards.chr}_qc --indep-pairwise 1000kb 50 0.2 --out 1000g/euro/qc/prune_whole/{wildcards.chr}"

rule all:
     input:
      ("gwas/pid_{imd}/plink/pruned/chr%d.bed" % x for x in range(1,23)),
      ("gwas/pid_{imd}/plink/pruned/chr%d.bim" % x for x in range(1,23)),
      ("gwas/pid_{imd}/plink/pruned/chr%d.fam" % x for x in range(1,23))
     output:
      "gwas/pid_{imd}/done.txt"
     shell:
      "touch {output}"

rule prune_subset:
     input:
      "gwas/pid_{imd}/plink/{chr}.bed",
      "gwas/pid_{imd}/plink/{chr}.bim",
      "gwas/pid_{imd}/plink/{chr}.fam",
      range_file = "gwas/pid_{imd}/prune/{chr}.prune.in"
     output:
      "gwas/pid_{imd}/plink/pruned/{chr}.bed",
      "gwas/pid_{imd}/plink/pruned/{chr}.bim",
      "gwas/pid_{imd}/plink/pruned/{chr}.fam",
     params:
      bfile = "gwas/pid_{imd}/plink/{chr}",
      out = "gwas/pid_{imd}/plink/pruned/{chr}"
     threads: 8
     shell:
      "plink --memory 16000 --threads 8 --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule subset_ld:
     input:
      ("1000g/euro/qc/ld/chr%d.ld" % x for x in range(1,23)),
      "1000g/euro/qc/ld/chrX.ld",
      "gwas/pid_{imd}/pid_{imd}.tsv.gz"
     output:
      ("gwas/pid_{imd}/ld/chr%d.ld" % x for x in range(1,23)),
      "gwas/pid_{imd}/ld/chrX.ld"
     shell:
      "Rscript scripts/subset_ld_with_merged_gwas.R -i gwas/pid_{wildcards.imd}  -ld 1000g/euro/qc/ld -o gwas/pid_{wildcards.imd}/ld"

rule prune_all:
     input:
      ("1000g/euro/qc/prune_whole/chr%d.prune.in" % x for x in range(1,23)),
      "1000g/euro/qc/prune_whole/chrX.prune.in"

rule ld_all:
     input:
      ("1000g/euro/qc/ld/chr%d.ld" % x for x in range(1,23)),
      "1000g/euro/qc/ld/chrX.ld"
"""