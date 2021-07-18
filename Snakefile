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
     benchmark:
      "benchmarks/{chr}_ld_benchmark.txt"
     threads: 8
     shell:
      "plink --threads 8 --bfile 1000g/euro/qc/{wildcards.chr}_qc --r2 --ld-window-r2 0.2 --ld-window-kb 1000 --out 1000g/euro/qc/ld/{wildcards.chr}"

rule all:
     input:
      ("1000g/euro/qc/ld/chr%d.ld" % x for x in range(1,23)),
      "1000g/euro/qc/ld/chrX.ld"