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
    params:
        out = "resources/1000g/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --vcf {input} --make-bed --out {params.out} --set-all-var-ids @:#:\$r:\$a --max-alleles 2 --new-id-max-allele-len 20 'truncate'"

rule make_euro_fam:
     input:
        "resources/1000g/20130606_g1k_3202_samples_ped_population.txt"
     output:
        "resources/1000g/euro.fam"
     shell:
        "Rscript workflow/scripts/get_euro_fam_file.R --ped_file {input} --output_file {output}"

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
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --keep {input.fam_file} --make-bed --silent --out {params.out}"

rule qc:
     input:
        ancient("resources/1000g/euro/{chr}_euro.bed"),
        ancient("resources/1000g/euro/{chr}_euro.bim"),
        ancient("resources/1000g/euro/{chr}_euro.fam")
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
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --make-bed --silent --out {params.out}"
