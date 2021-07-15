rule download_1000g:
     output:
      "1000g/{chr}.vcf.gz"
     run:
      if wildcards.chr == 'chrX':
         shell("wget -O 1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz")
      elif wildcards.chr == 'chrY':
           shell("wget -O 1000g/{wildcards.chr}.vcf.gz      http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz")
      else:
        shell("wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

rule vcf_to_bed:
     input:
      "1000g/{chr}.vcf.gz"
     output:
      "1000g/{chr}.bed",
      "1000g/{chr}.bim",
      "1000g/{chr}.fam"
     shell:
      "plink --vcf {wildcards.chr}.vcf.gz --make-bed --out 1000g/{wildcards.chr}"

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
      "plink --bfile 1000g/{wildcards.chr} --keep 1000g/euro.fam --make-bed --silent --out 1000g/euro/{wildcards.chr}_euro"

