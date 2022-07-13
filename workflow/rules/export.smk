rule run_gps_on_ukbb_traits:
    input:
        [f"results/pid_{trait}/pid_{trait}_all_gps_value.tsv" for trait in ukbb_traits]

rule run_naive_gps_on_finngen_traits:
    input:
        [f"results/pid_{trait}/pid_{trait}_all_gps_value_naive.tsv" for trait in finngen_traits]

finngen_pairs = list(chain(*[[f"{finngen_traits[i]}_{finngen_traits[j]}" for j in range(i+1,len(finngen_traits))] for i in range(len(finngen_traits))]))

rule run_sumher_on_finngen_traits:
    input:
        [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.cors" for trait_pair in finngen_pairs]

top_imds = ['uc-delange', 'sys-scl', 'lada', 't1d', 'psc', 'jia', 'addi', 'ms', 'igad', 'sle', 'pso', 'pbc', 'ra', 'igm', 'asthma', 'sjogren', 'hypothy', 'hyperthy', 'polymyal', 'still', 'ank-spond', 'spondylo', 'licp']

top_imd_pairs = list(chain(*[[f"{top_imds[i]}_{top_imds[j]}" for j in range(i+1,len(top_imds))] for i in range(len(top_imds))]))

rule estimate_rg_with_ldak_thin_for_top_imds:
    input:
        [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.cors" for trait_pair in top_imd_pairs]

rule run_gps_on_finngen_traits:
    input:
        [f"results/pid_{trait}/pid_{trait}_all_gps_value.tsv" for trait in finngen_traits]

rule run_naive_gps_on_all_traits:
    input:
        [f"results/pid_{trait}/pid_{trait}_all_gps_value_naive.tsv" for trait in finngen_traits+ukbb_traits+misc_traits]
    output:
        "results/compiled_all_gps_value_naive.tsv"
    shell:
        """
        echo -e "Trait_A\tTrait_B\tGPS" > {output}
        for x in {input}; do
            tail -n +2 $x >> {output}
        done
        """

rule compile_top_imd_sumher_results:
    input:
       cors_files = [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.cors.full" for trait_pair in top_imd_pairs],
       overlap_files = [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.overlap" for trait_pair in top_imd_pairs],
    output:
        "results/ldak/ldak-thin/compiled_top_imd_results.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A\ttrait_B\ttag_predictors\tsumstat_predictors\toverlap\th2.A\th2.A.se\th2.B\th2.B.se\tgcov\tgcov.se\trg\trg.se\trg.z\trg.p\n")
            for i, x in enumerate(input.cors_files):
                head, tail = os.path.split(x)

                trait_A, trait_B = re.match("([\w-]+)_([\w-]+)_all.cors.full", tail).groups()

                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                    # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
                    _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()
                    rg_z = float(rg)/float(rg_se)

                    rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)

                with open(input.overlap_files[i], 'r') as overlap_infile:
                    tag_predictors = overlap_infile.readline().split()[1]
                    sumstat_predictors = overlap_infile.readline().split()[1]
                    overlap = overlap_infile.readline().split()[1]
                    outfile.write(f"{trait_A}\t{trait_B}\t{tag_predictors}\t{sumstat_predictors}\t{overlap}\t{float(h2_A)}\t{float(h2_A_se)}\t{float(h2_B)}\t{float(h2_B_se)}\t{float(cov)}\t{float(cov_se)}\t{float(rg)}\t{float(rg_se)}\t{rg_z}\t{rg_p}\n")
