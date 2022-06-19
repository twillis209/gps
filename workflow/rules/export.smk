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

top_imds = ['uc-delange', 'sys-scl', 'lada', 't1d', 'psc', 'jia', 'addi', 'ms', 'igad', 'sle', 'pso', 'aster', 'pbc', 'ra', 'igm']

top_imd_pairs = list(chain(*[[f"{top_imds[i]}_{top_imds[j]}" for j in range(i+1,len(top_imds))] for i in range(len(top_imds))]))

rule estimate_rg_with_ldak_thin_for_top_imds:
    input:
        [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.cors" for trait_pair in top_imd_pairs]
#        [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_sans-mhc.cors" for trait_pair in top_imd_pairs]
