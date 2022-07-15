from itertools import chain

ex_traits = [
    "uc-delange",
    "t1d",
    "sle",
    "ra",
    "crohns",
    "hypothy",
    "derm-ecz",
    "rhin-ex",
    "asthma-ex",
    "cardiomyo-ex",
    "endomet-ex",
    "ibs-ex",
    "leio-ex",
    "cholelith-ex",
    "osteo-ex",
    "hyperchol-ex",
    "md-ex",
    "glaucoma-ex"
]

pan_ukbb_ex_traits = [
    "rhin-ex",
    "asthma-ex",
    "ibs-ex",
    "cholelith-ex",
    "osteo-ex",
    "hyperchol-ex",
    "md-ex"
]

existing_pairs = [x for x in top_imd_pairs if all([y in ex_traits for y in x.split('_')])]

ex_trait_pairs = list(chain(*[[f"{ex_traits[i]}_{ex_traits[j]}" for j in range(i+1,len(ex_traits))] for i in range(len(ex_traits))]))

rule create_pan_ukbb_meta_file:
    input:
        "resources/gwas/{trait}-ex.tsv.gz"
    output:
        "resources/gwas/{trait}-meta-ex.tsv.gz"
    params:
        columns_to_drop = ["af_cases_meta_hq", "af_controls_meta_hq", "pval_heterogeneity_hq", "af_cases_meta", "af_controls_meta", "beta_meta_hq", "se_meta_hq", "pval_meta_hq", "pval_heterogeneity", "af_cases_AFR", "af_cases_AMR", "af_cases_CSA", "af_cases_EAS", "af_cases_EUR", "af_cases_MID", "af_controls_AFR", "af_controls_AMR", "af_controls_CSA", "af_controls_EAS", "af_controls_EUR", "af_controls_MID", "beta_AFR", "beta_AMR", "beta_CSA", "beta_EAS", "beta_MID", "se_AFR", "se_AMR", "se_CSA", "se_EAS", "se_MID", "pval_AFR", "pval_AMR", "pval_CSA", "pval_EAS", "pval_MID", "low_confidence_AFR", "low_confidence_AMR", "low_confidence_CSA", "low_confidence_EAS", "low_confidence_EUR", "low_confidence_MID", "nearest_genes"],
        pval_col = "pval_meta",
        beta_col = "beta_meta",
        se_col = "se_meta"
    threads: 4
    script: "../scripts/create_pan_ukbb_meta_file.R"

rule estimate_rg_with_ldak_thin_for_gps_paper_traits:
    input:
        [y for y in [f"results/ldak/ldak-thin/{trait_pair}/{trait_pair}_all.cors" for trait_pair in ex_trait_pairs if trait_pair not in existing_pairs] if not os.path.exists(y)]

rule process_ex_datasets:
    input:
        [f"results/processed_gwas/{x}_recalculated_p.tsv.gz" for x in ex_traits]

rule process_pan_ukbb_datasets:
    input:
        [f"results/processed_gwas/{x}_recalculated_p.tsv.gz" for x in pan_ukbb_ex_traits]
