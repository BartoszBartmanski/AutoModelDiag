
rule all:
	input:
		expand(
			"Output/features_{mod_string}_i{iiv_scale}_e{error_scale}_f{frac_dense}.json", 
			mod_string=["correct", "struct_miss", "resid_miss"], 
			iiv_scale=[0.5, 1, 2, 4],  # 0.3, 0.6, 0.9
			error_scale=[0.5, 1, 2, 4],
			frac_dense=[0.25, 0.5, 0.75, 1.0]
		)
	
rule gen_fit_obj:
	output:
		"Output/fit_{mod_string}_i{iiv_scale}_e{error_scale}_f{frac_dense}.rds"
	conda:
		"Envs/nlmixr2.yaml"
	resources:
		runtime = lambda wc, attempt: 120 * attempt,
		mem_mb = lambda wc, attempt: 2000 * attempt
	script:
		"Scripts/gen_fit_obj.R"

rule extract_fit_features:
	input:
		"Output/fit_{mod_string}_i{iiv_scale}_e{error_scale}_f{frac_dense}.rds"
	output:
		"Output/features_{mod_string}_i{iiv_scale}_e{error_scale}_f{frac_dense}.json"
	conda:
		"Envs/nlmixr2.yaml"
	resources:
		runtime = lambda wc, attempt: 20 * attempt
	script:
		"Scripts/gen_fit_obj.R"
