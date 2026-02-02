
rule all:
	input:
		expand(
			"Output/fit_{mod_string}_i{iiv_scale}_e{error_scale}_f{frac_dense}.rds", 
			mod_string=["correct", "struct_miss", "resid_miss"], 
			iiv_scale=[0.5, 1, 2, 4],
			error_scale=[0.5, 1, 2, 4],
			frac_dense=[0.25, 0.5, 0.75, 1.0]
		)
	
rule gen_fit_obj:
	output:
		"Output/fit_{mod_string}_i{iiv_scale}_e{error_scale}_f{frac_dense}.rds"
	conda:
		"Envs/nlmixr2.yaml"
	script:
		"Scripts/gen_fit_obj.R"

