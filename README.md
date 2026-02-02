
# Automated  Model Diagnostic Tool

Links:
- [Google planning doc](https://docs.google.com/document/d/1cy7wuHZThLaz3uDDgJ35CUjHeqQdUK2HIDssJzwClz4/edit?usp=sharing)


## TODO: Bartosz

- [x] check extract_fit_features
- [x] add error_scale, iiv_scale, frac_dense inside create_sim_data
- [ ] test few cases
- [ ] add snakemake pipeline
- [ ] time since last dose
- [ ] how to fix number of iterations

- [ ] Set a variable with common initial estimates for all the models tested
- [ ] Check where the seed should be so everything is reproducible in this part

## TODO/Comments

* Check what are realistic parameters and that we are we considering everything in: <https://pmc.ncbi.nlm.nih.gov/articles/PMC7894400/?utm>
* Check if it makes sense to scale the variance of all parameters or just some
* If we want other diagnostic metrics we could use ggpmx
* Make it so that some subjects did not take some doses 
* add missings doses



