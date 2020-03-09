efakes
=========================

Script to obtain the electron to photon misidentification fake factors (FF)

## Loop

Config read from config.yaml

Creating rootfile with all information in 'output/tag' directory, e.g.:

	python main.py --input_loop /path/to/sample_1 /another/path/to/sample_2 /another/path/with/wildcards/sample_*_[3,4,5] --tag sample1_sample2_sample3

Optional flags:

	--nentries 20000 (to run only this number of events)
	--methods TT1 TT2 (to run only some methods, default = TP TT1 TT2)


## FFs

Config read from fit_config.yaml

	python main.py --input_FF /path/to/output/dir/


Optional flags:

	--methods TT1 TT2 (to run only some methods, default = TP TT1 TT2)
	--outputplots /dir/to/save/plots



## Systematic: reconstructed energy of the photon

Alternative loops are needed to obtain this systematic (you need first to run the nominal loop (to use the same config.yaml), and obtain FFs (to get the alpha from TT2.yaml)).

Runs up and down variation, with same config, only for TT2

	python main.py --syst_energy /path/to/output/dir/




