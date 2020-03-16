efakes
=========================

Script to obtain the electron to photon misidentification fake factors (FF).

The code runs in 2 stages. The first one is a loop in the ntuple, to obtain all the information needed to get the FFs. The second is the calculation of the FFs.

## Loop

Creating rootfile with all information in outputDir/tag directory (output_loop.root). Also saving the used configuration (used_config.yaml)

	python main.py --loop --inputFiles /path1/sample1.root /path2/samples/sample*.root --tag sample1_sample2

Optional flags:

	--nentries 20000 (to run only this number of events)
	--methods TT1 TT2 (to run only some methods, default = TP TT1 TT2)
	--config /path/to/config.yaml (default = ./config.yaml)
	--outputDir /output/dir (default = ./output)


## FFs

Using the used_config.yaml previously calculated, a YAML file is created for each method, with the FF in each bin. The fit configuration is read from configFit, and is also saved (used_fit_config.yaml)

	python main.py --getFF --config /path/to/used_config.yaml


Optional flags:

	--methods TT1 TT2 (to run only some methods, default = TP TT1 TT2)
	--outputDirPlots /dir/to/save/plots (default = /eos/user/g/goorella/plots/efakes/  (!!!))
	--configFit /path/to/fit_config.yaml (default = ./fit_config.yaml)


## Systematics:

Alternative loops are needed to obtain some systematics (you need first to run the nominal loop). Everything is read from the  used_config.yaml

Runs up and down variation of the photon energy, with same config, only for TT2. Alpha is considered as percentage, eg: if alpha=0.05, the energy variations will be E\*1.05 and E\*0.95

	python main.py --loop --systEnergy --config /path/to/used_config.yaml --alpha 0.05


Runs the mass window variation, with same config, only for TT1

	python main.py --loop --systMassWin --config /path/to/used_config.yaml


Then run:
	
	python main.py --getFF --config /path/to/used_config.yaml --systEnergy
	python main.py --getFF --config /path/to/used_config.yaml --systMassWin

This will modify the previous YAML files including now this systematic.