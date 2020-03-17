#! /usr/bin/env python

import os
import argparse
import yaml
import glob
import sys

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *


import lib.loop as lp
import lib.ff as ff


def main():

	with open(config, 'r') as f:
		cfg = yaml.safe_load(f)

	cfg['syst_energy'] = systEnergy
	cfg['syst_masswin'] = systMassWin

	if loop:


		if systEnergy:

			cfg['methods'] = ['TT2']
			cfg['alpha'] = 1+alpha

			print 'Running loop with:'
			print '\t Tag: %s' % cfg['tag']
			print '\t Output Dir: %s' % cfg['outputdir'] + cfg['tag']
			print '\t Number of events: %s' % cfg['nevents']
			print '\t Methods: %s' % ' '.join(cfg['methods'])
			print '\t Alpha: %1.2f' % cfg['alpha']
			print '\t Syts Energy: %s' % cfg['syst_energy']
			print '\t Syts Mass Win: %s' % cfg['syst_masswin']

			lp.init_loop(cfg)

			cfg['alpha'] = 1-alpha

			print '\n\t Alpha: %1.2f' % cfg['alpha']

			lp.init_loop(cfg)

		elif systMassWin:

			cfg['methods'] = ['TT1']
			cfg['TT1']['cut_mass'] = cfg['TT1']['cut_mass_up']

			print 'Running loop with:'
			print '\t Tag: %s' % cfg['tag']
			print '\t Output Dir: %s' % cfg['outputdir'] + cfg['tag']
			print '\t Number of events: %s' % cfg['nevents']
			print '\t Methods: %s' % ' '.join(cfg['methods'])
			print '\t Alpha: %1.2f' % cfg['alpha']
			print '\t Syts Energy: %s' % cfg['syst_energy']
			print '\t Syts Mass Win: %s' % cfg['syst_masswin']

			lp.init_loop(cfg)

		else:

			os.system('mkdir -p %s/%s' % (outputDir, tag))

			cfg['tag'] = tag
			cfg['path'] = inputFiles
			cfg['nevents'] = nevents
			cfg['alpha'] = 1.
			cfg['outputdir'] = outputDir
			cfg['methods'] = methods

			print 'Running loop with:'
			print '\t Tag: %s' % cfg['tag']
			print '\t Output Dir: %s' % cfg['outputdir'] + cfg['tag']
			print '\t Number of events: %s' % cfg['nevents']
			print '\t Methods: %s' % ' '.join(cfg['methods'])
			print '\t Alpha: %1.2f' % cfg['alpha']
			print '\t Syts Energy: %s' % cfg['syst_energy']
			print '\t Syts Mass Win: %s' % cfg['syst_masswin']

			lp.init_loop(cfg)



	elif getFF:

		cfg['fit_config'] = configFit

		if systEnergy:

			cfg['methods'] = 'TT2'
			file = glob.glob(cfg['outputdir'] + cfg['tag'] + '/output_loop_syst_energy_*.root')
			file.sort()

			h_grid, h_mass = ff.read_output(file[0], cfg)
			cfg['output_plots'] = outputDirPlots + cfg['tag'] + '_syst_energy_%s/' % file[0].replace('.root','').split('_')[-1]

			ff.get_FF(h_grid, h_mass, 'TT2', cfg)

			print ''

			h_grid, h_mass = ff.read_output(file[1], cfg)
			cfg['output_plots'] = outputDirPlots + cfg['tag'] + '_syst_energy_%s/' % file[1].replace('.root','').split('_')[-1]

			ff.get_FF(h_grid, h_mass, 'TT2', cfg)

		elif systMassWin:

			cfg['methods'] = 'TT1'

			h_grid, h_mass = ff.read_output(cfg['outputdir'] + cfg['tag'] +'/output_loop_syst_masswin.root', cfg)

			cfg['output_plots'] = outputDirPlots + cfg['tag'] + '_syst_masswin/'

			ff.get_FF(h_grid, h_mass, 'TT1', cfg)

		else:

			if len(methods)>0:
				cfg['methods'] = methods
				
			h_grid, h_mass = ff.read_output(cfg['outputdir'] + cfg['tag'] + '/output_loop.root', cfg)
			cfg['output_plots'] = outputDirPlots + cfg['tag'] + '/'

			for m in cfg['methods']:
				ff.get_FF(h_grid, h_mass, m, cfg)




def check_args(args):

	if args.loop and args.getFF:
		return 1
	if args.systEnergy and args.loop and args.alpha is None:
		return 1
	if args.loop and len(args.inputFiles) == 0:
		return 1
	if len(args.methods) > 0 and ('TP' not in args.methods or 'TT1' not in args.methods or 'TT2' not in args.methods):
		return 1


	return 0




parser = argparse.ArgumentParser()

parser.add_argument('--loop', dest='loop', action='store_true', default=False)
parser.add_argument('--inputFiles', dest='inputFiles', default=[], nargs='+')
parser.add_argument('--tag', dest='tag', type=str, default='default_tag')
parser.add_argument('--outputDir', dest='outputDir', type=str, default='output/')
parser.add_argument('--nevents', dest='nevents', type=int, default=None)

parser.add_argument('--getFF', dest='getFF', action='store_true', default=False)
parser.add_argument('--configFit', dest='configFit', type=str, default='fit_config.yaml')
parser.add_argument('--outputDirPlots', dest='outputDirPlots', type=str, default='/eos/user/g/goorella/plots/efakes/')

parser.add_argument('--config', dest='config', type=str, default='config.yaml')

parser.add_argument('--methods', dest='methods', default=['TP', 'TT1', 'TT2'], nargs='+')

parser.add_argument('--systEnergy', dest='systEnergy', action='store_true', default=False)
parser.add_argument('--alpha', dest='alpha', type=float, default=None)

parser.add_argument('--systMassWin', dest='systMassWin', action='store_true', default=False)

args = parser.parse_args()
check_args(args)

loop = args.loop
inputFiles = args.inputFiles

getFF = args.getFF
configFit = args.configFit
outputDirPlots = args.outputDirPlots

tag = args.tag

outputDir = args.outputDir

config = args.config

methods = args.methods

nevents = args.nevents

systEnergy = args.systEnergy
alpha = args.alpha

systMassWin = args.systMassWin


for i in inputFiles:
	if i[-1] != '/':
		i += '/'
if outputDirPlots[-1] != '/':
	outputDirPlots += '/'
if outputDir[-1] != '/':
	outputDir += '/'



if __name__ == '__main__':
	main()
#  




def print_msj(msj, level):

	color = ['[1;93mWARNING', '[1;91mERROR']

	print '\033%s\033[0m: %s' % (color[level], msj)