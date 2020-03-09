#! /usr/bin/env python

import os
import argparse
import yaml

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *


import lib.loop as lp
import lib.ff as ff


def main():

	if syst_energy is not None:

		with open(syst_energy+'used_config.yaml', 'r') as f:
			cfg = yaml.safe_load(f)

		with open(syst_energy+'TT2.yaml', 'r') as f:
			alph = yaml.safe_load(f)
		if 'alpha_fake' not in alph:
			print_msj('Alpha fake not found.', 1)
			return 1
		alpha_tmp = alph['alpha_fake']

		cfg['syst_energy'] = True
		cfg['methods'] = ['TT2']

		cfg['alpha'] = 1+abs(1-alpha_tmp)
		lp.init_loop(cfg)

		print ''

		cfg['alpha'] = 1-abs(1-alpha_tmp)
		lp.init_loop(cfg)

	elif len(input_loop)>0:

		os.system('mkdir -p output/%s' % tag)

		with open('config.yaml', 'r') as f:
			cfg = yaml.safe_load(f)

		cfg['path'] = input_loop
		cfg['nevents'] = nevents
		cfg['tag'] = tag
		cfg['syst_energy'] = False
		cfg['alpha'] = 1.
		if len(methods)>0:
			cfg['methods'] = methods
		lp.init_loop(cfg)

	elif input_FF is not None:

		with open(input_FF+'used_config.yaml', 'r') as f:
			cfg = yaml.safe_load(f)

		if len(methods)>0:
			cfg['methods'] = methods
			
		h_grid, h_mass = ff.read_output(input_FF+'output_loop.root', cfg)

		cfg['output'] = input_FF
		cfg['output_plots'] = outputplots + input_FF.split('/')[-2] + '/'

		for m in cfg['methods']:
			ff.get_FF(h_grid, h_mass, m, cfg)



def check_args(args):

	# improve this!

	return 0




parser = argparse.ArgumentParser()

parser.add_argument('--input_loop', dest='input_loop', default=[], nargs='+')
parser.add_argument('--input_FF', dest='input_FF', type=str, default=None)
parser.add_argument('--tag', dest='tag', type=str, default='test_tag')
parser.add_argument('--outputplots', dest='outputplots', type=str, default='/eos/user/g/goorella/plots/efakes/test_tag')
parser.add_argument('--nevents', dest='nevents', type=int, default=None)
parser.add_argument('--syst_energy', dest='syst_energy', type=str, default=None)
# parser.add_argument('--alpha', dest='alpha', type=str, default=None)
parser.add_argument('--methods', dest='methods', default=[], nargs='+')

args = parser.parse_args()
check_args(args)

input_loop = args.input_loop
input_FF = args.input_FF
tag = args.tag
outputplots = args.outputplots
nevents = args.nevents
syst_energy = args.syst_energy
# alpha = args.alpha
methods = args.methods


for i in input_loop:
	if i[-1] != '/':
		i += '/'
if input_FF is not None and input_FF[-1] != '/':
	input_FF += '/'
if outputplots[-1] != '/':
	outputplots += '/'
if syst_energy is not None and syst_energy[-1] != '/':
	syst_energy += '/'


if __name__ == '__main__':
	main()
#  




def print_msj(msj, level):

	color = ['[1;93mWARNING', '[1;91mERROR']

	print '\033%s\033[0m: %s' % (color[level], msj)