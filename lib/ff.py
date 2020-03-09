#! /usr/bin/env python


import os
import argparse
import glob
import math
import yaml
from array import array 


import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *


import lib.fitter as ft


def read_output(outputfile, cfg):

	h_grid = {}
	h_mass = {}

	file = TFile(outputfile)

	for t in ['ee', 'eg']:

		h_grid[t] = {}

		h_grid[t]['TP'] = file.Get('h_N_TP_'+t)
		h_grid[t]['TP'].SetDirectory(0)
		h_grid[t]['TT1_EE'] = file.Get('h_N_TT1_EE_'+t)
		h_grid[t]['TT1_EE'].SetDirectory(0)
		h_grid[t]['TT1_BE'] = file.Get('h_N_TT1_BE_'+t)
		h_grid[t]['TT1_BE'].SetDirectory(0)
		h_grid[t]['TT1_BB'] = file.Get('h_N_TT1_BB_'+t)
		h_grid[t]['TT1_BB'].SetDirectory(0)
		h_grid[t]['TT2'] = file.Get('h_N_TT2_'+t)
		h_grid[t]['TT2'].SetDirectory(0)

		h_mass['TT1_'+t] = {}
		for m in cfg['methods']:

			for r in ['EE', 'BE', 'BB']:

				h_mass['TT1_'+t][r] = file.Get('h_m_TT1_%s_%s'%(t, r))
				h_mass['TT1_'+t][r].SetDirectory(0)

		h_mass['TP_'+t] = []
		h_mass['TT2_'+t] = []
		for eta in xrange(len(cfg['eta_binning'])-1):

			h_mass['TP_'+t].append([])
			h_mass['TT2_'+t].append([])

			for pt in xrange(len(cfg['pt_binning'])-1):

				h_tmp = file.Get('h_m_TP_'+t+'_%i_%i'%(eta, pt))
				h_tmp.SetDirectory(0)
				h_mass['TP_'+t][eta].append(h_tmp)

				h_tmp = file.Get('h_m_TT2_'+t+'_%i_%i'%(eta, pt))
				h_tmp.SetDirectory(0)
				h_mass['TT2_'+t][eta].append(h_tmp)




	# h_mass['TT2_ee_all'] = file.Get('h_m_TT2_ee_all')
	# h_mass['TT2_ee_all'].SetDirectory(0)
	# h_mass['TT2_eg_all'] = file.Get('h_m_TT2_eg_all')
	# h_mass['TT2_eg_all'].SetDirectory(0)

	file.Close()

	return h_grid, h_mass



def get_FF(h_grid, h_mass, mthd, cfg):

	ff = {}


	with open('fit_config.yaml', 'r') as f:
		fit_config = yaml.safe_load(f)

	cl_fit = ft.FitClass()
	cl_fit.config = fit_config
	cl_fit.output = cfg['output_plots'] + 'fits/'
	cl_fit.fit_opt = 'SERQM0'

	if mthd == 'TT1':

		for t in ['ee', 'eg']:
			for reg in ['EE', 'BE', 'BB']:

				if h_mass['TT1_'+t][reg].GetEntries() == 0:
					print_msj('Empty histogram: h_m_TT1_'+t+'_'+reg, 0)
					ff['w_'+t+'_'+reg] = 0
					continue

				cl_fit.hist = h_mass['TT1_'+t][reg]
				cl_fit.rng = (70., 110.)
				cl_fit.sname = 'DSCB'
				cl_fit.bname = 'gauss'
				(sig, bkg, total, fitres) = cl_fit.generate_fit()
				p0 = total.GetParameter(0)
				p1 = total.GetParameter(1)

				ff['N_'+t+'_'+reg+'_total']     = total.Integral(p0-3*p1, p0+3*p1)
				ff['N_'+t+'_'+reg+'_total_err'] = total.IntegralError(p0-3*p1, p0+3*p1, fitres.GetParams(), fitres.GetCovarianceMatrix().GetMatrixArray())
				ff['N_'+t+'_'+reg+'_sig']     = sig.Integral(p0-3*p1, p0+3*p1)
				ff['N_'+t+'_'+reg+'_sig_err'] = sig.IntegralError(p0-3*p1, p0+3*p1, fitres.GetParams(), fitres.GetCovarianceMatrix().GetMatrixArray())

				ff['w_'+t+'_'+reg] = ff['N_'+t+'_'+reg+'_sig']/ff['N_'+t+'_'+reg+'_total']
				ff['w_'+t+'_'+reg+'_err'] = ff['w_'+t+'_'+reg]*math.sqrt( (ff['N_'+t+'_'+reg+'_sig_err']/ff['N_'+t+'_'+reg+'_sig'])**2 + (ff['N_'+t+'_'+reg+'_total_err']/ff['N_'+t+'_'+reg+'_total'])**2 )


		i_tmp = 0
		for eta in xrange(1, len(cfg['eta_binning'])):
			for pt in xrange(1, len(cfg['pt_binning'])):

				i_tmp+=1

				ff[i_tmp] = {}

				ff[i_tmp]['eta'] = '[%1.2f - %1.2f]' % (cfg['eta_binning'][eta-1], cfg['eta_binning'][eta])
				ff[i_tmp]['pT']  = '[%1.f - %1.f] GeV' % (cfg['pt_binning'][pt-1], cfg['pt_binning'][pt])

				for t in ['ee', 'eg']:
					ff[i_tmp]['N_'+t+'_all'] = 0.
					ff[i_tmp]['N_'+t+'_all_w1'] = 0.
					for reg in ['EE', 'BE', 'BB']:
						ff[i_tmp]['N_'+t+'_'+reg] = h_grid[t]['TT1_'+reg].GetBinContent(eta, pt)
						ff[i_tmp]['N_'+t+'_all'] += ff[i_tmp]['N_'+t+'_'+reg]*ff['w_'+t+'_'+reg]
						ff[i_tmp]['N_'+t+'_all_w1'] += ff[i_tmp]['N_'+t+'_'+reg]


				ff[i_tmp]['FF']    = safeDiv(ff[i_tmp]['N_eg_all'], ff[i_tmp]['N_ee_all'])
				ff[i_tmp]['FF_err'] = ff[i_tmp]['FF']*math.sqrt( (safeDiv(1, ff[i_tmp]['N_eg_all']))**2 + (safeDiv(1, ff[i_tmp]['N_ee_all']))**2 )
				ff[i_tmp]['FF_w1'] = safeDiv(ff[i_tmp]['N_eg_all_w1'], ff[i_tmp]['N_ee_all_w1'])
				ff[i_tmp]['FF_w1_err'] = ff[i_tmp]['FF_w1']*math.sqrt( (safeDiv(1, ff[i_tmp]['N_eg_all_w1']))**2 + (safeDiv(1, ff[i_tmp]['N_ee_all_w1']))**2 )


		f = open(cfg['output']+'TT1.yaml', 'w+')
		yaml.dump(ff, f, default_flow_style=False)
		f.close()
		print (cfg['output']+'TT1.yaml was created')


	else:

		# cl_fit.hist = h_mass['TT2_ee'][1][1]
		# cl_fit.rng = (60., 300.)
		# cl_fit.sname = 'DSCB+gauss'
		# cl_fit.bname = 'gauss'
		# (sig, bkg, total, fitres) = cl_fit.generate_fit()


		# cl_fit.rng = (60., 300.)
		# cl_fit.sname = 'DSCB'
		# cl_fit.bname = 'gauss'
		# cl_fit.hist = h_mass['TT2_ee_all']
		# (sig, bkg, total, fitres) = cl_fit.generate_fit()
		# p0 = total.GetParameter(0)
		# ff['m0_ee_all'] = p0
		# cl_fit.hist = h_mass['TT2_eg_all']
		# (sig, bkg, total, fitres) = cl_fit.generate_fit()
		# p0 = total.GetParameter(0)
		# ff['m0_eg_all'] = p0
		# ff['alpha_fake'] = ff['m0_eg_all']/ff['m0_ee_all']


		i_tmp = 0
		for eta in xrange(len(cfg['eta_binning'])-1):
			for pt in xrange(len(cfg['pt_binning'])-1):

				i_tmp+=1

				if h_mass[mthd+'_ee'][eta][pt].GetEntries() == 0 or h_mass[mthd+'_eg'][eta][pt].GetEntries() == 0:
					k_tmp = 'eta=[%1.2f-%1.2f]-pt=[%1.f-%1.f]' % (cfg['eta_binning'][eta], cfg['eta_binning'][eta+1], cfg['pt_binning'][pt], cfg['pt_binning'][pt+1])
					print_msj('Empty histograms: '+k_tmp, 0)
					continue

				ff[i_tmp] = {}

				ff[i_tmp]['eta'] = '[%1.2f - %1.2f]' % (cfg['eta_binning'][eta], cfg['eta_binning'][eta+1])
				ff[i_tmp]['pT']  = '[%1.f - %1.f] GeV' % (cfg['pt_binning'][pt], cfg['pt_binning'][pt+1])

				s_tmp = None
				for t in ['ee', 'eg']:

					cl_fit.hist = h_mass[mthd+'_'+t][eta][pt]
					cl_fit.rng = (60., 300.)
					cl_fit.sname = 'DSCB+gauss'
					cl_fit.bname = 'gauss'
					(sig, bkg, total, fitres) = cl_fit.generate_fit()

					p0 = total.GetParameter(0)
					p1 = total.GetParameter(1)
					ff[i_tmp]['m0_'+t] = p0
					ff[i_tmp]['s_'+t] = p1
					ff[i_tmp]['max_'+t+'_func'] = total.GetMaximum()
					ff[i_tmp]['max_'+t+'_hist'] = h_mass[mthd+'_'+t][eta][pt].GetMaximum()
					if t == 'ee':
						s_tmp = p1 # Nasella Thesis uses same sigma for ee and eg
					ff[i_tmp]['N_'+t+'_hist'] = h_mass[mthd+'_'+t][eta][pt].Integral(h_mass[mthd+'_'+t][eta][pt].FindFixBin(p0-3*s_tmp), h_mass[mthd+'_'+t][eta][pt].FindFixBin(p0+3*s_tmp), 'width')
					ff[i_tmp]['N_'+t+'_nobkg'] = total.Integral(p0-3*s_tmp, p0+3*s_tmp)
					ff[i_tmp]['N_'+t+'_sig']   = sig.Integral(p0-3*s_tmp, p0+3*s_tmp)
					ff[i_tmp]['N_'+t+'_sig_err'] = sig.IntegralError(p0-3*s_tmp, p0+3*s_tmp, fitres.GetParams(), fitres.GetCovarianceMatrix().GetMatrixArray())
					ff[i_tmp]['N_'+t+'_sig_w2']   = sig.Integral(p0-2*s_tmp, p0+2*s_tmp)
					ff[i_tmp]['N_'+t+'_sig_w4']   = sig.Integral(p0-4*s_tmp, p0+4*s_tmp)


				ff[i_tmp]['FF_hist']   = safeDiv(ff[i_tmp]['N_eg_hist'], ff[i_tmp]['N_ee_hist'])
				ff[i_tmp]['FF_nobkg']  = safeDiv(ff[i_tmp]['N_eg_nobkg'], ff[i_tmp]['N_ee_nobkg'])
				ff[i_tmp]['FF_w2']   = safeDiv(ff[i_tmp]['N_eg_sig_w2'], ff[i_tmp]['N_ee_sig_w2'])
				ff[i_tmp]['FF_w4']   = safeDiv(ff[i_tmp]['N_eg_sig_w4'], ff[i_tmp]['N_ee_sig_w4'])
				ff[i_tmp]['FF'] = safeDiv(ff[i_tmp]['N_eg_sig'], ff[i_tmp]['N_ee_sig'])
				ff[i_tmp]['FF_err_stat'] = ff[i_tmp]['FF']*math.sqrt( (ff[i_tmp]['N_eg_sig_err']/ff[i_tmp]['N_eg_sig'])**2 + (ff[i_tmp]['N_ee_sig_err']/ff[i_tmp]['N_ee_sig'])**2 )
				ff[i_tmp]['FF_err_syst_win']   = max([abs(ff[i_tmp]['FF']-ff[i_tmp]['FF_w2']), abs(ff[i_tmp]['FF']-ff[i_tmp]['FF_w2'])])
				ff[i_tmp]['FF_err_syst_nobkg'] = abs(ff[i_tmp]['FF']-ff[i_tmp]['FF_nobkg'])



		f = open(cfg['output']+mthd+'.yaml', 'w+')
		yaml.dump(ff, f, default_flow_style=False)
		f.close()
		print (cfg['output']+mthd+'.yaml was created')


	with open(cfg['output']+'used_fit_config.yaml', 'w+') as f:
		data = yaml.dump(cfg, f)
	print (cfg['output']+'used_fit_config.yaml was created')




def safeDiv(x, y):

	if y == 0.:
		return 0.
	return x/y




def print_msj(msj, level):

	color = ['[1;93mWARNING', '[1;91mERROR']

	print '\033%s\033[0m: %s' % (color[level], msj)