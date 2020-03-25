#! /usr/bin/env python


import os
import argparse
import glob
import math
import yaml
from array import array 

import sys
sys.path.append('/afs/cern.ch/user/g/goorella/harry_plotter')

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

import harry_plotter as hp



sty_range = {
	'pT' : (75., 300.),
	'eta' : (0., 2.5),
}

sty_axis = {
	'pT' : 'p_{T} [GeV]',
	'eta' : '|#eta|',
}

colors = [861, 810, 798, 878, 844, 404, 433, 607, 880, 425, 840, 820, 801]

def read_FF(cfg, ff_res):

	l_ff_pt = []
	l_ff_eta = []

	for j in xrange(len(cfg['eta_binning'])-1):

		x = array('f')
		y = array('f')
		dx = array('f')
		dy = array('f')

		if ff_res[j*(len(cfg['pt_binning'])-1)+1]['eta'][0] >= 1.37 and ff_res[j*(len(cfg['pt_binning'])-1)+1]['eta'][0] < 1.52: 			continue

		last_i = None
		for i in xrange(j*(len(cfg['pt_binning'])-1)+1, (j+1)*(len(cfg['pt_binning'])-1)+1):

			if ff_res[i]['FF']<1E-07: continue

			dx_tmp = (ff_res[i]['pT'][1]-ff_res[i]['pT'][0])/2.
			dx.append(dx_tmp)
			x_tmp = ff_res[i]['pT'][0] + dx_tmp
			x.append(x_tmp)
			y_tmp = ff_res[i]['FF']
			y.append(y_tmp)
			dy_tmp = ff_res[i]['FF_err_total']
			dy.append(dy_tmp)

			last_i = i

		if len(x)>0:
			l_ff_pt.append( (ROOT.TGraphErrors(len(cfg['pt_binning'])-1, x, y, dx, dy), '#eta = [%1.2f-%1.2f]' % (ff_res[last_i]['eta'][0], ff_res[last_i]['eta'][1]) ) )


	for j in xrange(1, len(cfg['pt_binning'])):

		x = array('f')
		y = array('f')
		dx = array('f')
		dy = array('f')

		last_i = None
		for i in xrange(j, j+1+(len(cfg['eta_binning'])-2)*(len(cfg['pt_binning'])-1), len(cfg['pt_binning'])-1):

			if ff_res[i]['FF']<1E-07: continue

			dx_tmp = (ff_res[i]['eta'][1]-ff_res[i]['eta'][0])/2.
			dx.append(dx_tmp)
			x_tmp = ff_res[i]['eta'][0] + dx_tmp
			x.append(x_tmp)
			y_tmp = ff_res[i]['FF']
			y.append(y_tmp)
			dy_tmp = ff_res[i]['FF_err_total']
			dy.append(dy_tmp)

			last_i = i

		l_ff_eta.append( (ROOT.TGraphErrors(len(cfg['eta_binning'])-2, x, y, dx, dy), 'p_{T} = [%1.f-%1.f] GeV' % (ff_res[last_i]['pT'][0], ff_res[last_i]['pT'][1]) ) )


	return l_ff_pt, l_ff_eta




def plot_FF(l_ff, comp = None, **kwargs):

	pconfig = hp.PlotConfig()
	
	# pconfig.Logy = 1

	pconfig.Output = kwargs['output']

	pconfig.OptStat = 0
	pconfig.OptTitle = 0
	pconfig.RightMargin = 0.05
	pconfig.LeftMargin = 0.11
	pconfig.TopMargin = 0.06
	pconfig.BottomMargin = 0.12

	pconfig.XTitle = sty_axis[kwargs['var']]
	pconfig.XLabelSize = 0.05
	pconfig.XTitleSize = 0.05
	pconfig.XRangeUser = sty_range[kwargs['var']]

	pconfig.YTitle = 'FF'
	pconfig.YLabelSize = 0.05
	pconfig.YTitleSize = 0.05
	pconfig.YRangeUser = (0., 0.15)
	pconfig.YTitleOffset = 1.1

	pconfig.LegendX1 = 0.52
	pconfig.LegendX2 = 0.93
	pconfig.LegendY1 = 0.65
	pconfig.LegendY2 = 0.93
	pconfig.LegendTextSize = 0.033

	pconfig.LabelData = '#sqrt{s} = 13 TeV, 139.0 fb^{-1}'
	pconfig.LabelX = 0.2
	pconfig.LabelY = 0.8

	if comp is not None:
		pconfig.LegendNColumns = 2

	l_po = []

	for i,p in enumerate(l_ff):

		po = hp.PlotObject()
		po.Object = p[0]
		po.LineColor = colors[i]
		po.MarkerColor = colors[i]
		po.FillColor = colors[i]
		po.Legend = p[1]
		if comp is not None:
			po.Legend = kwargs['comp_label'][0]
		po.LegendFill = 'lp'
		po.Draw = 'EPZ same'
		po.LineWidth = 2 
		po.MarkerStyle = 20
		po.MarkerSize = 1.2

		l_po.append(po)

		if comp is not None:

			po = hp.PlotObject()
			po.Object = comp[i][0]
			po.LineColor = colors[i]
			po.MarkerColor = colors[i]
			po.FillColor = colors[i]
			po.Legend = kwargs['comp_label'][1] + '  ' + comp[i][1]
			po.LegendFill = 'lp'
			po.Draw = 'EPZ same'
			po.LineWidth = 2 
			po.LineStyle = 2
			po.MarkerStyle = 24
			po.MarkerSize = 1.2

			l_po.append(po)




	hp.plot_main(l_po, pconfig)


def latex_table(cfg, ff_res, output = 'latex_table.tex'):

	f = open(output, 'w+')

	f.write('\\begin{table} \n')
	f.write('\\begin{tabular}{ l l | c | c c c c | c } \n')
	f.write('\\hline \n')
	f.write('\\hline \n')
	f.write('\\multirow{2}{*}{$|\\eta|$} & \\multirow{2}{*}{$p_{T}$[GeV]} & \\multirow{2}{*}{Fake factor} & Statistical & \\multicolumn{3}{c |}{Systematics Unc} & Total \\\\ \n')
	f.write('\\cline{5-7} \n')
	f.write(' & & & Unc & Integrating win & No bkg sustr & Energy bias & Unc\\\\ \n')
	f.write('\\hline \n')
	f.write('\\hline \n')

	for j in xrange(len(cfg['eta_binning'])-1):

		if ff_res[j*(len(cfg['pt_binning'])-1)+1]['eta'][0] >= 1.37 and ff_res[j*(len(cfg['pt_binning'])-1)+1]['eta'][0] < 1.52: 			continue

		i_tmp = 0
		for i in xrange(j*(len(cfg['pt_binning'])-1)+1, (j+1)*(len(cfg['pt_binning'])-1)+1):

			line_tmp = ''
			if i_tmp == 0:
				line_tmp = '\\multirow{%i}{*}{%1.2f-%1.2f}' % (len(cfg['pt_binning'])-1, ff_res[i]['eta'][0], ff_res[i]['eta'][1])
			line = '%s & %1.f-%1.f & %1.4f & %1.4f (%1.1f\\%%) & %1.4f (%1.1f\\%%) & %1.4f (%1.1f\\%%) & %1.4f (%1.1f\\%%) & %1.4f (%1.1f\\%%) \\\\ \n' % (line_tmp, ff_res[i]['pT'][0], ff_res[i]['pT'][1], ff_res[i]['FF'], ff_res[i]['FF_err_stat_2'], 100.*ff_res[i]['FF_err_stat_2']/ff_res[i]['FF'], ff_res[i]['FF_err_syst_win'], 100.*ff_res[i]['FF_err_syst_win']/ff_res[i]['FF'], ff_res[i]['FF_err_syst_nobkg'], 100.*ff_res[i]['FF_err_syst_nobkg']/ff_res[i]['FF'], ff_res[i]['FF_err_syst_energy'], 100.*ff_res[i]['FF_err_syst_energy']/ff_res[i]['FF'], ff_res[i]['FF_err_total'], 100.*ff_res[i]['FF_err_total']/ff_res[i]['FF'])
			f.write(line)

			i_tmp+=1

		f.write('\\hline \n')

	f.write('\\hline \n')
	f.write('\\end{tabular} \n')
	f.write('\\end{table} \n')

	print output+' has been created'
