#! /usr/bin/env python


import os
import argparse
import glob
import math
from array import array 

import sys
sys.path.append('/afs/cern.ch/user/g/goorella/projects/harry_plotter')
import harry_plotter as hp

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *





class FitClass:

	YRange = None
	XRange = None

	extra = None

	def init(self):

		if self.sname == 'DSCB+gauss':

			self.par_sig_n = 10
			self.par_sig_1 = 7
			self.par_sig_2 = 3

		if self.sname == 'DSCB':

			self.par_sig_n = 7
			self.par_sig_1 = 7

		if self.bname == 'exp_poly':

			self.par_bkg_n = 3

		if self.bname == 'poly':

			self.par_bkg_n = 3

		if self.bname == 'gauss':

			self.par_bkg_n = 3


	def f_signal(self, x, par):

		l_par = []
		for i in xrange(len(par)):
			l_par.append(par[i])

		if self.sname == 'DSCB':
			return f_DoubleSidedCrystalBall(x,l_par[:self.par_sig_1])
		elif self.sname == 'DSCB+gauss':
			return f_DoubleSidedCrystalBall(x,l_par[:self.par_sig_1]) + f_gauss(x,l_par[self.par_sig_1:self.par_sig_1+self.par_sig_2])


	def f_background(self, x, par):

		if self.bname == 'exp_poly':
			return f_exp_poly(x, par)
		elif self.bname == 'poly':
			return f_poly(x, par)
		elif self.bname == 'gauss':
			return f_gauss(x, par)


	def f_total(self, x, par):

		l_par = []
		for i in xrange(len(par)):
			l_par.append(par[i])

		return self.f_signal(x, l_par[:self.par_sig_n]) + self.f_background(x, l_par[self.par_sig_n:self.par_sig_n+self.par_bkg_n])



	def generate_fit(self):

		self.init()

		ZCanvas = ROOT.TCanvas('tmp_canvas')

		pycallable = self.f_total
		fitFcn = ROOT.TF1('fitFcn', pycallable, self.rng[0], self.rng[1], self.par_sig_n+self.par_bkg_n)


		name = self.hist.GetTitle()
		if not name in self.config:
			print_msj('Using default fit config: '+name, 0)
			name = 'default'
			
		for p in self.config[name]:

			par = self.config[name][p]

			fitFcn.SetParName(par[0], p)
			fitFcn.SetParameter(par[0], par[1])
			fitFcn.SetParLimits(par[0], par[2], par[3])

		ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(10000)

		fitRes = self.hist.Fit(fitFcn, self.fit_opt, 'ep')

		params = fitFcn.GetParameters()
		signalFcn = ROOT.TF1('signalFcn', self.f_signal, self.rng[0], self.rng[1], self.par_sig_n)
		for i in xrange(self.par_sig_n):
			signalFcn.SetParameter(i, params[i])
		backFcn = ROOT.TF1('backFcn', self.f_background, self.rng[0], self.rng[1], self.par_bkg_n)
		for i in xrange(self.par_bkg_n):
			backFcn.SetParameter(i, params[i+self.par_sig_n])



		## plot


		os.system('mkdir -p %s' % self.output)
		pconfig = hp.PlotConfig()
		
		pconfig.Logy = 1

		pconfig.Output = self.output + self.hist.GetTitle() + '.pdf'

		pconfig.OptStat = 0
		pconfig.OptTitle = 0
		pconfig.RightMargin = 0.05
		# pconfig.LeftMargin = 0.11
		pconfig.TopMargin = 0.05
		# pconfig.BottomMargin = 0.12

		pconfig.XTitle = 'm_{ee} [GeV]'
		pconfig.XLabelSize = 0.05
		pconfig.XTitleSize = 0.05
		pconfig.XRange = (self.rng[0] , self.rng[1])
		if self.XRange is not None:
			pconfig.XRange = (self.XRange[0] , self.XRange[1])
		pconfig.XTitleOffset = 0.8
		pconfig.YTitle = 'Events'
		pconfig.YLabelSize = 0.05
		pconfig.YTitleSize = 0.05
		if self.YRange is not None:
			pconfig.YRange = (self.YRange[0], self.YRange[1])
		pconfig.YTitleOffset = 0.95

		pconfig.LegendX1 = 0.69
		pconfig.LegendX2 = 0.94
		pconfig.LegendY1 = 0.72
		pconfig.LegendY2 = 0.91
		pconfig.LegendTextSize = 0.035

		pconfig.LabelData = 'Data, #sqrt{s} = 13 TeV, 139.0 fb^{-1}'
		pconfig.LabelX = 0.32
		pconfig.LabelY = 0.8


		pair_tpye = 'ee'
		if '_eg_' in self.hist.GetTitle() :
			pconfig.XTitle = 'm_{e\gamma} [GeV]'
			pair_tpye = 'e\gamma'
		tmp_char = ['all', 'EE', 'BE', 'BB']
		if not any(x in self.hist.GetTitle() for x in tmp_char):
			bin_eta = int(self.hist.GetTitle().split('.')[0].split('_')[4])
			bin_pt = int(self.hist.GetTitle().split('.')[0].split('_')[5])
			pconfig.LabelData2 = '%s  |#eta|: %1.2f-%1.2f  p_{T}: %1.f-%1.f GeV' % (pair_tpye, self.extra[0][bin_eta], self.extra[0][bin_eta+1], self.extra[1][bin_pt], self.extra[1][bin_pt+1])


		l_po = []

		po = hp.PlotObject()

		po.Object = self.hist
		po.LineColor = ROOT.kBlack
		po.MarkerColor = ROOT.kBlack
		po.Legend = 'Data'
		po.LegendFill = 'lpe'
		po.Draw = 'same ep'
		po.LineWidth = 2 
		po.MarkerStyle = 20
		po.MarkerSize = 1.2
		l_po.append(po)


		po = hp.PlotObject()
		backFcn.SetNpx(10000)
		po.Object = backFcn
		po.LineColor = ROOT.kSpring+9
		po.Legend = 'Background'
		po.LegendFill = 'l'
		po.Draw = 'same'
		l_po.append(po)



		po = hp.PlotObject()
		signalFcn.SetNpx(10000)
		po.Object = signalFcn
		po.LineColor = ROOT.kPink
		po.Legend = 'Signal'
		po.LegendFill = 'l'
		po.Draw = 'same'
		l_po.append(po)



		po = hp.PlotObject()
		fitFcn.SetNpx(10000)
		po.Object = fitFcn
		po.LineColor = ROOT.kAzure-3
		po.Legend = 'Signal+Bkg'
		po.LegendFill = 'l'
		po.Draw = 'same'
		l_po.append(po)

		po = hp.PlotObject()
		fitFcn_cpy = fitFcn.Clone()
		fitFcn_cpy.SetNpx(10000)
		p0 = fitFcn.GetParameter(0)
		p1 = fitFcn.GetParameter(1)
		fitFcn_cpy.SetRange(p0-3*p1, p0+3*p1)
		po.Object = fitFcn_cpy
		po.LineColor = ROOT.kAzure-3
		po.FillColorAlpha = (ROOT.kAzure-3, 0.4)
		po.FillStyle = 1001
		po.Draw = 'fc same'
		l_po.append(po)

		hp.plot_main(l_po, pconfig)


		pconfig.Output = self.output + self.hist.GetTitle() + '_params.pdf'
		l_param = ''
		for i in xrange(len(self.config[name])):
			l_param += '#splitline{%s: %f}{' % ( fitFcn.GetParName(i), fitFcn.GetParameter(i) )
		l_param += ' }'*len(self.config[name])
		pconfig.LabelCustom2 = l_param
		pconfig.LabelCustom2X = 0.7
		pconfig.LabelCustom2Y = 0.62

		hp.plot_main(l_po, pconfig)

		return (signalFcn, backFcn, fitFcn, fitRes)







def f_exp_poly(x, par):

	# 3 parameters

	return math.exp(par[0]+par[1]*(x[0]-par[2])**2)


def f_poly(x, par, n_pam = 3):

	# 'n_pam' number of parameters

	y = 0.

	for i in xrange(n_pam):

		y += par[i]*(x[0]**i)

	return y

def f_gauss(x, par):

	# 3 parameters

	t = (x[0]-par[1])/par[2]

	return par[0]*math.exp(-0.5*t*t)


def f_DoubleSidedCrystalBall(x, par):

	# 7 parameters

	m0 = par[0]
	sigma = par[1]
	alphaLo = par[2]
	alphaHi = par[3]
	nLo = par[4]
	nHi = par[5]
	N = par[6]

	t = (x[0]-m0)/sigma

	if t < -alphaLo:
		a = math.exp(-0.5*alphaLo*alphaLo)
		b = nLo/alphaLo - alphaLo 
		return N*a/ROOT.TMath.Power(alphaLo/nLo*(b - t), nLo)
	
	elif t > alphaHi:
		a = math.exp(-0.5*alphaHi*alphaHi)
		b = nHi/alphaHi - alphaHi 
		return N*a/ROOT.TMath.Power(alphaHi/nHi*(b + t), nHi)

	return N*math.exp(-0.5*t*t)



def print_msj(msj, level):

	color = ['[1;93mWARNING', '[1;91mERROR']

	print '\033%s\033[0m: %s' % (color[level], msj)