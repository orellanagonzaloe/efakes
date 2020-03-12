#! /usr/bin/env python


import os
import argparse
import glob
import math
import yaml
import copy
from array import array 


import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *


cfg = {}


def init_loop(cfg_tmp):

	global cfg 
	cfg = copy.deepcopy(cfg_tmp)

	### Initialize histograms ###

	h_grid = {}
	h_mass = {}

	for t in ['ee', 'eg']:

		h_grid[t] = {}

		h_grid[t]['TP'] = ROOT.TH2D( 'h_N_TP_'+t, 'h_N_TP_'+t, len(cfg['eta_binning'])-1, array('f',cfg['eta_binning']), len(cfg['pt_binning'])-1, array('f',cfg['pt_binning']))
		h_grid[t]['TP'].Sumw2(ROOT.kTRUE)
		h_grid[t]['TT1_EE'] = ROOT.TH2D( 'h_N_TT1_EE_'+t, 'h_N_TT1_EE_'+t, len(cfg['eta_binning'])-1, array('f',cfg['eta_binning']), len(cfg['pt_binning'])-1, array('f',cfg['pt_binning']))
		h_grid[t]['TT1_EE'].Sumw2(ROOT.kTRUE)
		h_grid[t]['TT1_BE'] = ROOT.TH2D( 'h_N_TT1_BE_'+t, 'h_N_TT1_BE_'+t, len(cfg['eta_binning'])-1, array('f',cfg['eta_binning']), len(cfg['pt_binning'])-1, array('f',cfg['pt_binning']))
		h_grid[t]['TT1_BE'].Sumw2(ROOT.kTRUE)
		h_grid[t]['TT1_BB'] = ROOT.TH2D( 'h_N_TT1_BB_'+t, 'h_N_TT1_BB_'+t, len(cfg['eta_binning'])-1, array('f',cfg['eta_binning']), len(cfg['pt_binning'])-1, array('f',cfg['pt_binning']))
		h_grid[t]['TT1_BB'].Sumw2(ROOT.kTRUE)
		h_grid[t]['TT2'] = ROOT.TH2D( 'h_N_TT2_'+t, 'h_N_TT2_'+t, len(cfg['eta_binning'])-1, array('f',cfg['eta_binning']), len(cfg['pt_binning'])-1, array('f',cfg['pt_binning']))
		h_grid[t]['TT2'].Sumw2(ROOT.kTRUE)


		h_mass['TT1_'+t] = {}
		for r in ['EE', 'BE', 'BB']:
			h_mass['TT1_'+t][r] = ROOT.TH1D('h_m_TT1_%s_%s'%(t, r), 'h_m_TT1_%s_%s'%(t, r), 400, 0., 400.)
			h_mass['TT1_'+t][r].Sumw2(ROOT.kTRUE)

		h_mass['TP_'+t] = []
		h_mass['TT2_'+t] = []

		for eta in xrange(len(cfg['eta_binning'])-1):

			h_mass['TP_'+t].append([])
			h_mass['TT2_'+t].append([])

			for pt in xrange(len(cfg['pt_binning'])-1):

				h_tmp = ROOT.TH1D('h_m_TP_'+t+'_%i_%i'%(eta, pt), 'h_m_TP_'+t+'_%i_%i'%(eta, pt), 400, 0., 400.)
				h_tmp.Sumw2(ROOT.kTRUE)
				h_mass['TP_'+t][eta].append(h_tmp)
				h_tmp = ROOT.TH1D('h_m_TT2_'+t+'_%i_%i'%(eta, pt), 'h_m_TT2_'+t+'_%i_%i'%(eta, pt), 400, 0., 400.)
				h_tmp.Sumw2(ROOT.kTRUE)
				h_mass['TT2_'+t][eta].append(h_tmp)

	h_mass['TT2_ee_all'] = ROOT.TH1D('h_m_TT2_ee_all', 'h_m_TT2_ee_all', 400, 0., 400.)
	h_mass['TT2_ee_all'].Sumw2(ROOT.kTRUE)
	h_mass['TT2_eg_all'] = ROOT.TH1D('h_m_TT2_eg_all', 'h_m_TT2_eg_all', 400, 0., 400.)
	h_mass['TT2_eg_all'].Sumw2(ROOT.kTRUE)


	loop(h_grid, h_mass)



def loop(h_grid, h_mass):

	chain = TChain('mini')

	all_files = []

	for file in cfg['path']:
		chain.Add(file)

	nevents = cfg['nevents']
	if nevents is None:
		nevents = chain.GetEntries()
	print('Total events =  %i' % nevents)

	n_TP_ee = 0
	n_TP_eg = 0

	n_TT1_ee = 0
	n_TT1_eg = 0

	n_TT2_ee = 0
	n_TT2_eg = 0


	for entry in xrange(nevents):


		if entry % int(nevents/10) == 0:
			print ('%i%%') % round(100.*entry/nevents)

		chain.GetEntry(entry)

		el_n = chain.el_medium_n
		ph_n = chain.ph_noniso_n


		# small skim
		if el_n < 1 or (el_n < 2 and ph_n < 1):
			continue

		# momentary using 'medium' branches
		el_pt = chain.el_medium_pt
		el_etas2 = chain.el_medium_etas2
		el_phi = chain.el_medium_phi
		el_ch = chain.el_medium_ch

		# momentary using 'noniso' branches
		ph_pt = chain.ph_noniso_pt
		ph_etas2 = chain.ph_noniso_etas2
		ph_phi = chain.ph_noniso_phi
		ph_iso = chain.ph_noniso_iso
		ph_trackiso = chain.ph_noniso_trackiso

		met_et = chain.met_et

		# systematic variation (missing MET variation)
		for i in xrange(ph_n):
			ph_pt[i] *= cfg['alpha'] # proportional to E
			ph_iso[i] *= cfg['alpha'] # topoetcone is proportional to E
			ph_trackiso[i] /= cfg['alpha'] # inverse proportional to E


		pair = get_pair(el_n, el_pt, el_etas2, el_phi, el_ch, ph_n, ph_pt, ph_etas2, ph_phi, ph_iso, ph_trackiso, met_et)


		if pair['TP'][0] is not None:

			eta_bin = get_bin(abs(pair['TP'][4]), cfg['eta_binning'])
			pt_bin = get_bin(pair['TP'][5], cfg['pt_binning'])

			if eta_bin is not None and pt_bin is not None:
				h_mass['TP_'+pair['TP'][1]][eta_bin][pt_bin].Fill(pair['TP'][0])

			if pair['TP'][1] == 'ee':
				n_TP_ee+=1
			elif pair['TP'][1] == 'eg':
				n_TP_eg+=1


		if pair['TT1'][0] is not None:

			reg = get_pair_region(pair['TT1'][2], pair['TT1'][4])

			if pair['TT1'][1] == 'ee':
				n_TT1_ee+=1
				if pair['TT1'][0]>cfg['TT1']['cut_mass'][0] and pair['TT1'][0]<cfg['TT1']['cut_mass'][1]:
					h_grid['ee']['TT1_'+reg].Fill(abs(pair['TT1'][2]), pair['TT1'][3])
					h_grid['ee']['TT1_'+reg].Fill(abs(pair['TT1'][4]), pair['TT1'][5])
			if pair['TT1'][1] == 'eg':
				n_TT1_eg+=1
				if pair['TT1'][0]>cfg['TT1']['cut_mass'][0] and pair['TT1'][0]<cfg['TT1']['cut_mass'][1]:
					h_grid['eg']['TT1_'+reg].Fill(abs(pair['TT1'][4]), pair['TT1'][5])

			h_mass['TT1_'+pair['TT1'][1]][reg].Fill(pair['TT1'][0])


		if pair['TT2'][0] is not None:

			eta_bin = get_bin(abs(pair['TT2'][4]), cfg['eta_binning'])
			pt_bin = get_bin(pair['TT2'][5], cfg['pt_binning'])

			if eta_bin is not None and pt_bin is not None:
				h_mass['TT2_'+pair['TT2'][1]][eta_bin][pt_bin].Fill(pair['TT2'][0])

			if pair['TT2'][1] == 'ee':

				n_TT2_ee+=1

				eta_bin = get_bin(abs(pair['TT2'][2]), cfg['eta_binning'])
				pt_bin = get_bin(pair['TT2'][3], cfg['pt_binning'])

				if eta_bin is not None and pt_bin is not None:
					h_mass['TT2_ee'][eta_bin][pt_bin].Fill(pair['TT2'][0])

			elif pair['TT2'][1] == 'eg':

				n_TT2_eg+=1


	h_mass['TT2_ee_all'].Add(h_mass['TT1_ee']['EE'])
	h_mass['TT2_ee_all'].Add(h_mass['TT1_ee']['BE'])
	h_mass['TT2_ee_all'].Add(h_mass['TT1_ee']['BB'])

	h_mass['TT2_eg_all'].Add(h_mass['TT1_eg']['EE'])
	h_mass['TT2_eg_all'].Add(h_mass['TT1_eg']['BE'])
	h_mass['TT2_eg_all'].Add(h_mass['TT1_eg']['BB'])


	cfg['n_TP_ee'] = n_TP_ee
	cfg['n_TP_eg'] = n_TP_eg
	cfg['n_TT1_ee'] = n_TT1_ee
	cfg['n_TT1_eg'] = n_TT1_eg
	cfg['n_TT2_ee'] = n_TT2_ee
	cfg['n_TT2_eg'] = n_TT2_eg

	output(h_grid, h_mass)
	

def get_pair(el_n, el_pt, el_etas2, el_phi, el_ch, ph_n, ph_pt, ph_etas2, ph_phi, ph_iso, ph_trackiso, met_et):

	pair = {}

	Z_mass = 91.1876 # PDG

	for m in cfg['methods']:

		pair_mass = -999999999.

		pair[m] = []
		pair[m].append(None) 	# mass
		pair[m].append(None)	# pair type
		pair[m].append(None)	# part1 eta
		pair[m].append(None)	# part1 pt
		pair[m].append(None)	# part2 eta
		pair[m].append(None)	# part2 pt

		if met_et < cfg[m]['cut_met_et'][0] or met_et > cfg[m]['cut_met_et'][1]:
			continue

		for i in xrange(el_n):

			if el_pt[i] < cfg[m]['cut_el_pt_1'][0] or el_pt[i] > cfg[m]['cut_el_pt_1'][1]:
				continue

			for j in xrange(el_n):

				if i < j:
					continue

				if el_ch[i] == el_ch[j]:
					continue

				if el_pt[j] < cfg[m]['cut_el_pt_2'][0] or el_pt[j] > cfg[m]['cut_el_pt_2'][1]:
					continue

				tmp_mass = inv_mass(el_pt[i], el_etas2[i], el_phi[i], el_pt[j], el_etas2[j], el_phi[j], 0.000511)

				if abs(tmp_mass - Z_mass) < abs(pair_mass - Z_mass):

					pair_mass = tmp_mass

					pair[m][0] = pair_mass
					pair[m][1] = 'ee'
					pair[m][2] = el_etas2[i]
					pair[m][3] = el_pt[i]
					pair[m][4] = el_etas2[j]
					pair[m][5] = el_pt[j]

			for j in xrange(ph_n):


				if ph_pt[j] < cfg[m]['cut_ph_pt'][0] or ph_pt[j] > cfg[m]['cut_ph_pt'][1]:
					continue

				if ph_iso[j] > cfg[m]['cut_ph_iso'] or ph_trackiso[j] > cfg[m]['cut_ph_trackiso']:
					continue

				tmp_mass = inv_mass(el_pt[i], el_etas2[i], el_phi[i], ph_pt[j], ph_etas2[j], ph_phi[j], 0.)

				if abs(tmp_mass - Z_mass) < abs(pair_mass - Z_mass):

					pair_mass = tmp_mass

					pair[m][0] = pair_mass
					pair[m][1] = 'eg'
					pair[m][2] = el_etas2[i]
					pair[m][3] = el_pt[i]
					pair[m][4] = ph_etas2[j]
					pair[m][5] = ph_pt[j]

	return pair

def get_pair_region(eta1, eta2):

	if eta1 < 1.475 and eta2 < 1.475:
		return 'BB'
	elif eta1 > 1.475 and eta2 > 1.475:
		return 'EE'
	else:
		return 'BE'


def get_bin(val, values):

	for i,v in enumerate(values[:-1]):

		if val >= v and val < values[i+1]:

			return i

	return None



		
def inv_mass(pt1, eta1, phi1, pt2, eta2, phi2, mass2):

	v1 = ROOT.TLorentzVector() 
	v2 = ROOT.TLorentzVector()

	v1.SetPtEtaPhiM( pt1, eta1, phi1, 0.000511)
	v2.SetPtEtaPhiM( pt2, eta2, phi2, mass2)

	return (v1 + v2).M()



def output(h_grid, h_mass):

	filename = 'output/%s/output_loop.root' % (cfg['tag'])
	if cfg['syst_energy']:
		filename = 'output/%s/output_loop_syst_energy_%s.root' % (cfg['tag'], str(cfg['alpha']))
	if cfg['syst_masswin']:
		filename = 'output/%s/output_loop_syst_masswin.root' % cfg['tag']
	
	op_file = TFile(filename, 'RECREATE')

	gStyle.SetOptStat(111111)


	for t in ['ee', 'eg']:

		h_grid[t]['TP'].Write('h_N_TP_'+t)
		h_grid[t]['TT1_EE'].Write('h_N_TT1_EE_'+t)
		h_grid[t]['TT1_BE'].Write('h_N_TT1_BE_'+t)
		h_grid[t]['TT1_BB'].Write('h_N_TT1_BB_'+t)
		h_grid[t]['TT2'].Write('h_N_TT2_'+t)

		for r in ['EE', 'BE', 'BB']:
			h_mass['TT1_'+t][r].Write('h_m_TT1_%s_%s'%(t, r))

		for eta in xrange(len(cfg['eta_binning'])-1):
			for pt in xrange(len(cfg['pt_binning'])-1):

				h_mass['TP_'+t][eta][pt].Write('h_m_TP_'+t+'_%i_%i'%(eta, pt))
				h_mass['TT2_'+t][eta][pt].Write('h_m_TT2_'+t+'_%i_%i'%(eta, pt))

	h_mass['TT2_ee_all'].Write('h_m_TT2_ee_all')
	h_mass['TT2_eg_all'].Write('h_m_TT2_eg_all')


	op_file.Close()


	print ('%s was created' % filename)

	if not cfg['syst_energy']:

		with open('output/' + cfg['tag'] +'/used_config.yaml', 'w+') as f:
			data = yaml.dump(cfg, f)
		print ('output/' + cfg['tag'] +'/used_config.yaml was created')










def print_event(chain, entry):

	chain.GetEntry(entry)

	jet_n = chain.jet_n

	el_n = chain.el_medium_n
	ph_n = chain.ph_noniso_n
	el_pt = chain.el_medium_pt
	el_etas2 = chain.el_medium_etas2
	el_phi = chain.el_medium_phi
	el_ch = chain.el_medium_ch

	ph_pt = chain.ph_noniso_pt
	ph_etas2 = chain.ph_noniso_etas2
	ph_phi = chain.ph_noniso_phi
	ph_iso = chain.ph_noniso_iso
	ph_trackiso = chain.ph_noniso_trackiso

	met_et = chain.met_et
	ht = chain.ht

	print 'Entry: %i' % entry
	print 'el_n: %i' % el_n
	for i in xrange(el_n):
		print 'el %i:' % i,
		print '\t el_pt: %f' % el_pt[i]
		print '\t el_etas2: %f' % el_etas2[i]
		print '\t el_phi: %f' % el_phi[i]
		print '\t el_ch: %f' % el_ch[i]
	print 'ph_n: %i' % ph_n
	for i in xrange(ph_n):
		print 'el %i:' % i,
		print '\t ph_pt: %f' % ph_pt[i]
		print '\t ph_etas2: %f' % ph_etas2[i]
		print '\t ph_phi: %f' % ph_phi[i]
		print '\t ph_iso: %f' % ph_iso[i]
		print '\t ph_trackiso: %f' % ph_trackiso[i]
	print 'jet_n: %i' % jet_n
	print 'met_et: %f' % met_et
	print 'ht: %f' % ht
	print ''

