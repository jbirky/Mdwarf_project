import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from operator import itemgetter
import os
import math

from astropy.table import Table
from astropy.io import fits, ascii
from astropy import units as u

AP_PATH = os.environ['APOGEE_DATA']

import apogee_tools as ap

from plot import *


def compareToModels(**kwargs):

	"""
	This function compares a spectrum to two models and outputs chi-squared values
	The default M dwarf model is 2M19213157+4317347 and 
	the default M giant model is 2M00001653+5540107
	"""

	default_spec1 = ap.Spectrum(id='2M19213157+4317347', type='aspcap')
	default_spec2 = ap.Spectrum(id='2M00001653+5540107', type='aspcap')

	default_params = [['4079','4.75','-0.20'],['4080','1.57','-0.15']]
	default_uncert = [['238', '0.15', '0.10 *'],['70', '0.08', '0.02']]

	spec   = kwargs.get('spec')
	slbl   = kwargs.get('slbl')
	models = kwargs.get('models', [default_spec1, default_spec2])
	params = kwargs.get('params', default_params)
	uncert = kwargs.get('unc', default_uncert)
	regs   = kwargs.get('regs', [[15650,15780], [16150,16280]]) #for plotting
	plot   = kwargs.get('plot', False)

	mdl1, mdl2 = models[0], models[1]
	par1, par2 = params[0], params[1]
	unc1, unc2 = uncert[0], uncert[1]

	chi1, spec, mdl1 = ap.compareSpectra(spec, mdl1)
	chi2, spec, mdl2 = ap.compareSpectra(spec, mdl2)

	if plot == True:
		plotChiComparison(spec=spec, models=models, params=params, unc=uncert, \
			chis=[chi1, chi2], regs=regs, save=True, slbl=slbl)

	print(spec.name, chi1, chi2)

	return chi1, chi2


def sortMStars(**kwargs):

	"""
	Enter up to 3 lists of spectrum apogee IDs, return list of chi-squared values 
	for two models and plot (optional).

	Input:  Enter 3 separate lists to plot in different colors:
			1) dwarf list = blue, 2) giant list = green, 3) unknown list = black

	Output: 3 lists of chi values : d_chivals, g_chivals, u_chivals
			where each element of these lists is [chi_dwarf_template, chi_giant_template]
	"""

	dnames = kwargs.get('dnames')
	gnames = kwargs.get('gnames')
	unames = kwargs.get('unames')

	if not os.path.exists('sort_output/'):
		os.makedirs('sort_output/')

	# Compute chi values of list of dwarfs
	if 'dnames' in kwargs:

		d_chivals = []
		for name in dnames:
			try:
				sp = ap.Spectrum(id=name, type='aspcap')
				chis = compareToModels(spec=sp)
			except:
				print('Error:' + name + ' not in downloaded files.')

			d_chivals.append(chis)

		# Separate chi value pairs
		dd_chis = list(map(itemgetter(0), d_chivals))
		dg_chis = list(map(itemgetter(1), d_chivals))

		# Save chi values to 3 separate csv file
		d_dict = {'ID': dnames, 'dchi': dd_chis, 'gchi': dg_chis}
		d_dataframe = pd.DataFrame(d_dict)
		d_dataframe.to_csv('sort_output/known_dwarf_chivals.csv')

	else:
		d_chivals = None

	# Compute chi values of list of giants
	if 'gnames' in kwargs:

		g_chivals = []		
		for name in gnames:
			try:
				sp = ap.Spectrum(id=name, type='aspcap')
				chis = compareToModels(spec=sp)
			except:
				print('Error:' + name + ' not in downloaded files.')

			g_chivals.append(chis)

		gd_chis = list(map(itemgetter(0), g_chivals))
		gg_chis = list(map(itemgetter(1), g_chivals))

		g_dict = {'ID': gnames, 'dchi': gd_chis, 'gchi': gg_chis}
		g_dataframe = pd.DataFrame(g_dict)
		g_dataframe.to_csv('sort_output/known_giant_chivals.csv')

	else:
		g_chivals = None

	# Compute chi values of list of unknown M stars

	if 'unames' in kwargs:

		u_chivals = []		
		for name in unames:
			try:
				sp = ap.Spectrum(id=name, type='aspcap')
				chis = compareToModels(spec=sp)
			except:
				print('Error:' + name + ' not in downloaded files.')

			u_chivals.append(chis)

		ud_chis = list(map(itemgetter(0), u_chivals))
		ug_chis = list(map(itemgetter(1), u_chivals))

		u_dict = {'ID': unames, 'dchi': ud_chis, 'gchi': ug_chis}
		u_dataframe = pd.DataFrame(u_dict)
		u_dataframe.to_csv('sort_output/unkown_star_chivals.csv')

	else:
		u_chivals = None

	# Now plot the scatter; plot automatically saved to plots directory as 'Chi_Squared_Scatter.pdf'
	plotChiComparisonScatter(dlist=d_chivals, glist=g_chivals, ulist=u_chivals)

	return d_chivals, g_chivals, u_chivals


def readSortedMStars(**kwargs):

	d_chivals, g_chivals, u_chivals = [], [], []

	plotChiComparisonScatter(dlist=d_chivals, glist=g_chivals, ulist=u_chivals)

	return d_chivals, g_chivals, u_chivals


def spectralSequence(md, **kwargs):

	"""
	Format spectra to plot spectral sequence.

	Input:  'md' : Cannon model
	"""



