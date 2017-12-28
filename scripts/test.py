import apogee_tools as ap
from spec import *
from plot import *
from search import *
from run_cannon import *
import os

os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2016/bin/x86_64-darwin'

from datetime import datetime
startTime = datetime.now()


if __name__ == '__main__':

# Search for list of IDs to see if they are in APOGEE

	# simbad_specs = pd.read_csv('tables/teff_simbad_ids.csv')
	# in_apogee = []
	# for spec in simbad_specs['ID']:
	# 	if len(spec) != 0:
	# 		in_apogee.append(ap.searchStars(id_name=spec)[0])
	# print(in_apogee)

#--------------------------------------------------------------------
# Download sources

	# ap.download('2M21354275-0113120', type='aspcap')

	# data = ap.Spectrum(id='2M21354275-0113120', type='aspcap')
	# print(data.params)
	# data.plot(save=True)

#--------------------------------------------------------------------
# Mass download

	# source_list = pd.read_csv('/Users/admin/Desktop/Research/Mdwarf_project/tables/2MASS_ID_and_SPT.csv')['ID']

	# for source in source_list:
	# 	ap.download(source, type='aspcap', dir='/Users/admin/Desktop/Research/cannon/3_west_sample/spec_sequence')

	# ap_names, ap_spts = inAPOGEE('tables/2MASS_ID_and_SPT.csv')

#--------------------------------------------------------------------
# Plot spectra in regions

	# spc1 = ap.Spectrum(id='2M19213157+4317347', type='aspcap')
	# spc2 = ap.Spectrum(id='2M00001653+5540107', type='aspcap')

	par1 = ['4079','4.75','-0.20']
	par2 = ['4080','1.57','-0.15']

	unc1 = ['238', '0.15', '0.10 *']
	unc2 = ['70', '0.08', '0.02']

	xlim = [15190, 15200]
	# out  = 'FeI_comparison.pdf'

	# plotRegion(specs=[spc1,spc2], params=[par1,par2], unc=[unc1,unc2], save=True, xlim=xlim)
	
	# plotMultiRegion(specs=[spc1,spc2], params=[par1,par2], unc=[unc1,unc2], save=True, \
		# dir='plots/', nplots=2, regs=[[15650,15780], [16150,16280]])

	# plotMultiSpec(specs=[[spc1,spc2]], params=[[par1,par2]], unc=[[unc1,unc2]], save=True, \
	# 	dir='plots/', dim=[1,2], regs=[[15650,15780], [16150,16280]])

	# dnames = ['2M11052903+4331357','2M11032023+3558117','2M19543665+4357180','2M17362594+6820220','2M19211069+4533525','2M18543080+4823277','2M20531977+6209156','2M09142485+5241118','2M09142298+5241125','2M19343286+4249298','2M19213157+4317347','2M19320191+4218259','2M18444674+4729496','2M00182549+4401376','2M19300081+4304593','2M19283288+4225459','2M19312949+4103513','2M05312734-0340356','2M13454354+1453317','2M23315244+1956138','2M00182256+4401222','2M02001278+1303112','2M02361535+0652191','2M04311147+5858375','2M05420897+1229252','2M06000351+0242236','2M06011106+5935508','2M06544902+3316058','2M07272450+0513329','2M07444018+0333089','2M08585633+0828259','2M10121768-0344441','2M10193634+1952122','2M10285555+0050275','2M10505201+0648292','2M11000432+2249592','2M11474440+0048164','2M13295979+1022376','2M16252459+5418148','2M17093153+4340531','2M22464980+4420030','2M22563497+1633130','2M23315208+1956142','2M23491255+0224037']
	# gnames = ['2M00000317+5821383', '2M00002227+6223341', '2M00002903+6403157', '2M00004297+5739255', '2M00010832+6307323', '2M00011022+6346190', '2M00012412+6427175', '2M00014487-0001080', '2M00015350+6459174', '2M00021937+1545415', '2M00022256+0113360', '2M00024579+7352118', '2M00025558+7034422', '2M00031070+6221346', '2M00032059+7546359', '2M00032797+0114040', '2M00032999+7451590', '2M00033334+5917055', '2M00033942+7335499', '2M00040831+7440077', '2M00042226+7113282', '2M00042543+7127105', '2M00051500+7131319', '2M00054450+7028061', '2M00055073+7023018', '2M00060972+7516165', '2M00062020-0003149', '2M00065212+5821138', '2M00065881+7453184', '2M00070426+5838505', '2M00071394+7533482', '2M00072200+7400230', '2M00072283+5904090', '2M00074235+7536162', '2M00074572+5750059', '2M00075614+5851490', '2M00080066+5825585', '2M00081288+5825054', '2M00081458+5826219', '2M00082377+5827043', '2M00090845+6927361', '2M00092683+7323498', '2M00093033+5832102', '2M00093509+7440151', '2M00094558+7559454']
	# unames = pd.read_csv('tables/mstars.csv')['ID'][0:1000]

	# sortMStars(dnames=dnames, gnames=gnames, unames=unames)
	# print(datetime.now() - startTime)

	# d_slbl = [['4411', '4.56', '-0.19'], ['..', '..', '.. *']]
	# d_spec = ap.Spectrum(id='2M18444674+4729496', type='aspcap')
	# g_slbl = [['3817', '0.97', '-0.27'], ['..', '..', '..']]
	# g_spec = ap.Spectrum(id='2M00004297+5739255', type='aspcap')
	# compareToModels(spec=d_spec, plot=True, slbl=d_slbl)

#--------------------------------------------------------------------
# Search for sources

	params = ['TEFF']
	ranges = [[-10000,4500]]
	rel = ['dr13']

	# source_table = multiParamSearch(par=params, select=ranges, rel=rel)

	# data = pd.read_csv('tables/mgiants_apogee_search.csv')
	# print(data['2MASS_ID'][0:45])

	# source_table = searchDatabase(par='TEFF', select=[-10000,4000])

#--------------------------------------------------------------------
# Create table of aspcap data table info

	# test_sample = pd.read_csv('/Users/admin/Desktop/Research/Mdwarf_project/tables/mstars.csv')['ID'][0:1000]
	# returnAspcapTable(source_list)

#--------------------------------------------------------------------
# Simbad search

	# from astroquery.simbad import Simbad
	# result = Simbad.query_criteria('sptypes >= M0 & sptypes <= M9 & fe_h.teff <= 4000 & fe_h.metidx <= 5')
	# names = list(result['MAIN_ID'])

	# TM_ID, SPTS = returnTMIDandSPT('/Users/admin/Desktop/Research/Mdwarf_project/mstar_searches/msources_with_feh.csv')
	
	# TM_ID = returnTMID(names)

	# dtable = returnSimbadParams('/Users/admin/Desktop/Research/MPIA/simbad_search/muirhead2014.csv')

#--------------------------------------------------------------------
# Run the Cannon

	data_path = '/Users/admin/Desktop/Research/cannon/3_west_sample/Data/'
	ref_path  = '/Users/admin/Desktop/Research/cannon/3_west_sample/west_ref_labels.csv'

	# ds = initializeTrainingSet(data=data_path, ref=ref_path, label_names=['SPT'])

	# trn_labels, crv_labels = crossValidate(ds, label_names=['SPT'])

	# print(trn_labels)
	# print(crv_labels)

	# trn_labels = [[  3.55100000e+03,  -2.80000000e-01],[  3.63000000e+03,  -8.80000000e-01],[  3.27600000e+03,  -1.50000000e-01],[  3.24700000e+03,  -2.40000000e-01],[  3.27200000e+03,   3.60000000e-01],[  3.62600000e+03,   6.00000000e-01],[  3.16700000e+03,  -2.50000000e-01],[  3.17800000e+03,   8.00000000e-02],[  3.34000000e+03,  -9.00000000e-02],[  3.44800000e+03,  -2.00000000e-02],[  3.31700000e+03,  -1.10000000e-01],[  3.04500000e+03,   4.00000000e-01],[  3.26300000e+03,  -1.10000000e-01],[  3.86800000e+03,  -3.80000000e-01],[  3.76900000e+03,  -4.00000000e-01],[  3.61600000e+03,   1.70000000e-01],[  3.37000000e+03,   1.50000000e-01],[  3.52900000e+03,  -2.80000000e-01],[  3.31400000e+03,   2.90000000e-01],[  3.52600000e+03,  -9.00000000e-02],[  3.55000000e+03,  -4.00000000e-01],[  3.66400000e+03,  -8.30000000e-01],[  3.28800000e+03,  -8.00000000e-02],[  3.89600000e+03,  -1.20000000e-01],[  3.62300000e+03,  -1.00000000e-01],[  3.47500000e+03,  -3.50000000e-01],[  3.21800000e+03,   3.30000000e-01],[  3.58300000e+03,   0.00000000e+00],[  4.41100000e+03,  -1.90000000e-01],[  4.23700000e+03,  -2.30000000e-01],[  4.26100000e+03,  -1.80000000e-01],[  4.07900000e+03,  -2.00000000e-01],[  4.59000000e+03,   3.40000000e-01],[  4.55900000e+03,   7.00000000e-02],[  3.82000000e+03,   2.00000000e-01],[  4.11800000e+03,  -7.10000000e-01],[  4.50000000e+03,  -1.80000000e-01],[  4.04800000e+03,  -2.10000000e-01],[  3.66000000e+03,  -1.30000000e-01],[  3.10600000e+03,  -7.00000000e-02],[  3.76400000e+03,   3.00000000e-01],[  3.35300000e+03,   3.00000000e-02],[  3.20000000e+03,   9.00000000e-02],[  3.64600000e+03,  -4.50000000e-01]]
	# crv_labels = [[  3.43104991e+03,   4.83612287e-01], [  3.63922866e+03,  -4.99350443e-01], [  3.04815024e+03,   6.08638279e-02], [  3.07772087e+03,   2.66969666e-02], [  3.68154067e+03,  -7.44568934e-01], [  3.28234139e+03,  -1.14931034e-01], [  4.05446683e+03,   3.30979194e-02], [  3.12156122e+03,  -8.78026876e-02], [  3.14203861e+03,   4.16061299e-02], [  3.22496554e+03,   9.16231171e-02], [  3.40024025e+03,  -7.67238720e-03], [  3.22582072e+03,  -1.39836620e-01], [  3.13354355e+03,   6.86924068e-03], [  3.11375278e+03,   2.76132945e-01], [  4.10167459e+03,  -2.72643886e-02], [  4.30947078e+03,  -4.65434775e-02], [  3.81491826e+03,  -1.99211257e-03], [  3.34794188e+03,  -4.20072612e-02], [  3.59850135e+03,  -4.68917409e-01], [  3.27027780e+03,   9.02255452e-02], [  3.42336764e+03,   1.08702842e-01], [  3.51385810e+03,  -2.89298123e-01], [  3.64757360e+03,  -7.84459625e-01], [  3.09375430e+03,   5.93125364e-02], [  3.74671886e+03,  -5.73158000e-01], [  3.66082284e+03,  -8.14267249e-01], [  3.57648896e+03,  -5.17530035e-01], [  3.25923008e+03,   6.77584766e-02], [  4.68462772e+03,   5.92506977e-04], [  3.76331520e+03,  -2.67081833e-01], [  4.38619839e+03,  -2.81642547e-02], [  3.91635482e+03,  -4.15530185e-01], [  4.43742402e+03,  -1.22815459e-01], [  4.45499116e+03,  -1.73711238e-01], [  4.35534103e+03,   2.86814598e-02], [  4.54148322e+03,   4.34857946e-02], [  4.29803256e+03,  -3.10781253e-02], [  3.84618960e+03,  -6.86506940e-01], [  3.90240498e+03,  -3.72110119e-03], [  3.16627892e+03,  -4.25012530e-02], [  3.93806336e+03,  -8.01379106e-02], [  3.21567960e+03,   7.36737530e-02], [  3.30043750e+03,   2.47389253e-02], [  3.60238114e+03,  -8.46320992e-01]]

	# plotCrossValidation(trn_labels, crv_labels, label_names=['SPT'])

#--------------------------------------------------------------------
# Color cross validation by parameter value

	default_params = [	'APOGEE_ID', 'TEFF', 'LOGG', 'M_H', \
						'J', 'J_ERR', 'H', 'H_ERR', 'K', 'K_ERR', \
						'WASH_M', 'WASH_M_ERR', 'WASH_T2', 'WASH_T2_ERR', \
						'DDO51', 'DDO51_ERR', \
						'WASH_DDO51_GIANT_FLAG', 'WASH_DDO51_STAR_FLAG', \
						'RA', 'DEC', 'SNR' ]	

	# chival_table = pd.read_csv('sort_output/unknown_star_chivals_1000.csv')
	# tm_ids = chival_table['ID']

	# aspcap_table = returnAspcapTable(tm_ids, save=True, out='west_sample_info.csv', params=default_params)

	# plot_param = aspcap_table['SNR']

	# plotChiComparisonScatterByParam(dchi=chival_table['dchi'], gchi=chival_table['gchi'], \
	# 		cparam=plot_param)

#--------------------------------------------------------------------
# 4 Spectra Plots

	ids = np.array(['2M19543665+4357180','2M02361535+0652191','2M05312734-0340356','2M23491255+0224037'])

	sample1 = ap.Spectrum(id=ids[0], type='aspcap')
	sample2 = ap.Spectrum(id=ids[1], type='aspcap')
	sample3 = ap.Spectrum(id=ids[2], type='aspcap')
	sample4 = ap.Spectrum(id=ids[3], type='aspcap')	

	par1 = ['4048', '4.75', '-0.21']
	par2 = ['3247', '4.5', '-0.24']
	par3 = ['3626', '4.8', '0.6']
	par4 = ['3646', '4.5', '-0.45']

	# plotFourSpecComparison(specs=[sample1, sample2, sample3, sample4], \
	# 	params=[par1, par2, par3, par4])

	# ds = initializeTrainingSet(data=data_path, ref=ref_path)
	# md, synth_fluxes, test_labels = runCannon(ds)

	# index = [i for i, item in enumerate(ids) if item in set(ds.tr_ID)]

	# model1 = ap.Spectrum(wave=ds.wl, flux=synth_fluxes[index[0]])
	# model2 = ap.Spectrum(wave=ds.wl, flux=synth_fluxes[index[1]])
	# model3 = ap.Spectrum(wave=ds.wl, flux=synth_fluxes[index[2]])
	# model4 = ap.Spectrum(wave=ds.wl, flux=synth_fluxes[index[3]])

	# lbl1 = test_labels[index[0]]
	# lbl2 = test_labels[index[1]]
	# lbl3 = test_labels[index[2]]
	# lbl4 = test_labels[index[3]]

	# plotFourSpecComparisonWithModels(\
	# 	specs =[sample1, sample2, sample3, sample4], \
	# 	params=[par1, par2, par3, par4], \
	# 	models=[model1, model2, model3, model4], \
	# 	labels=[lbl1, lbl2, lbl3, lbl4])

#--------------------------------------------------------------------
# Plot Parameter sensitive regions with BTSettl/PHOENIX models

	# source_index = 2
	# data_path = '/Users/admin/Desktop/Research/cannon/rajpurohit_sample/Data/'
	# ref_path  = '/Users/admin/Desktop/Research/cannon/rajpurohit_sample/rajpurohit_ref_labels.csv'
	# source_info = pd.read_csv(ref_path)

	# data = ap.Spectrum(id=source_info['ID'][source_index], type='aspcap')

	# ds = initializeTrainingSet(data=data_path, ref=ref_path)
	# md, synth_fluxes, test_labels = runCannon(ds)
	# cannon_mdl = ap.Spectrum(wave=ds.wl, flux=synth_fluxes[source_index], params=test_labels[source_index])

	# bts_mdl = ap.getModel(params=[3100, 5.5, 0.0], grid='BTSETTLb', xrange=[15200,16940])

	# plotOHBands(specs=[data, cannon_mdl, bts_mdl])

#--------------------------------------------------------------------
# Fit Cannon models to test set; Plot spectral sequence

	# west_data_path = '/Users/admin/Desktop/Research/cannon/3_west_sample/Data/'
	# west_ref_path  = '/Users/admin/Desktop/Research/cannon/3_west_sample/west_ref_labels.csv'

	# seq_data_path = '/Users/admin/Desktop/Research/cannon/3_west_sample/spec_sequence/Data/'
	# seq_ref_path  = None

	# tr_ds = initializeTrainingSet(data=west_data_path, ref=west_ref_path, label_names=['SPT'])
	# te_ds = initializeTrainingSet(data=seq_data_path, ref=seq_ref_path, label_names=['SPT'])

	# sptypes = [0,1,2,3,4,5,6,7,8,9]
	# te_ds.tr_label = sptypes

	# test_specs, test_models = fitCannonModels(tr_ds, te_ds)
	
	# plotSpectralSequence(specs=test_specs, models=test_models, labels=sptypes)

#--------------------------------------------------------------------
# Plot figure from Rajpurohit paper with model(s)

	spec_id = '2M18451027+0620158'
	raj_data_path = '/Users/admin/Desktop/Research/cannon/2_rajpurohit_sample/Data/'
	raj_ref_path  = '/Users/admin/Desktop/Research/cannon/2_rajpurohit_sample/rajpurohit_ref_labels.csv'

	source_info  = pd.read_csv(raj_ref_path)
	source_index = np.where(source_info['ID'] == spec_id)[0][0]
	data = ap.Spectrum(id=source_info['ID'][source_index], type='aspcap')

	# my_model = ap.getModel(params=[3100,5.5,0.0])
	# plotRajBands(specs=[data, my_model], title='Data vs PHOENIX')

	# ds = initializeTrainingSet(data=raj_data_path, ref=raj_ref_path, lbl_names=['TEFF', 'M_H'], save_base='raj')
	# md, synth_fluxes, test_labels = runCannon(ds)
	# cannon_mdl = ap.Spectrum(wave=ds.wl, flux=synth_fluxes[source_index], params=test_labels[source_index])
	# plotRajBands(specs=[data, cannon_mdl], title='Rajpurohit-trained Cannon Model')

	# phn_mdl = interpolateGrids(irange=[[3000,4000],[5.5,5.5],[0.0,0.5]], grid='PHOENIX', step=[100, .1, .1])[0]
	# bts_mdl = interpolateGrids(irange=[[3100,3100],[5.5,5.5],[0.0,0.0]], grid='BTSETTLb', step=[100, .1, .1])[0]

	wave1 = np.load('normalized_specs/PHOENIX_wl.npy')
	flux1 = np.load('normalized_specs/PHOENIX_t3000_4000_l4.0_5.5_z-0.5_0.5_flux.npy')
	model_norm = ap.Spectrum(wave=wave1, flux=flux1[0], name=[3000,5.5,0.0])

	wave2 = np.load('normalized_specs/PHOENIX_synth_wave.npy')
	flux2 = np.load('normalized_specs/PHOENIX_synth_fluxes.npy')
	phn_mdl = ap.Spectrum(wave=wave2, flux=flux2[0], name=[3000,5.5,0.0])

	# plotRajBands(specs=[model_norm, phn_mdl], title='PHOENIX-trained Cannon Model')

	# plotRajBands(specs=[data, phn_mdl], title='BTSETTL-trained Cannon Model')


	


