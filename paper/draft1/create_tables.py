import pandas as pd
import numpy as np
import os

if __name__ == 'main':

	file = 'data_files/mann_results.csv'
	cont = pd.read_csv(file)

	des = cont['ID']
	ra  = cont['RA']
	dec = cont['DEC']
	trn_teff = cont['TRAIN_TEFF']
	tst_teff = cont['TEST_TEFF']
	crv_teff = cont['CROSS_TEFF']
	trn_feh = cont['TRAIN_FEH']
	tst_feh = cont['TEST_FEH']
	crv_feh = cont['CROSS_FEH']
		
	with open('tables/mann_results.txt', 'w') as wf:
		for i in range(len(des)):
			line = "{} & {} & {} & {} & {} & {} & {} & {} & {} \\".format(des[i], ra[i], dec[i], \
				trn_teff[i], tst_teff[i], crv_teff[i], \
				trn_feh[i], tst_feh[i], crv_feh[i])

			print(line)
			wf.write(line)
