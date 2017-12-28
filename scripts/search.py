import os
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits, ascii
from astroquery.simbad import Simbad

db_path = os.environ['APOGEE_DATA'] 

def searchDatabase(**kwargs):

    """
    Searching DR13 Database for stars
    Keywords are explained in the Datamodel file
    https://goo.gl/5hcygL
    Author: Christian Aganze
    """

    database = db_path + '/allStar-l30e.2.fits'

    search_par = kwargs.get('par', 'TEFF')
    select     = kwargs.get('select', [3000, 4500])
    save       = kwargs.get('save', True)
    output     = kwargs.get('out', str(select) + '.csv')

    data = Table(fits.open(database)[1].data)
    p = data[search_par]

    if isinstance(select[0], type('k')):
        # matching by strings
        print("Matching by ", search_par)
        condition = [i for i in range(0, len(data[search_par]))  if data[search_par][i] in select]

    else:
        print(" ALL stars with ", search_par, "between", select)
        condition = np.where((p>=select[0] ) & (p<=select[1])) [0]

    select_data = data[condition]

    data_dict  = {'2MASS_ID':select_data['APOGEE_ID'], 'TEFF':select_data['TEFF'], 'LOGG':select_data['LOGG'], 'M_H':select_data['M_H']}
    data_table = pd.DataFrame(data=data_dict)

    if save == True:
        data_table.to_csv('tables/'+output)
   
    return data[condition]


def multiParamSearch(**kwargs):

    """
    Search APOGEE database.
    Return table with source names and aspcap fit parameters
    """

    search_par = kwargs.get('par', ['TEFF'])
    select     = kwargs.get('select', [[3000, 4500]])
    save       = kwargs.get('save', True)
    output     = kwargs.get('out', str(select) + '.csv')
    release    = kwargs.get('rel', ['dr13'])

    data_releases = {'dr10':'allStar-v304.fits', 'dr11':'allStar-v402.fits', \
    'dr12':'allStar-v603.fits', 'dr13':'allStar-l30e.2.fits', 'dr14':'allStar-l31c.2.fits'}

    database = db_path + '/allStar-l30e.2.fits'
    data = Table(fits.open(database)[1].data)

    for i in range(len(search_par)):
        p = data[search_par[i]]
        if isinstance(select[i][0], type('k')):
            # matching by strings
            condition = [j for j in range(0, len(data[search_par[i]]))  if data[search_par[i]][j] in select[i]]

        else:
            print(" ALL stars with ", search_par[i], "between", select[i])
            condition = np.where((p>=select[i][0] ) & (p<=select[i][1]))[0]
        data = data[condition]

    # Put only interesting values into the data table: 2MASS name and aspcap parameters
    data_dict  = {'2MASS_ID':data['APOGEE_ID'], 'TEFF':data['TEFF'], 'LOGG':data['LOGG'], 'M_H':data['M_H']}
    # data_dict = {'2MASS_ID':data['APOGEE_ID'], 'RA':data['RA'], 'DEC':data['DEC']}
    data_table = pd.DataFrame(data=data_dict)

    # Concatenate frames of all data release searches and save
    if save == True:
        data_table.to_csv('tables/'+output)

    return data


def returnAspcapTable(tm_ids, **kwargs):

    """
    Return a table of parameters (in csv format) for a list of APOGEE spectra by 2MASS name
    """

    default_params = [  'APOGEE_ID', 'TEFF', 'LOGG', 'M_H', \
                        'J', 'J_ERR', 'H', 'H_ERR', 'K', 'K_ERR', \
                        'WASH_M', 'WASH_M_ERR', 'WASH_T2', 'WASH_T2_ERR', \
                        'DDO51', 'DDO51_ERR', \
                        'WASH_DDO51_GIANT_FLAG', 'WASH_DDO51_STAR_FLAG', \
                        'RA', 'DEC', 'SNR' ]

    # optional
    params = kwargs.get('par', default_params)
    save   = kwargs.get('save', True)
    output = kwargs.get('out', 'aspcap_table.csv')

    database = db_path + '/allStar-l30e.2.fits'
    data = Table(fits.open(database)[1].data)

    index = [np.where(data['APOGEE_ID'] == TMID)[0][0] for TMID in tm_ids]

    data_dict  = {params[i]: data[params[i]][index] for i in range(len(params))}
    data_table = pd.DataFrame(data=data_dict)

    if save == True:
        data_table.to_csv('tables/'+output)

    return data_table


def readSavedChiValues(**kwargs):

    """
    Read saved list of chi-squared comparison values created from sortMStars() in spec.py
    """

    # optional
    path = kwargs.get('file', 'sort_output/unkown_star_chivals_1000.csv')

    values = pd.read_csv(path)



def returnSimbadParams(id_list):

    customSimbad = Simbad()
    customSimbad.add_votable_fields('fe_h', 'rv_value', 'rvz_error', 'sptype', 'rot')
    customSimbad.get_votable_fields()

    sample = pd.read_csv(id_list)
    objects = []

    no_entry = 0
    for s in sample['ID']:
        ID = str(s)
        obj = customSimbad.query_object(ID)

        try:
            spt = str(obj['SP_TYPE']).split('-\n')[1]
            d = {'ID':str(s), 'SPT':spt, 'TEFF':obj['Fe_H_Teff'], 'LOGG':obj['Fe_H_log_g'], \
                 'FE_H':obj['Fe_H_Fe_H']} 
            df = pd.DataFrame(data=d)
            objects.append(df)
        except:
            no_entry += 1

    print('No simbad values for ' + no_entry + ' objects.')

    result = pd.concat(objects)

    # Query:
    # sptypes >= M0 & sptypes <= M9 & fe_h.teff >=2000 & fe_h.logg >= 3 & fe_h.metidx <= 1

    return result


def returnTMID(id_list):

    """
    Input list of sources with miscellaneous identifiers; 
    output those with 2MASS IDs with their spectral types
    """

    TM_ID = []

    for identifier in id_list:

        all_names = Simbad.query_objectids(identifier)
        
        for name in all_names['ID']:

            if name[0:7] == b'2MASS J':
                TM_ID.append(name)

    d = {'ID': TM_ID}
    df = pd.DataFrame(data=d)

    df.to_csv('tables/2MASS_IDs.csv')

    return TM_ID


def returnTMIDandSPT(source_list):

    """
    Input file of list of sources with miscellaneous identifiers; 
    output those with 2MASS IDs with their spectral types
    """

    simbad_sample = pd.read_csv(source_list)

    id_list  = simbad_sample['ID']
    spt_list = simbad_sample['SPT']

    TM_ID = []
    SPTS  = []

    for i in range(len(id_list)):

        all_names = Simbad.query_objectids(id_list[i])
        spt = spt_list[i]
        
        for j in range(len(all_names['ID'])):

            name = all_names['ID'][j]

            if name[0:7] == b'2MASS J':
                TM_ID.append(name)
                SPTS.append(spt)

                print(name, spt)

    d = {'ID': TM_ID, 'SPT': SPTS}
    df = pd.DataFrame(data=d)

    df.to_csv('tables/2MASS_ID_and_SPT.csv')

    return TM_ID, SPTS


def inAPOGEE(source_list):

    """
    Intake list of 2MASS names. Return csv file of those with spectra in APOGEE
    """

    source_info = pd.read_csv(source_list)

    id_list = source_info['ID']

    ap_names, ap_spts = [], []
    for i in range(len(id_list)):

        name = source_info['ID'][i]
        spt = source_info['SPT'][i]

        try:
            ap.Spectrum(id=name, type='aspcap')
            ap_names.append(name)
            ap_spts.append(spt)

        except:
            print(name, 'not in APOGEE.')

    d = {'ID': ap_names, 'SPT': ap_spts}
    df = pd.DataFrame(data=d)

    df.to_csv('tables/APOGEE_specs_with_SPT.csv')

    return ap_names, ap_spts

