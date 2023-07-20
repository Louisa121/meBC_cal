import os
import pickle
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pltl
import seaborn as sns
import matplotlib.dates as mdate
from scipy import stats
import math
import PyMieScatt as ps
from matplotlib import ticker, cm
from scipy.optimize import curve_fit
from scipy.integrate import trapz
from sympy.solvers import solve
from sympy import Symbol
from matplotlib.ticker import LogFormatter 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,NullFormatter,ScalarFormatter,
                               AutoMinorLocator)
import cal_models.Models as models
# from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange 


fileid = 4
fBrC = 0.999


UHSAS_filename = ['RF02_1','RF02_2','RF03','RF04','RF05_1','RF05_2','RF05_3','RF06_1','RF06_2','RF07_1','RF07_2','RF08_1','RF08_2','RF09','RF10','RF11','RF12','RF13']
BC_filename = [
'ORACLES_20180930.txt',
'ORACLES_20181002.txt',
'ORACLES_20181003.txt',
'ORACLES_20181005.txt',
'ORACLES_20181007.txt',
'ORACLES_20181010.txt',
'ORACLES_20181012.txt',
'ORACLES_20181015.txt',
'ORACLES_20181017.txt',
'ORACLES_20181019.txt',
'ORACLES_20181021.txt',
'ORACLES_20181023.txt',
]
data_dir = r'/msegalrolab/luzhang/ORACLES/2018/'
out_dir = r'//msegalrolab/luzhang/Mie_simulations/Figures/ORACLES/2018/'
dataout_dir = '/msegalrolab/luzhang/Mie_simulations/Data/ORACLES/2018/'
BC_dir = r'/msegalrolab/luzhang/ORACLES/SP2_data_for_filters'
Shell_info_dir = r'/msegalrolab/luzhang/ORACLES/SP2_data_for_filters/Coatingthicknesses_sizemodes_consistent.xlsx'
# data_dir = r'/Users/zhanglu/Postdoc/TAU/Research/ORACLES/2018/'
# out_dir = r'/Users/zhanglu/Postdoc/TAU/Research/Mie_simulations/Figures/ORACLE/2018/'
# dataout_dir = '/Users/zhanglu/Postdoc/TAU/Research/Mie_simulations/Data/ORACLE/2018/'
# BC_dir = r'/Users/zhanglu/Postdoc/TAU/Research/ORACLES/SP2_data_for_filters'
# Shell_info_dir = r'/Users/zhanglu/Dropbox/ORACLES/SP2_data_for_filters/Coatingthicknesses_sizemodes_consistent.xlsx'
# Shell_info = pd.read_excel(Shell_info_dir)
full_path_uhsasbin = os.path.join(data_dir,'GITUHSAS_bin.csv')
uhsasbin = pd.read_csv(full_path_uhsasbin)
filter_path = os.path.join(data_dir,'part_data_all_shortEDX_minus_RF2_2_w_new_traj_type_withBLtime.csv')
# filter_path = os.path.join(dataout_dir,'part_data_all_2.csv')
filter_info = pd.read_csv(filter_path)

offline_inf = pd.ExcelFile(r'/msegalrolab/luzhang/ORACLES/filterID_oracles.xlsx')
# offline_inf = pd.ExcelFile(r'/Users/zhanglu/Postdoc/TAU/Research/ORACLES/filterID_oracles.xlsx')
filter2017 = offline_inf.parse('2017')
filter2018 = offline_inf.parse('2018')

# a = 7.4
# b = 0
# c = -0.005

def read_data(i):
    print('Loading data ...')

    full_path_gituhsas = os.path.join(data_dir,filter2018.Filter[i],'newpsd.csv')
    gituhsas = pd.read_csv(full_path_gituhsas, index_col = 0, parse_dates=['time'])
    gituhsas[gituhsas<0]=0
    dp = gituhsas.columns.astype('float').values
    y = ((dp+22.1)/3.63)**(1/0.756) 
#     y = ((dp+22.1+a)/(3.63+b))**(1/(0.756+c)) 
    gituhsas.columns = y 
    full_path_uhsas = os.path.join(data_dir,filter2018.Filter[i],'psd.csv')
    uhsas = pd.read_csv(full_path_uhsas, index_col = 0, parse_dates=['time'])
    uhsas[uhsas<0]=0
    full_path_pcasp = os.path.join(data_dir,filter2018.Filter[i],'pcasp.csv')
    pcasp = pd.read_csv(full_path_pcasp, index_col = 0, parse_dates=['time'])
    pcasp[pcasp<0]=0
    dp_pcasp = pcasp.columns.astype('float').values
    span = np.diff(dp_pcasp)/2
    dp_min = -span+dp_pcasp[:-1]
    dp_min = np.concatenate([dp_min,[dp_pcasp[-1]-span[-1]]])
    dp_max = span+dp_pcasp[:-1]
    dp_max = np.concatenate([dp_max,[dp_pcasp[-1]+span[-1]]])
    dlogdp = np.log10(dp_max)-np.log10(dp_min)
    dpcasp = pcasp.div(dlogdp)
    
    full_path_ang = os.path.join(data_dir,filter2018.Filter[i],'Angstrom_exponent.csv')
    ang = pd.read_csv(full_path_ang, index_col = 0, parse_dates=['time'])

    full_path_neph = os.path.join(data_dir,filter2018.Filter[i], 'Neph.csv')
    neph = pd.read_csv(full_path_neph,index_col = 0, parse_dates=['time'])
    full_path_m903 = os.path.join(data_dir,filter2018.Filter[i], 'M903.csv')
    m903 = pd.read_csv(full_path_m903,index_col = 0, parse_dates=['time'])
    full_path_sp2 = os.path.join(data_dir,filter2018.Filter[i], 'sp2.csv')
    sp2 = pd.read_csv(full_path_sp2, index_col = 0, parse_dates=['time'])
    full_path_BCdis = os.path.join(data_dir,filter2018.Filter[i], 'BC_dis_timeseries.xlsx')
    BC_dis = pd.read_excel(full_path_BCdis,parse_dates=['Start_UTC','End_UTC','Mid_UTC'],index_col=1)
    ams_path = os.path.join(data_dir,filter2018.Filter[i],'ams.csv')
    ams = pd.read_csv(ams_path, index_col = 0, parse_dates=['time'])
    ams[ams<0]=np.nan
#     full_path_absorp_spec = full_path_absorp_spec = os.path.join(data_dir,filter2018.Filter[i],'PSAP_Spec.csv')  ##Specific rear
    full_path_absorp_spec = os.path.join(data_dir,filter2018.Filter[i],'PSAP_Avg.csv')  ##Specific rear
    absorp = pd.read_csv(full_path_absorp_spec, index_col = 0, parse_dates=['time'])
    full_path_aps = os.path.join(data_dir,filter2018.Filter[i],'aps.csv')
    aps = pd.read_csv(full_path_aps, index_col= 0, parse_dates=['time'])
    dp_aps = aps.columns.astype('float').values
    dp_aps_ve = dp_aps/1.5**.5 #1.5 was the average of RF10
    aps_ve = aps.copy()
    aps_ve.columns = dp_aps_ve    
    aps = aps_ve
    
    full_path_Coating = os.path.join(data_dir,filter2018.Filter[i],'2Dcoating.txt')
    coating = pd.read_csv(full_path_Coating,sep='\t') if os.path.exists(full_path_Coating) else np.nan
    
    full_path_coef = os.path.join(dataout_dir,filter2018.Filter[i],'coef_mea_interpolate_absavg.csv')
    coef_interp = pd.read_csv(full_path_coef, index_col= 0, parse_dates=['time'])
    coef_interp[coef_interp<1.5]=np.nan
#     full_path_apsTSP = os.path.join(data_dir,filter2018.Filter[i],'aps_TSP.csv')    
#     apsTSP = pd.read_csv(full_path_apsTSP, index_col= 0, parse_dates=['time'])
    day = filter2018.Date[i].day
#     print(day)
#     """in case mutiple episode in one filter"""
#     minute1 = filter2018.iloc[i,3].minute if type(filter2018.iloc[i,3]) != str else int(filter2018.iloc[i,3].split(',')[0].split(':')[1])
#     minute2 = filter2018.iloc[i,4].minute if type(filter2018.iloc[i,4]) != str else int(filter2018.iloc[i,4].split(',')[0].split(':')[1])

#     for fls in os.listdir(BC_dir):
#         if str("{0:0=2d}.txt".format(day)) in fls:      
#             print(fls)
#             BC_core = pd.read_csv('{}/{}'.format(BC_dir,fls),sep='\t')
#     for ii,BCi in enumerate(BC_core.columns[1:]):
#         a1 = int(BCi.split('_')[0][-4:-2])
#         a2 = int(BCi.split('_')[1][-4:-2])
#         if a1 == minute1 or a2 == minute2:
# #             print('Finished!')
#             BCcoredis = BC_core.iloc[:,ii+1:]
    print('Data loaded!')
    return(gituhsas,ang,aps,pcasp,dpcasp,neph,m903,sp2,BC_dis,coating,ams,absorp,coef_interp)

def vol_inorg(ams):
    # density
    Rho_NH4NO3 = 1.72                                       # From lide, 2001, Fierz et al., 2010
    Rho_NH42SO4 = 1.77                                      # From Lide, 2001, Fierz et al., 2010
    Rho_Na2SO4 = 2.66                                      # Wiki
    Rho_NH4HSO4 = 1.78                                      # From Fierz et al., 2010
    Rho_NH4Cl = 1.53
    
    n = ams.NH4/18
    s = ams.SO4/96
    n_mNO3 = n - ams.NO3/62 - ams.Chl/35.5
    n_NH4NO3 = np.zeros(len(n))
    n_NH42SO4 = np.zeros(len(n))
    n_NH4HSO4 = np.zeros(len(n))
    n_NH4Cl = np.zeros(len(n))
    n_residual_N = np.zeros(len(n))
    n_residual_S = np.zeros(len(n))
    for i in range(len(n)):
        if n_mNO3[i] - 2*s[i]>0:
            n_NH4NO3[i] = ams.NO3[i]/62  
            n_NH4Cl[i] = ams.Chl[i]/35.5
            n_NH42SO4[i] = s[i]
            n_residual_N[i] = n_mNO3[i] - 2*s[i]
    #                 print('2')
        elif n_mNO3[i]-s[i]>0:
            n_NH4NO3[i] = ams.NO3[i]/62
            n_NH4Cl[i] = ams.Chl[i]/35.5
            n_NH42SO4[i] = n_mNO3[i] - s[i] 
            n_NH4HSO4[i] = 2*s[i] - n_mNO3[i] 
    #                 n_residual_N[i] = n_mNO3[i] - 2*n_NH42SO4[i] - n_NH4HSO4[i]
    #                 print('3')
        elif n_mNO3[i]>0:
            n_NH4NO3[i] = ams.NO3[i]/62
            n_NH4Cl[i] = ams.Chl[i]/35.5
            n_NH4HSO4[i] = n_mNO3[i]
            n_residual_S[i] = s[i] - n_mNO3[i]
        else:
            n_NH4NO3[i] = ams.NH4[i]/18
            n_NH4Cl[i] = ams.Chl[i]/35.5

    v_NH42SO4 = (n_NH42SO4*132)/Rho_NH42SO4
    v_Na2SO4 = (n_residual_S*142)/Rho_Na2SO4
    v_NH4NO3 = (n_NH4NO3*80)/Rho_NH4NO3
    v_NH4HSO4 = (n_NH4HSO4*115)/Rho_NH4HSO4
    v_NH4Cl = n_NH4Cl*53.5/Rho_NH4Cl
       
    return(v_NH42SO4,v_Na2SO4,v_NH4NO3,v_NH4HSO4,v_NH4Cl)

gituhsas,ang,aps,pcasp,dpcasp,neph,m903,sp2,BC_dis,coating,ams,absorp,coef_interp = read_data(fileid)
ams = ams.resample('10S').mean()
sp2 = sp2.resample('10S').mean()
gituhsas = gituhsas.resample('10S').mean()
aps = aps.resample('10S').mean()

# fpathRI = "{}/{}/RI_BrC Copy.csv".format(dataout_dir,filter2018.Filter[fileid])
fpathRI = "{}/{}/RI_BrC_vn_absAvg.csv".format(dataout_dir,filter2018.Filter[fileid])
RI = pd.read_csv(fpathRI)

d_BC_info = pd.read_csv('/msegalrolab/luzhang/ORACLES/BC_dis_processed/D.csv',index_col=0)
dlogdbc = np.log10(d_BC_info.loc['Dmax']/d_BC_info.loc['Dmin'])

dp_gituhsas = gituhsas.columns.astype('float')
bin_path = os.path.join(data_dir,'GITUHSAS_bin.csv')
gituhsas_bin = pd.read_csv(bin_path)
gituhsas_bin.GITUHSASdiameters = ((gituhsas_bin.GITUHSASdiameters+22.1)/3.63)**(1/0.756)
# gituhsas_bin.GITUHSASdiameters = ((gituhsas_bin.GITUHSASdiameters+22.1+a)/(3.63+b))**(1/(0.756+c))

bins = gituhsas_bin[::2].append(gituhsas_bin.iloc[-1,:])
bins = bins.GITUHSASdiameters
dlogDp = np.diff(np.log10(bins))


###### read BCNSD and nonBCNSD ######
filepath_PNSD = os.path.join(dataout_dir,filter2018.Filter[fileid],'BCNSD.csv')
BCNSD = pd.read_csv(filepath_PNSD,index_col = 0, parse_dates=['End_UTC'])
filepath_PNSD_log = os.path.join(dataout_dir,filter2018.Filter[fileid],'BCNSD_log.csv')
BCNSD_log = pd.read_csv(filepath_PNSD_log,index_col = 0, parse_dates=['End_UTC'])
filepath_PNSD_log = os.path.join(dataout_dir,filter2018.Filter[fileid],'nonBCNSD_log.csv')
nonBCNSD_log = pd.read_csv(filepath_PNSD_log,index_col = 0, parse_dates=['time'])

# dCore is np.arange(10,1005,5)
dCore = 1000*BC_dis.columns[2:].astype('float')
coating = coating.fillna(0)
coating_idx = [(abs(coating.rBCCore_dia_nm-dCorei)).idxmin() for dCorei in dCore]
coating_new = coating.iloc[coating_idx,:]
coating_new = coating_new.reset_index(drop=True)

dCore = BC_dis.columns[2:].astype('float')*1000
vBCi = np.pi*dCore**3/6
vcoatedBCi = np.pi*dp_gituhsas**3/6
vBC = (BC_dis.iloc[:,2:]*dlogdbc.values*vBCi).sum(axis=1)
vcoatedBC = (BCNSD.mul(vcoatedBCi)).sum(axis=1)
vcoating = vcoatedBC-vBC
f1 = (vBC/vcoatedBC).values

# coating_new
Dp_BC = coating_new.WeightedCoatingThickness_nm*2+dCore
Dp3 = np.sum(Dp_BC.values**3*BC_dis.iloc[:,2:],axis=1)
Dc3 = np.sum((1*dCore)**3*BC_dis.iloc[:,2:],axis=1)
bulk_ratio = (Dp3/Dc3)**(1/3)
print(bulk_ratio.mean())

# Refractive index
# mCore = 1.95+0.79j  # Zanatta et al., 2018 Bond and Bergstrom, 2006
m_NH4Cl= {'450': 1.642, '550': 1.642, '650': 1.642, '700': 1.642}
m_NH4NO3 = {'450': 1.559, '550': 1.556, '650': 1.553, '700': 1.553}   # From Fierz et al., 2010
# m_NH42SO4 = {'450': 1.536, '550': 1.530, '700': 1.524}  # From Fierz et al., 2010
# m_NH4HSO4 = {'450': 1.462, '550': 1.442, '700': 1.43}  # From Fierz et al., 2010
m_NH42SO4 = {'450': 1.46279, '550': 1.45163, '650': 1.44522, '700': 1.44522}  # From Cotterell et al., 2017
m_Na2SO4 = {'450': 1.471, '550': 1.434, '650': 1.413,'700': 1.413}  # From Cotterell et al., 2017
m_K2SO4 = {'450': 1.495, '550': 1.495, '650': 1.495,'700': 1.495}  # https://pubchem.ncbi.nlm.nih.gov/compound/Potassium-sulfate#section=Stability-Shelf-Life
m_NH4HSO4 = {'450': 1.46183, '550': 1.44207, '650': 1.4307,'700': 1.4307}   # From Cotterell et al., 2017 http://eodg.atm.ox.ac.uk/ARIA/data?Salts/Sodium_Sulphate/60%25_(Cotterell_et_al._2017)/Na2SO4_60_Cotterell_2017.ri
# m_BC = {'450': 1.75 + 0.46j, '550': 1.75 + 0.44j, '650': 1.75 + 0.43j, '700': 1.75 + 0.43j}       # From Fierz et al., 2010
m_H2SO4 = {'450': 1.427, '550': 1.427, '650': 1.427, '700': 1.427}
m_BC = {'450': 2.26+1.26j, '550': 2.26+1.26j, '650': 2.26+1.26j, '700': 2.26+1.26j}  
# m_BC = {'450': 1.9 + 0.4j, '550': 1.9 + 0.4j, '700': 1.9 + 0.4j}       # assumption
# m_BC = {'450': 1.8 + 0.74j, '550': 1.8 + 0.74j, '700': 1.8 + 0.74j}           # From Feng et al., 2013, not suitable here
# m_Org = {'450': 1.48, '550': 1.48, '650': 1.48}                             # From Fierz et al., 2010
m_Org = {'450': 1.65, '550': 1.65, '650': 1.65,'700': 1.65}                               # Feng et al., 2013 only real part
m_BrC_m = {'450': 1.65 +0.02j, '550': 1.65 + 0.003j, '650': 1.65 + 0.0003j,'700': 1.65 + 0.0003j}   # Feng et al., 2013
m_BrC_s = {'450': 1.65 +0.063j, '550': 1.65 + 0.03j, '650': 1.65 + 0.005j,'700': 1.65 + 0.005j}    # Feng et al., 2013
m_H20 = {'450': 1.337, '550': 1.333, '650': 1.331, '700': 1.331} #https://refractiveindex.info/?shelf=main&book=H2O&page=Hale
# density
Rho_NH4NO3 = 1.72                                       # From lide, 2001, Fierz et al., 2010
Rho_NH42SO4 = 1.77                                      # From Lide, 2001, Fierz et al., 2010
Rho_Na2SO4 = 2.66                                      # Wiki
Rho_K2SO4 = 2.66                                      # Wiki
Rho_NH4HSO4 = 1.78                                      # From Fierz et al., 2010
Rho_NH4Cl = 1.53
# Rho_BC = 1.7                                            # From Fierz et al., 2010
Rho_BC = 1.8
# Rho_BC = 1.35                                           # From Bond 2006
Rho_Org = 1.4                                           # From Lide, 2001, Fierz et al., 2010
# Rho_Org = 1.77                                          # Cross et al., 2007
# Rho_Inorg = 1.2                                         # Cross et al., 2007
Rho_BrC = 1.569                                         # Feng et al., 2013
Rho_H2SO4 = 1.83                                        # From Sjogren et al., 2008
# growth factor g(85%)
g_NH4NO3 = 1.59                                         # From Sjogren et al., 2008
g_NH42SO4 = 1.56                                        # From Sjogren et al., 2008
g_NH4HSO4 = 1.62                                        # From Sjogren et al., 2008
g_BC = 1.0                                              # From Sjogren et al., 2008
g_Org = 1.2                                             # From Sjogren et al., 2008
g_H2SO4 = 1.88                                          # From Sjogren et al., 2008
g_BrC = 1.0
# kappa to be determined
k_NH4NO3 = 1.59                                       # From Fierz et al., 2010
k_NH42SO4 = 1.56                                      # From Fierz et al., 2010
k_NH4HSO4 = 1.62                                      # From Fierz et al., 2010
k_BC = 1.0                                            # From Fierz et al., 2010
k_Org = 1.2                                           # From Fierz et al., 2010
k_BrC = 1.0


def Ion_balance(wavelength,ams,BC_dis,m_BrC_test,mCore):
    print('Testing ion balance ...')
    
    n = ams.NH4/18
    s = ams.SO4/96
    n_mNO3 = n - ams.NO3/62 - ams.Chl/35.5
    n_NH4NO3 = np.zeros(len(n))
    n_NH42SO4 = np.zeros(len(n))
    n_NH4HSO4 = np.zeros(len(n))
    n_NH4Cl = np.zeros(len(n))
    n_residual_N = np.zeros(len(n))
    n_residual_S = np.zeros(len(n))
    for i in range(len(n)):
        if n_mNO3[i] - 2*s[i]>0:
            n_NH4NO3[i] = ams.NO3[i]/62  
            n_NH4Cl[i] = ams.Chl[i]/35.5
            n_NH42SO4[i] = s[i]
            n_residual_N[i] = n_mNO3[i] - 2*s[i]
    #                 print('2')
        elif n_mNO3[i]-s[i]>0:
            n_NH4NO3[i] = ams.NO3[i]/62
            n_NH4Cl[i] = ams.Chl[i]/35.5
            n_NH42SO4[i] = n_mNO3[i] - s[i] 
            n_NH4HSO4[i] = 2*s[i] - n_mNO3[i] 
    #                 n_residual_N[i] = n_mNO3[i] - 2*n_NH42SO4[i] - n_NH4HSO4[i]
    #                 print('3')
        elif n_mNO3[i]>0:
            n_NH4NO3[i] = ams.NO3[i]/62
            n_NH4Cl[i] = ams.Chl[i]/35.5
            n_NH4HSO4[i] = n_mNO3[i]
            n_residual_S[i] = s[i] - n_mNO3[i]
        else:
            n_NH4NO3[i] = ams.NH4[i]/18
            n_NH4Cl[i] = ams.Chl[i]/35.5

    df = pd.DataFrame({
        '$NH_4NO_3$':n_NH4NO3,
        '$NH_4Cl$':n_NH4Cl,
        '$(NH_4)_2SO_4$': 2*n_NH42SO4,
        '$NH_4HSO_4$': n_NH4HSO4,
        '$residual\, NH_4^+$': n_residual_N,
        '$residual\, SO_4^{2-}$': n_residual_S})
    df.index = ams.index
    df_N = pd.DataFrame({
        '$NH_4NO_3$':n_NH4NO3/n*100,
        '$NH_4Cl$':n_NH4Cl/n*100,
        '$(NH_4)_2SO_4$': 2*n_NH42SO4/n*100,
        '$NH_4HSO_4$': n_NH4HSO4/n*100,
        '$residual\, NH_4^+$': n_residual_N/n*100})
    df_S = pd.DataFrame({
        '$(NH_4)_2SO_4$': n_NH42SO4/s*100,
        '$NH_4HSO_4$': n_NH4HSO4/s*100,
        '$residual\, SO_4^{2-}$': n_residual_S/s*100})
    df_N.index = ams.index
    df_S.index = ams.index
    
    # Assuming chemical form and calculating respective volumn
    v_NH42SO4 = (n_NH42SO4*132)/Rho_NH42SO4
    # v_Na2SO4 = (n_residual_S*142)/Rho_Na2SO4
    v_K2SO4 = (n_residual_S*142)/Rho_K2SO4
    v_NH4NO3 = (n_NH4NO3*80)/Rho_NH4NO3
    v_NH4HSO4 = (n_NH4HSO4*115)/Rho_NH4HSO4
    v_NH4Cl = n_NH4Cl*53.5/Rho_NH4Cl
    v_Org = ams.ORG/Rho_Org
    v_OrgexceptBrC = ams.ORG*(1-fBrC)/Rho_Org
    v_BrC = ams.ORG*fBrC/Rho_BrC

    d_BC_info = pd.read_csv('/msegalrolab/luzhang/ORACLES/BC_dis_processed/D.csv',index_col=0)
    dlogdbc = np.log10(d_BC_info.loc['Dmax']/d_BC_info.loc['Dmin'])
    dBC = BC_dis.columns[2:].astype('float').values
    vBC_s = np.pi*dBC**3/6
    mBC = np.asarray([np.sum(vBC_s*BC_dis.iloc[i,2:]*dlogdbc.values)*1.8/100 for i in range(len(BC_dis.index))])
    v_BC = np.asarray([np.sum(vBC_s*BC_dis.iloc[i,2:]*dlogdbc.values)/100 for i in range(len(BC_dis.index))])


    
    Rho = (Rho_NH42SO4*v_NH42SO4\
          +Rho_Na2SO4*v_K2SO4\
          +Rho_NH4NO3*v_NH4NO3\
          +Rho_NH4HSO4*v_NH4HSO4\
          +Rho_NH4Cl*v_NH4Cl\
          +Rho_BC*v_BC\
          +Rho_Org*v_Org)/(v_NH42SO4+v_K2SO4+v_NH4NO3+v_NH4HSO4+v_NH4Cl+v_Org+v_BC)
    Rho_coating = (Rho_NH42SO4*v_NH42SO4\
          +Rho_Na2SO4*v_K2SO4\
          +Rho_NH4NO3*v_NH4NO3\
          +Rho_NH4HSO4*v_NH4HSO4\
          +Rho_NH4Cl*v_NH4Cl\
          +Rho_Org*v_Org)/(v_NH42SO4+v_K2SO4+v_NH4NO3+v_NH4HSO4+v_NH4Cl+v_Org)


    m = []
    ## OC:OM = 1:1.4 ##
    ## BrC:OC = 92% for biofuel combustion and biomass burning sources ##
    ## BrC:OM = 66% (92%/1.4) ##
    ## Another 34% of OM is assumed to be non-absorbing ##
    def calc_m(x):
        return (m_NH42SO4[wavelength]*v_NH42SO4
                +m_Na2SO4[wavelength]*v_K2SO4
                +m_NH4NO3[wavelength]*v_NH4NO3
                +m_NH4HSO4[wavelength]*v_NH4HSO4
                +m_NH4Cl[wavelength]*v_NH4Cl
                +mCore*v_BC
                +m_Org[wavelength]*v_OrgexceptBrC+x*v_BrC)/(v_NH42SO4+v_K2SO4+v_NH4NO3+v_NH4HSO4+v_NH4Cl+v_OrgexceptBrC+v_BrC+v_BC)
    def calc_mShell(x):
        return (m_NH42SO4[wavelength]*v_NH42SO4
                +m_Na2SO4[wavelength]*v_K2SO4
                +m_NH4NO3[wavelength]*v_NH4NO3
                +m_NH4HSO4[wavelength]*v_NH4HSO4
                +m_NH4Cl[wavelength]*v_NH4Cl
                +m_Org[wavelength]*v_OrgexceptBrC+x*v_BrC)/(v_NH42SO4+v_K2SO4+v_NH4NO3+v_NH4HSO4+v_NH4Cl+v_OrgexceptBrC+v_BrC)


    mShell1 = (list(map(calc_mShell, m_BrC_test)))
    mShell = pd.DataFrame(mShell1).T
    m1 = (list(map(calc_m, m_BrC_test)))
    m = pd.DataFrame(m1).T

    k = (k_NH42SO4*v_NH42SO4
        +k_NH4NO3*v_NH4NO3
        +k_Org*v_Org
        +k_NH4HSO4*v_NH4HSO4
        +k_BC*v_BC)/(v_NH42SO4+v_NH4NO3+v_Org+v_NH4HSO4+v_BC)
    
#     m.to_csv('{}{}/m_{}.csv'.format(dataout_dir,filter2018.Filter[fileid],wavelength))
#     mShell.to_csv('{}{}/mShell_{}.csv'.format(dataout_dir,filter2018.Filter[fileid],wavelength))
#     Rho.to_csv('{}{}/Rho.csv'.format(dataout_dir,filter2018.Filter[fileid]))
#     print('Refractive index saved!')
    return(Rho,Rho_coating,m,mShell,k)



def BG_best_mCore(m_BC_test,m_Org,wavelength,fileid):
    mea_coef = coef_interp.loc[:,['Ext660','Sca660','Abs660']]
    mea_coef.columns = ['Ext','Sca','Abs']
    mBC_1 = pd.DataFrame(columns=['mBC1'],index = ams.index)
    mBC_2 = pd.DataFrame(columns=['mBC2'],index = ams.index)
    m_ORG_sca = pd.DataFrame(columns=['m_ORG'],index = ams.index)
    m_Org_test = np.arange(1.35,1.71,0.1)
    for i in range(len(ams.index)):
        print(ams.index[i])
        if ams.index[i] not in mea_coef.dropna().index or ams.index[i] not in ams.dropna().index:
            continue
        else:
            ############## refractive index ############
            wavelength_m = '650'
            v_NH42SO4,v_Na2SO4,v_NH4NO3,v_NH4HSO4,v_NH4Cl = vol_inorg(ams.iloc[[i],:])
            v_Org = ams.iloc[i,:].ORG/1.4
            vShell = v_Org+v_NH4Cl+v_NH4HSO4+v_NH4NO3+v_Na2SO4+v_NH42SO4
            mShell = (v_NH42SO4 * m_NH42SO4[wavelength_m]+
            v_Na2SO4 * m_Na2SO4[wavelength_m]+
            v_NH4NO3 * m_NH4NO3[wavelength_m]+
            v_NH4HSO4 * m_NH4HSO4[wavelength_m]+
            v_NH4Cl * m_NH4Cl[wavelength_m]+
            v_Org * m_Org)/vShell
            print('mShell',mShell)

            ############## RI ####################
            cal_coefi = pd.DataFrame()
            for j in m_BC_test:
                mCore = j
                perm1 = mCore**2
                perm2 = mShell**2
                bb = 3*f1[i]*(perm1-perm2)+2*perm2-perm1
                m_BG = ((bb+(bb**2+8*perm1*perm2)**.5)/4)**.5    
                mcoatedBC = pd.DataFrame(m_BG)
                mShell = pd.DataFrame(mShell)

                cal_coef_cs = models.Mie_aps(BCNSD_log.iloc[[i],:],mcoatedBC,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')   
                cal_coef_mie = models.Mie_aps(nonBCNSD_log.iloc[[i],:],mShell,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')   
                cal_coef_aps = models.Mie_aps(aps.iloc[[i],:],mShell,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')
                cal_coefij = cal_coef_cs+cal_coef_aps+cal_coef_mie
    #             print('cal_coefij',cal_coefij)
                cal_coefi = pd.concat([cal_coefi,cal_coefij])

            abs_cor = ((cal_coefi.Abs-mea_coef.iloc[i,:].Abs)/mea_coef.iloc[i,:].Abs)**2
            var = abs_cor
            po = np.argmin(abs(var))
            mBC_1.iloc[i,:] = m_BC_test[po] if pd.isna(var[0])==False else np.nan
            print('po1',po,'var',var[po],'var_abs',abs_cor[po],
                  'abs_cor',abs_cor,'mBC',m_BC_test[po])
    return(mBC_1)






def VM_best_mCore(m_BC_test,m_Org,wavelength,fileid):
    mea_coef = coef_interp.loc[:,['Ext660','Sca660','Abs660']]
    mea_coef.columns = ['Ext','Sca','Abs']
    mBC_1 = pd.DataFrame(columns=['mBC1'],index = ams.index)
    mBC_2 = pd.DataFrame(columns=['mBC2'],index = ams.index)
    m_ORG_sca = pd.DataFrame(columns=['m_ORG'],index = ams.index)
    m_Org_test = np.arange(1.35,1.71,0.1)
    for i in range(len(ams.index)):
        print(ams.index[i])
        if ams.index[i] not in mea_coef.dropna().index or ams.index[i] not in ams.dropna().index:
            continue
        else:
            ############## refractive index ############
            wavelength_m = '650'
            v_NH42SO4,v_Na2SO4,v_NH4NO3,v_NH4HSO4,v_NH4Cl = vol_inorg(ams.iloc[[i],:])
            v_Org = ams.iloc[i,:].ORG/1.4
            vShell = v_Org+v_NH4Cl+v_NH4HSO4+v_NH4NO3+v_Na2SO4+v_NH42SO4
            mShell = (v_NH42SO4 * m_NH42SO4[wavelength_m]+
            v_Na2SO4 * m_Na2SO4[wavelength_m]+
            v_NH4NO3 * m_NH4NO3[wavelength_m]+
            v_NH4HSO4 * m_NH4HSO4[wavelength_m]+
            v_NH4Cl * m_NH4Cl[wavelength_m]+
            v_Org * m_Org)/vShell


            ############## RI ####################
            cal_coefi = pd.DataFrame()
            for j in m_BC_test:
                mCore = j
                mcoatedBC = (mCore*vBC[i]+vcoating[i]*(mShell))/vcoatedBC[i]
                mcoatedBC = pd.DataFrame(mcoatedBC)
                mShell = pd.DataFrame(mShell)

                cal_coef_cs = models.Mie_aps(BCNSD_log.iloc[[i],:],mcoatedBC,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')   
                cal_coef_mie = models.Mie_aps(nonBCNSD_log.iloc[[i],:],mShell,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')   
                cal_coef_aps = models.Mie_aps(aps.iloc[[i],:],mShell,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')
                cal_coefij = cal_coef_cs+cal_coef_aps+cal_coef_mie
    #             print('cal_coefij',cal_coefij)
                cal_coefi = pd.concat([cal_coefi,cal_coefij])

            abs_cor = ((cal_coefi.Abs-mea_coef.iloc[i,:].Abs)/mea_coef.iloc[i,:].Abs)**2
            var = abs_cor
            po = np.argmin(abs(var))
            mBC_1.iloc[i,:] = m_BC_test[po] if pd.isna(var[0])==False else np.nan
            print('po1',po,'var',var[po],'var_abs',abs_cor[po],
                  'abs_cor',abs_cor,'mBC',m_BC_test[po])
    return(mBC_1)







def MG_best_mCore(m_BC_test,m_Org,wavelength,fileid):
    mea_coef = coef_interp.loc[:,['Ext660','Sca660','Abs660']]
    mea_coef.columns = ['Ext','Sca','Abs']
    mBC_1 = pd.DataFrame(columns=['mBC1'],index = ams.index)
    mBC_2 = pd.DataFrame(columns=['mBC2'],index = ams.index)
    m_ORG_sca = pd.DataFrame(columns=['m_ORG'],index = ams.index)
    m_Org_test = np.arange(1.35,1.71,0.1)
    for i in range(len(ams.index)):
        print(ams.index[i])
        if ams.index[i] not in mea_coef.dropna().index or ams.index[i] not in ams.dropna().index:
            continue
        else:
            ############## refractive index ############
            wavelength_m = '650'
            v_NH42SO4,v_Na2SO4,v_NH4NO3,v_NH4HSO4,v_NH4Cl = vol_inorg(ams.iloc[[i],:])
            v_Org = ams.iloc[i,:].ORG/1.4
            vShell = v_Org+v_NH4Cl+v_NH4HSO4+v_NH4NO3+v_Na2SO4+v_NH42SO4
            mShell = (v_NH42SO4 * m_NH42SO4[wavelength_m]+
            v_Na2SO4 * m_Na2SO4[wavelength_m]+
            v_NH4NO3 * m_NH4NO3[wavelength_m]+
            v_NH4HSO4 * m_NH4HSO4[wavelength_m]+
            v_NH4Cl * m_NH4Cl[wavelength_m]+
            v_Org * m_Org)/vShell

            ############## RI ####################
            cal_coefi = pd.DataFrame()
            for j in m_BC_test:
                mCore = j
                perm1 = mCore**2
                perm2 = mShell**2
                m_MG = ((perm1+2*perm2+2*f1[i]*(perm1-perm2))/(perm1+2*perm2-f1[i]*(perm1-perm2))*perm2)**0.5
                mcoatedBC = pd.DataFrame(m_MG)
                mShell = pd.DataFrame(mShell)

                cal_coef_cs = models.Mie_aps(BCNSD_log.iloc[[i],:],mcoatedBC,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')   
                cal_coef_mie = models.Mie_aps(nonBCNSD_log.iloc[[i],:],mShell,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')   
                cal_coef_aps = models.Mie_aps(aps.iloc[[i],:],mShell,wavelength,ams.iloc[[i],:],fileid,'testRRI.csv')
                cal_coefij = cal_coef_cs+cal_coef_aps+cal_coef_mie
    #             print('cal_coefij',cal_coefij)
                cal_coefi = pd.concat([cal_coefi,cal_coefij])

            abs_cor = ((cal_coefi.Abs-mea_coef.iloc[i,:].Abs)/mea_coef.iloc[i,:].Abs)**2
            var = abs_cor
            po = np.argmin(abs(var))
            mBC_1.iloc[i,:] = m_BC_test[po] if pd.isna(var[0])==False else np.nan
            print('po1',po,'var',var[po],'var_abs',abs_cor[po],
                  'abs_cor',abs_cor,'mBC',m_BC_test[po])
    return(mBC_1)


def CS_best_mCore(mCore,m_Org,wavelength,fileid):
    mea_coef = coef_interp.loc[:,['Ext660','Sca660','Abs660']]
    mea_coef.columns = ['Ext','Sca','Abs']
    mBC_1 = pd.DataFrame(columns=['mBC1'],index = ams.index)
    mBC_2 = pd.DataFrame(columns=['mBC2'],index = ams.index)
    
    Abs = np.zeros(len(BC_dis.index))
    Sca = np.zeros(len(BC_dis.index))
    Ext = np.zeros(len(BC_dis.index))
    dShell = np.arange(0,222,2)
    _length = np.size(dShell)
    _length_B = np.size(coating_new.rBCCore_dia_nm)

    for k in range(len(BC_dis.index)):
        print(ams.index[k])
        if ams.index[k] not in mea_coef.dropna().index or ams.index[k] not in ams.dropna().index:
            continue
        else:
            Q_ext = np.zeros(_length)
            Q_sca = np.zeros(_length)
            Q_abs = np.zeros(_length)
            Q_pr = np.zeros(_length)
            Q_back = np.zeros(_length)
            Q_ratio = np.zeros(_length)
            g = np.zeros(_length)
            aSDn = np.zeros(_length)
            Bext= np.zeros((_length_B))
            Bsca= np.zeros((_length_B))
            Babs= np.zeros((_length_B))
            Bback= np.zeros((_length_B))
            Bratio= np.zeros((_length_B))
            bigG= np.zeros((_length_B))
            Bpr = np.zeros((_length_B))

            cal_coefi = pd.DataFrame()
            for mCore in m_BC_test:
                wavelength_m = '650'
                v_NH42SO4,v_Na2SO4,v_NH4NO3,v_NH4HSO4,v_NH4Cl = vol_inorg(ams.iloc[[k],:])
                v_Org = ams.iloc[k,:].ORG/1.4
                vShell = v_Org+v_NH4Cl+v_NH4HSO4+v_NH4NO3+v_Na2SO4+v_NH42SO4
                mShell = (v_NH42SO4 * m_NH42SO4[wavelength_m]+
                v_Na2SO4 * m_Na2SO4[wavelength_m]+
                v_NH4NO3 * m_NH4NO3[wavelength_m]+
                v_NH4HSO4 * m_NH4HSO4[wavelength_m]+
                v_NH4Cl * m_NH4Cl[wavelength_m]+
                v_Org * m_Org)/vShell
                mShell = pd.DataFrame(mShell)
                for i in coating_new.index:
                    Dc = coating_new.rBCCore_dia_nm[i]
                    coating_Dc = coating_new.iloc[i,3:]
                    coating_Dc_sum = np.sum(coating_Dc)
                    for j in range(len(dShell)):  
                        Pj = coating_Dc[j]/coating_Dc_sum if coating_Dc_sum>0 else 0
                        Dp = Dc+2*dShell[j]
                        aSDn[j] = np.pi*((Dp/2)**2)*BC_dis.iloc[k,2+i]*Pj*(1e-6)
                        if aSDn[j]>0:
                            Q_ext[j],Q_sca[j],Q_abs[j],g[j],Q_pr[j],Q_back[j],Q_ratio[j]=ps.MieQCoreShell(mCore, mShell.iloc[0,0],wavelength,Dc,Dp)    
                                                                        # Wavelength and particle diameters: nanometers, 
                                                                        # efficiencies: unitless, cross-sections: nm2, 
                                                                        # coefficients: Mm-1, 
                                                                        # size distribution concentration: cm-3
                        else:
                            Q_ext[j] = 0
                            Q_sca[j] = 0
                            Q_abs[j] = 0
                    Bext[i] = np.sum(Q_ext*aSDn)
                    Bsca[i] = np.sum(Q_sca*aSDn)
                    Babs[i] = Bext[i]-Bsca[i]

                Ext[k] = np.nansum(Bext*dlogdbc.values)
                Abs[k] = np.nansum(Babs*dlogdbc.values)
                Sca[k] = Ext[k]-Abs[k]
                cal_coef_cs = pd.DataFrame([Ext[k],Sca[k],Abs[k]]).T
                cal_coef_cs.index = ams.index[[k]]
                cal_coef_cs.columns = ['Ext','Sca','Abs']
#                 print('cal_coef_cs',cal_coef_cs)
                cal_coef_mie = models.Mie_aps(nonBCNSD_log.iloc[[k],:],mShell,wavelength,ams.iloc[[k],:],fileid,'testRRI.csv')   
                cal_coef_aps = models.Mie_aps(aps.iloc[[k],:],mShell,wavelength,ams.iloc[[k],:],fileid,'testRRI.csv')
                cal_coefij = cal_coef_cs+cal_coef_aps+cal_coef_mie
#                 print(cal_coefij)
                cal_coefi = pd.concat([cal_coefi,cal_coefij])
            print('cal_coefi',cal_coefi)
            abs_cor = ((cal_coefi.Abs-mea_coef.iloc[k,:].Abs)/mea_coef.iloc[k,:].Abs)**2
            var = abs_cor
            po = np.argmin(abs(var))
            mBC_1.iloc[k,:] = m_BC_test[po] if pd.isna(var[0])==False else np.nan
            print('po1',po,'var',var[po],'var_abs',abs_cor[po],
                  'abs_cor',abs_cor,'mBC',m_BC_test[po])
    return(mBC_1)




wavelength = 660
RRI = 1.6 # RRI of OA, set to 1.6, Liu DT, 2021, est

# imag_mCore = np.arange(0.8,1.31,.01)
# m_BC_test = [1.95+keBCi*1j for keBCi in imag_mCore]

# mBC_1 = BG_best_mCore(m_BC_test,RRI,wavelength,fileid)
# fpath = "{}/{}/{}/RI_absAvg_meBC_each_abs.csv".format(dataout_dir,filter2018.Filter[fileid],'BG')
# outdata = pd.DataFrame(mBC_1)
# outdata.to_csv(fpath)

# mBC_1 = MG_best_mCore(m_BC_test,RRI,wavelength,fileid)
# fpath = "{}/{}/{}/RI_absAvg_meBC_each_abs.csv".format(dataout_dir,filter2018.Filter[fileid],'MG')
# outdata = pd.DataFrame(mBC_1)
# outdata.to_csv(fpath)

# mBC_1 = VM_best_mCore(m_BC_test,RRI,wavelength,fileid)
# fpath = "{}/{}/{}/RI_absAvg_meBC_each_abs.csv".format(dataout_dir,filter2018.Filter[fileid],'VM')
# outdata = pd.DataFrame(mBC_1)
# outdata.to_csv(fpath)



imag_mCore = np.arange(1.1,1.81,.01)
m_BC_test = [1.95+keBCi*1j for keBCi in imag_mCore]
mBC_1 = CS_best_mCore(m_BC_test,RRI,wavelength,fileid)
fpath = "{}/{}/{}/RI_absAvg_meBC_each_abs.csv".format(dataout_dir,filter2018.Filter[fileid],'CS')
outdata = pd.DataFrame(mBC_1)
outdata.to_csv(fpath)






