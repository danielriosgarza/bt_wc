# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 10:57:56 2021

@author: u0139894
"""

import streamlit as stl
import matplotlib.pyplot as plt

from MonodCom_simulator import *




###################Experimental Data#######################

time_exp = np.array([0, 4, 6, 8, 10, 12, 14, 16, 48])
biomass_exp = np.array([0.0098,0.0440,0.1477,0.4025,0.7105,0.8433,0.9160,0.9275,0.8223])
glucose_exp = np.array([1.1276,1.1139,1.0581,0.9161,0.6495,0.3921,0.1578,0.0548,0.0157])
pyruvate_exp = np.array([0.8929,0.8754,0.8587,0.8506,0.7990,0.6386,0.3929,0.2121,0.0220])
succinate_exp = np.array([0.0526,0.0598,0.0889,0.1849,0.3446,0.5134,0.6298,0.6871,0.8417])
lactate_exp = np.array([0.1703,0.1701,0.1703,0.1753,0.2128,0.3164,0.3751,0.3950,0.3683])
acetate_exp = np.array([0.0664,0.0694,0.0961,0.2626,1.4017,1.7584,1.8904,1.9114,1.2832])

##########################################################    

fitted_data = np.array([ 1.50000e-03, 1.93870e+00, 8.00000e-03, 9.53010e+00, 7.55900e-01,
 1.13670e+00, 1.06065e+01, -4.74630e+00, -6.14500e-01, -2.01780e+00, 2.79443e+01, 8.68121e+01, 3.44627e+01, 2.34176e+01, 3.41022e+01, 0.00000e+00, 0.00000e+00, 1.13000e-02, 0.00000e+00, 1.49000e-02])

# fitted_data[0] = lag_param
# fitted_data[1] = mu_max
# fitted_data[2] = death_rate (dr)
# fitted_data[3] = weight growth rate term1 (glucose)
# fitted_data[4] = weight growth rate term2 (pyruvate)
# fitted_data[5] = glucose v
# fitted_data[6] = pyruvate v
# fitted_data[7] = lactate v
# fitted_data[8] = succinate v
# fitted_data[9] = acetate v
# fitted_data[10] = glucose k
# fitted_data[11] = pyruvate k
# fitted_data[12] = lactate k
# fitted_data[13] = succinate k
# fitted_data[14] = acetate k
# fitted_data[15] = glucose deg
# fitted_data[16] = pyruvate deg
# fitted_data[17] = lactate deg
# fitted_data[18] = succinate deg
# fitted_data[19] = acetate deg



name = 'BT'

lphase_param = fitted_data[0]
mu_max =  fitted_data[1]
dr = fitted_data[2]
w_BT = [fitted_data[3], fitted_data[4]] 
metabs_BT = ['Glucose', 'Pyruvate', 'Lactate', 'Succinate', 'Acetate']
metabs_BT_v = [fitted_data[5], fitted_data[6], fitted_data[7], fitted_data[8], fitted_data[9]]
metabs_BT_k = [fitted_data[10], fitted_data[11], fitted_data[12], fitted_data[13], fitted_data[14]]
growth_BT = [(w_BT[0], 'Glucose'), (w_BT[1], 'Pyruvate')]
feeding_BT = [(1,0), (0,1), (0,1), (1,1), (1,1)]
        
    

    
metabolome = ['Glucose', 'Pyruvate', 'Lactate', 'Succinate', 'Acetate']






######User control#######

stl.sidebar.write('### Parameters')

x_0_exp = biomass_exp[0]
x_0 = stl.sidebar.slider('Initial bacterial concentration (OD)', min_value = 0.0, max_value = 1.0, value = float(x_0_exp), step=0.001, format='%.4f')



stl.sidebar.write('*-----------------------------------------*\n\n\n')
stl.sidebar.write('Glucose')


glucose_s0_exp = glucose_exp[0]
glucose_s0 = stl.sidebar.slider('Initial glucose concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(glucose_s0_exp), step=0.25, format='%.2f')



glucose_feed = stl.sidebar.slider('Glucose feed concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(glucose_s0_exp), step=0.25, format='%.2f')


stl.sidebar.write('*-----------------------------------------*\n\n\n')
stl.sidebar.write('Pyruvate')

pyruvate_s0_exp = pyruvate_exp[0]
pyruvate_s0 = stl.sidebar.slider('Initial pyruvate concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(pyruvate_s0_exp), step=0.25, format='%.2f')


pyruvate_feed = stl.sidebar.slider('Pyruvate feed concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(pyruvate_s0_exp), step=0.25, format='%.2f')

stl.sidebar.write('*-----------------------------------------*\n\n\n')
stl.sidebar.write('Lactate')

lactate_s0_exp = lactate_exp[0]
lactate_s0 = stl.sidebar.slider('Initial lactate concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(lactate_s0_exp), step=0.25, format='%.2f')

lactate_feed = stl.sidebar.slider('Lactate feed concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(lactate_s0_exp), step=0.25, format='%.2f')


stl.sidebar.write('*-----------------------------------------*\n\n\n')
stl.sidebar.write('Succinate')

succinate_s0_exp = succinate_exp[0]
succinate_s0 = stl.sidebar.slider('Initial succinate concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(succinate_s0_exp), step=0.25, format='%.2f')

succinate_feed = stl.sidebar.slider('Succinate feed concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(succinate_s0_exp), step=0.25, format='%.2f')

stl.sidebar.write('*-----------------------------------------*\n\n\n')
stl.sidebar.write('Acetate')

acetate_s0_exp = acetate_exp[0]
acetate_s0 = stl.sidebar.slider('Initial acetate concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(acetate_s0_exp), step=0.25, format='%.2f')

acetate_feed = stl.sidebar.slider('Acetate feed concentration (g/L)', min_value = 0.0, max_value = 25.0, value = float(acetate_s0_exp), step=0.25, format='%.2f')



stl.sidebar.write('*-----------------------------------------*\n\n\n')
stl.sidebar.write('Culture conditions')

dparam = stl.sidebar.slider('Dilution parameter (1/h)', min_value=0.0, max_value = 1.0, value=0.0, step=0.005, format='%.4f')

stl.sidebar.write('the maximum dilution allowed by the AMBR is approx. 0.07 (considering a vol of 10 mL and inflow of 0.7 ml/h)')

tend = stl.sidebar.slider('Culture time (h)', min_value=5.0, max_value = 1000.0, value=float(time_exp[-1]), step=0.5, format='%.2f')


metabolome_c = [glucose_s0, pyruvate_s0, lactate_s0, succinate_s0, acetate_s0]
metabolome_c_exp = [glucose_s0_exp, pyruvate_s0_exp, lactate_s0_exp, succinate_s0_exp, acetate_s0_exp]
metabolome_f = [glucose_feed, pyruvate_feed, lactate_feed, succinate_feed, acetate_feed]
dilution = dparam



######################

metabolome_deg = [fitted_data[15], fitted_data[16], fitted_data[17], fitted_data[18], fitted_data[19]]

BT = Strain(name=name, x_0 = x_0, q_0 = lphase_param, mu=mu_max, dr=dr, metab = metabs_BT, metab_v = metabs_BT_v, metab_k = metabs_BT_k, growth_model = growth_BT, feeding = feeding_BT)

BT_exp = Strain(name=name, x_0 = x_0_exp, q_0 = lphase_param, mu=mu_max, dr=dr, metab = metabs_BT, metab_v = metabs_BT_v, metab_k = metabs_BT_k, growth_model = growth_BT, feeding = feeding_BT)

culture = Culture(strains = [BT], metabolome = metabolome, metabolome_c =metabolome_c, metab_deg = metabolome_deg, dilution=dilution, feed_c = metabolome_f)




culture.simulate(0, tend, nsteps=1000)

stl.markdown("<h1 style='text-align: center; color: red;'>BT (WC media)</h1>", unsafe_allow_html=True)
stl.line_chart(culture.community_dyn)
stl.line_chart(culture.environment_dyn)


@stl.cache
def exp():
    culture_exp = Culture(strains = [BT_exp], metabolome = metabolome, metabolome_c =metabolome_c_exp, metab_deg = metabolome_deg)
    
    
    culture_exp.simulate(time_exp[0], time_exp[-1])

    return culture_exp

stl.markdown("<h1 style='text-align: center; color: red;'>Parameter Fit</h1>", unsafe_allow_html=True)

culture_exp = exp()
fig1, ax1 = plt.subplots()
ax1.plot(time_exp, biomass_exp, 'o', label='OD600', color='b')
ax1.plot(culture_exp.system_time, culture_exp.community_dyn['BT'].values, '-', label='BT Model', color = 'b')
ax1.legend()

fig2, ax2 = plt.subplots()
ax2.plot(time_exp, glucose_exp, 'o', label = 'glucose (HPLC)', color='orange')
ax2.plot(culture_exp.system_time, culture_exp.environment_dyn['Glucose'], label = 'glucoseModel', color = 'orange')

ax2.plot(time_exp, pyruvate_exp, 'o', label = 'pyruvate (HPLC)', color='lightblue')
ax2.plot(culture_exp.system_time, culture_exp.environment_dyn['Pyruvate'], label = 'pyruvateModel', color = 'lightblue')

ax2.plot(time_exp, lactate_exp, 'o', label = 'lactate (HPLC)', color='red')
ax2.plot(culture_exp.system_time, culture_exp.environment_dyn['Lactate'], label = 'lactateModel', color = 'red')

ax2.plot(time_exp, succinate_exp, 'o', label = 'succinate (HPLC)', color='green')
ax2.plot(culture_exp.system_time, culture_exp.environment_dyn['Succinate'], label = 'succinateModel', color = 'green')

ax2.plot(time_exp, acetate_exp, 'o', label = 'acetate (HPLC)', color='b')
ax2.plot(culture_exp.system_time, culture_exp.environment_dyn['Acetate'], label = 'acetateModel', color = 'b')
ax2.legend()



stl.pyplot(fig1)
stl.pyplot(fig2)
#show()
   
