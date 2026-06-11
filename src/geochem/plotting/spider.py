#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:56:39 2023

@author: Shavin Kaluthantri
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi
from scipy.stats import sem, norm
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.interpolate import LinearNDInterpolator
from src.utils.molecular import MolecularWeightCalculator
from src.geochem.processing.censored import gausscensor
from matplotlib.figure import Figure
import re

def plot_spider_norm(data, ref_data, norm_ref_data, layer,el_list=None, style='Quanta', quantiles=[0.25, 0.5, 0.75], ax=None, c='b', label=None):
    
    el_list_lower = [re.sub(r'\d', '', el).lower() for el in el_list]
    # Filter ref_dataerence data
    ref_data_chem = ref_data.loc[(ref_data['reference'] == norm_ref_data) | (ref_data['model'] == norm_ref_data) & (ref_data['layer'] == layer) & (ref_data['sigma'] == 0)]
    plot_inf = {}
    if el_list is None:
        el_list = ['Cs', 'Rb', 'Ba', 'Th', 
        'U', 'Nb', 'Ta', 'K', 
        'La', 'Ce', 'Pb', 'Mo',
        'Pr', 'Sr', 'Ga', 'Zr', 'Hf',
        'Nd', 'Sm', 'Eu', 'Li', 
        'Ti', 'Gd', 'Tb', 'Dy', 
        'Y', 'Ho', 'Er',
        'Tm', 'Yb', 'Lu', 'Zn', 
        'Mn', 'V', 'Sc', 'Co', 
        'Cu', 'Ni', 'Cr']
    
    
    # If no axes are provided, create new axes
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    else:
        fig = ax.figure


    # Number of variables we're plotting
    num_vars = len(el_list)
    
    if not any('ppm' in col for col in data.columns):
        # data.columns = [col[:2].lower() for col in data.columns]
        if not ref_data.empty:
            ref_data.columns = [col.lower() for col in ref_data.columns]
        ref_data_chem.columns = [col.replace('_ppm', '') for col in ref_data_chem.columns]
    else:
        ref_data = ox2ppm(el_list, ref_data)
        
        
    # ref_data = ox2ppm(el_list, ref_data)
    # el_list = [elemts.lower()]
    
    ymodel = {'mu': np.nan, 'sigma': np.nan}
    
    results = []
    for el in el_list:
        
        y= data[el].dropna().values
        if len(y)==0:
            print (el)
        if all(y>0): 
            if style == 'MeanSD' or style == 'MeanSE':
                ymodel['mu'] = np.mean(np.log10(y))
                ymodel['sigma'] = np.std(np.log10(y))
            elif style=='Quanta':
                y_q = np.quantile(np.log10(y), q=quantiles)
        else: #censorred data
            ymodel, y_q = gausscensor(y,scale = 'log',q=quantiles); 
            y_q = y_q[:,1]
        
        
        if style == 'MeanSD' or style == 'MeanSE':
            results.append([ymodel['mu'], ymodel['sigma']])
        elif style=='Quanta':
            results.append(np.append(y_q,len(y)))
    
    result_df = pd.DataFrame(results, columns=[ 'mu','sigma'] if style != 'Quanta' else ['Q' + str(int(q*100)) for q in quantiles] + ['N'], index = el_list_lower)
    
    ref_series = ref_data_chem[el_list_lower].squeeze()
    
    
    if style in ['MeanSD', 'MeanSE']:
            
        result_df['mu_norm'] = np.log10(10**result_df['mu'].div(ref_series, axis=0))
        result_df['sigma_norm'] = result_df['sigma'].div(np.log10(ref_series))
        
        result_df['se'] = result_df['sigma']/len(result_df['sigma'])
        result_df['se_norm'] = result_df['sigma_norm']
        
        plot_inf['dimention_names'] = ['Element', 'Concentration Data'] 
        plot_inf ['variable_units'] = ['log10 (ppm)']
        
    elif style == 'Quanta':
        result_df[result_df.columns[:-1]] = np.log10((10**result_df[result_df.columns[:-1]]).div(ref_series, axis=0))
        
        plot_inf ['variable_names'] = result_df.columns[:-1]
    
    result_df['elid'] = range(0, len(result_df))
    
    # fig = Figure()
    # ax = fig.add_subplot(111)
    ax, yl = plot_data(ax,result_df,c, style,el_list,Q=quantiles, label = label) 
    
    # add element list to results to add more information to export
    result_df['element'] = el_list
    
    return ax, yl, result_df
    # df = pd.DataFrame(results, columns=['N', 'mu', 'mu_norm', 'sigma', 'sigma_norm'] if style != 'Quanta' else ['Q' + str(int(q*100)) for q in Q] + ['N'])
    
def logax(ax, lim, axis='y', label='', tick_label_rotation=0):
    """
    Produces log-axes limits and labels.

    Parameters:
    ax (matplotlib.axes.Axes): The axes to modify.
    lim (list): The log10 values of the axes limits.
    axis (str): 'x' or 'y' to add ticks to x- or y-axis, default is 'y'.
    label (str): Label for the axis.
    tick_label_rotation (float): Angle of text rotation, default is 0.
    """
    # Create tick marks and labels
    mt = np.log10(np.arange(1, 10))
    ticks = []
    tick_labels = []
    for i in range(int(lim[0]), int(lim[1]) + 1):
        ticks.extend([i + m for m in mt])
        tick_labels.extend([f'{10 ** i}'] + [''] * (len(mt) - 1))

    # Apply settings based on the axis
    if axis.lower() == 'x':
        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels, rotation=tick_label_rotation)
        ax.set_xlim([10**lim[0], 10**lim[1]])
        if label:
            ax.set_xlabel(label)
    elif axis.lower() == 'y':
        ax.set_yticks(ticks)
        ax.set_yticklabels(tick_labels, rotation=tick_label_rotation)
        ax.set_ylim([10**lim[0], 10**lim[1]])
        if label:
            ax.set_ylabel(label)
    else:
        print('Incorrect axis argument. Please use "x" or "y".')
    
    
    
def fillinterval(ax, x, y, color, alpha):


    # Flattening and combining x and y for filling the interval
    # x[:,1] =  x[:,1][::-1]
    # y[:,1] =  y[:,1][::-1]

    # Filling the interval in the plot
    ax.fill_between(x[:,0], y[:,0],y[:,1], color=color, alpha=alpha, edgecolor=None)

def plot_data(ax, t, C, style,el_list, Q=None, label= None):
    x = np.hstack([t.loc[:,['elid']].values, t.loc[:,['elid']].values])
    if style == 'MeanSD':
        ind = ~np.isnan(t.loc[:,'mu_norm'])
        x = t.loc[ind, 'elid']
        y1 = np.vstack([t.loc[ind, 'mu_norm'] - t.loc[ind, 'sigma_norm'], 
                        t.loc[ind, 'mu_norm'] + t.loc[ind, 'sigma_norm']])
        y2 = np.vstack([t.loc[ind, 'mu_norm'] - 2 * t.loc[ind, 'sigma_norm'], 
                        t.loc[ind, 'mu_norm'] + 2 * t.loc[ind, 'sigma_norm']])

        fillinterval(ax, x, y2, C, 0.15)
        fillinterval(ax, x, y1, C, 0.3)
        ax.plot(x, t.loc[ind, 'mu_norm'], color=C, linewidth=0.75)

    
    elif style == 'MeanSE':
        ind = ~np.isnan(t.loc[:,'mu_norm'])
        y1 = np.vstack([t.loc[ind,'mu_norm'] - t.loc[ind,'se_norm'], t.loc[ind,'mu_norm'] + t.loc[ind,'se_norm']]).T
        y2 = np.vstack([t.loc[ind,'mu_norm'] - 2*t.loc[ind,'se_norm'], t.loc[ind,'mu_norm'] + 2*t.loc[ind,'se_norm']]).T 

        fillinterval(ax, x, y2, C, 0.15)
        fillinterval(ax, x, y1, C, 0.3)
        ax.plot(x, t.loc[ind,'mu_norm'], color=C, linewidth=0.75)
    
    elif style == 'Quanta':
        Q = t.columns[:-2]
        if len(Q) == 1:
            ind = ~t.loc[:, Q[0]].isna()
            ax.plot(t.loc[ind, 'elid'], t.loc[ind, Q[0]], color=C, linewidth=0.75, label=label)
            yl = [np.floor(np.nanmin(t.iloc[:, 0])), np.ceil(np.nanmax(t.iloc[:, 0]))]
        elif len(Q) == 2:
            ind = ~t.loc[:, Q[0]].isna()
            y = t.loc[ind, [Q[0], Q[1]]].values
            fillinterval(ax, x[ind], y, C, 0.3)
            yl = [np.floor(np.nanmin(t.iloc[:, 0])), np.ceil(np.nanmax(t.iloc[:, 1]))]
        elif len(Q) == 3:
            ind = ~t.loc[:, Q[1]].isna()
            y = t.loc[ind, [Q[0], Q[2]]].values
            fillinterval(ax, x[ind], y, C, 0.3)
            ax.plot(x[ind,0], t.loc[ind, Q[1]], color=C, linewidth=0.75, label=label)
            yl = [np.floor(np.nanmin(t.iloc[:, 0])), np.ceil(np.nanmax(t.iloc[:, 2]))]
        elif len(Q) == 5:
            ind = ~t.loc[:, Q[2]].isna()
            y = t.loc[ind, [Q[0], Q[4]]].values
            ax.plot(x[ind], y, color=C, linewidth=0.25)
            y = t.loc[ind, [Q[1], Q[3]]].values
            fillinterval(ax, x[ind], y, C, 0.3)
            ax.plot(x[ind], t.loc[ind, Q[2]], color=C, linewidth=0.75, label=label)
            yl = [np.floor(np.nanmin(t.iloc[:, 0])), np.ceil(np.nanmax(t.iloc[:, 4]))]
        else:
            for i in Q:
                ax.plot(t['elid'], t.loc[:, i], color=C)
    
 # Set y-limits for the axes
    # yl = [np.floor(t.min().min()), np.ceil(t.max().max())]
    # ax.set_ylim(yl)   
    ax.set_xticks(x[:,1])
    ax.set_xticklabels(el_list, rotation=45)
    logax(ax, yl, 'y')
    ax.set_ylim(yl) 
    return ax, yl


# # Repeat first value to close the circle in the plot
# norm_data.loc[len(norm_data)] = norm_data.iloc[0]

# # Draw one axe per variable and add labels
# labels = el_list + [el_list[0]]
# ax.set_xticks(angles)
# ax.set_xticklabels(labels)

# if style == 'MeanSD':
#     means = norm_data.mean()
#     std_devs = norm_data.std()
#     means = np.concatenate((means.values, [means.values[0]]))
#     std_devs = np.concatenate((std_devs.values, [std_devs.values[0]]))
#     ax.plot(angles, means, 'b-', linewidth=1.5, label='Mean')
#     ax.fill_between(angles, means - std_devs, means + std_devs, color='blue', alpha=0.3, label='SD')

# elif style == 'MeanSE':
#     means = norm_data.mean()
#     errors = sem(norm_data, axis=0)
#     means = np.concatenate((means.values, [means.values[0]]))
#     errors = np.concatenate((errors, [errors[0]]))
#     ax.plot(angles, means, 'r-', linewidth=1.5, label='Mean')
#     ax.fill_between(angles, means - errors, means + errors, color='red', alpha=0.3, label='SE')

# elif style == 'Quanta':
#     for q in quantiles:
#         quantile_values = norm_data.quantile(q)
#         quantile_values = np.concatenate((quantile_values.values, [quantile_values.values[0]]))
#         q = int(q*100)
#         ax.plot(angles, quantile_values, linewidth=1, linestyle='solid', label=f'{q}th Percentile')

# else:
#     raise ValueError('Unknown style specified')

# ax.legend(loc='upper right')


def ox2ppm(element_list, data):
    mol_cal = MolecularWeightCalculator()
    for element in element_list:
        oxide = element + "O"
        if oxide in data.columns:
            molar_mass_element = mol_cal.molecular_weight(element)
            molar_mass_oxide = mol_cal.molecular_weight(oxide)
            conversion_factor = {
                'K': 2, 'Na': 2, 'Al': 2, 'Fe': 1, 'Ni': 1, 'Mn': 1, 'Cr': 2/3, 'P': 2
            }.get(element, 1)
            data[element.lower() + '_ppm'] = conversion_factor * molar_mass_element / molar_mass_oxide * data[oxide] * 10000
    return data

# ref_data = pd.read_excel('earthref.xlsx')

# # data = pd.read_csv('/Users/a1904121/LaserMapExplorer/laser_mapping/Alex_garnet_maps/processed data/RM01.csv')
# data = pd.read_csv('/Users/shavinkalu/Library/CloudStorage/GoogleDrive-a1904121@adelaide.edu.au/.shortcut-targets-by-id/1r_MeSExALnv9lHE58GoG7pbtC8TOwSk4/laser_mapping/Alex_garnet_maps/processed data/RM01.csv')
# el_list = data.columns[5:10]


# i = 1


# # # plot_spider_norm(data, ref_data, norm_ref_data, layer,el_list=None, style='Quanta', quantiles=[0.05, 0.25, 0.5, 0.75, 0.95], ref_data_field='sio2', ref_data_val='median', ax=None)
# fig, ax = plot_spider_norm(data = data[el_list.values],ref_data = ref_data,norm_ref_data =  ref_data['model'][i], layer = ref_data['layer'][i],el_list= el_list ,style = 'Quanta')

# plt.show()
