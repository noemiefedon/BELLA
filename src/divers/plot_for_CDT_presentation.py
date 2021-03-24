# -*- coding: utf-8 -*-
"""
scipt to prepare plots for CDT presentation
April 2019
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
sys.path.append(r'C:\BELLA')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from src.CLA.lampam_functions import calc_lampam
from src.CLA.ABD import D_from_lampam
from src.BELLA.materials import Material
from src.buckling.buckling import buckling_factor

#==============================================================================
# Material properties
#==============================================================================
# Elastic modulus in the fibre direction in Pa
E11 = 20.5/1.45038e-10  # 141 GPa
# Elastic modulus in the transverse direction in Pa
E22 = 1.31/1.45038e-10  # 9.03 GPa
# Poisson's ratio relating transverse deformation and axial loading (-)
nu12 = 0.32
# In-plane shear modulus in Pa
G12 = 0.62/1.45038e-10 # 4.27 GPa
# Density in g/m2
density_area = 300.5
# Ply thickness in m
ply_t = (25.40/1000)*0.0075 # 0.191 mmm

mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
               density_area=density_area, ply_t=ply_t)
#==============================================================================
# Stacking sequences
#==============================================================================
n_values = 100

theta = np.linspace(-90, 90, n_values, endpoint=True)

lampam = np.zeros((n_values, 12), float)
D11 = np.zeros((n_values,), float)
D22 = np.zeros((n_values,), float)
D12 = np.zeros((n_values,), float)
D66 = np.zeros((n_values,), float)
buck = np.zeros((n_values,), float)

for index in range(n_values):
    ss = np.zeros((1,), float)
    ss[0] = theta[index]
    lampam[index] = calc_lampam(ss)
    D = D_from_lampam(lampam[index], mat)
    D11[index] = D[0, 0]
    D22[index] = D[1, 1]
    D12[index] = D[0, 1]
    D66[index] = D[2, 2]
    buck[index] = buckling_factor(lampam[index],
                   mat,
                   n_plies = 1,
                   N_x= 10,
                   N_y= 10,
                   length_x=0.1,
                   length_y=0.1)

fig, ax = plt.subplots(1,1, sharex=True, figsize=(8,5.8))
# sharex=True to align x-axes of the graphs
# figsize size of the combines sub-plots
my_labelsize = 20
my_titlesize = 40
my_axissize = 26
my_font = 'Arial'
# color of line of the graphs, boxes
ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')
ax.spines['left'].set_color('black')
# set the ticks
ax.tick_params(
    direction='out',
    length = 0, width = 4,
    colors = 'black',
    labelsize = my_labelsize)
ax.yaxis.set_ticks([0,00, 0.04, 0.08])
# plot
ax.plot(theta, D11, 'b-')
# set axes' labels
ax.set_ylabel(
  r'$\mathrm{Flexural\ stiffnesses\ [Pa.m^{-3}]}$',
  fontweight="normal",
  fontsize = my_axissize,
  fontname=my_font)
ax.plot(theta, D22, 'g-')
ax.yaxis.set_ticks([0,00, 0.04, 0.08])
ax.plot(theta, D12 + 2*D66, 'r-')
ax.set_xlabel(
    r'$\mathrm{Fibre\ orientation\ \theta\  [^{\circ}]}$',
    fontsize =my_axissize,
    fontname=my_font)
ax.yaxis.set_ticks([0,00, 0.02, 0.04, 0.06, 0.08])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0f'))
ax.xaxis.set_ticks([-90, -75,-60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90])
# format labels
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.2f'))
# white background
ax.set_facecolor((1, 1, 1))
# labels fonts
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname(my_font) for label in labels]
# axes' ranges
ax.set_xlim([-90, 90]) #ax.set_ylim([-1, 1])
# grid
ax.grid(False)
# aspect ratio ??
ax.set_aspect('auto')
# padding when several subplots
plt.tight_layout(pad=1, w_pad=1, h_pad=3)
# legends
ax.legend((r'$\mathrm{D_{11}}$',
           r'$\mathrm{D_{22}}$',
           r'$\mathrm{D_{12} + 2*D_{66}}$'),
          handletextpad=0,
          loc='right',
          fontsize =my_labelsize)
# padding for the labels
ax.tick_params(pad = 10)
# ???
#plt.tight_layout(pad=20, w_pad=20, h_pad=20)
# save image
plt.savefig('plot_Ds.svg')

fig, ax = plt.subplots(1,1, sharex=True, figsize=(9,3))
my_labelsize = 20
my_titlesize = 40
my_axissize = 26
my_font = 'Arial'
ax.spines['bottom'].set_color('black')
ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')
ax.spines['left'].set_color('black')
ax.tick_params(direction='out',
               length = 0, width = 4,
               colors = 'black', labelsize = my_labelsize,
               pad = 10)
ax.plot(theta, buck)
ax.set_xlabel(
            r'$\mathrm{\theta}$',
            fontsize =my_axissize,
            fontname=my_font) # X label
ax.set_ylabel(
  r'$\mathrm{Buckling\ \ factors}$', fontweight="normal",
  fontsize = my_axissize,
  fontname=my_font) # Y label
ax.yaxis.set_ticks([4, 5, 6, 7, 8, 9])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0f'))
ax.xaxis.set_ticks([-90, -75,-60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90])
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0f'))
ax.set_facecolor((1, 1, 1))
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname(my_font) for label in labels]
ax.set_xlim([-90, 90])
ax.grid(False)
plt.savefig('plot_buck.svg')
