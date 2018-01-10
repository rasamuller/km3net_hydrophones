# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 17:43:52 2017

@author: mullerrs
"""

from matplotlib import pyplot as plt
import numpy as np
import sys

###############################################################################

def plot_ERES_N(col, X_MAX, ERES, CL, p_005, p_3s, p_5s, CL_err, p_005_err, p_3s_err, p_5s_err):
    X_MIN = 0
    #X_MAX = max(max(p_005[:], p_3s[:], p_5s[:]))
    Y_MIN = min(ERES)
    Y_MAX = max(ERES)
    
    # omweg om te bereiken wat ik wil met kleuren in de plot
    ERES = np.insert(ERES, 0, Y_MIN)
    ERES = np.append(ERES, Y_MAX)

    CL = np.insert(CL, 0, X_MIN)
    CL = np.append(CL, X_MAX)
    p_005 = np.insert(p_005, 0, X_MIN)
    p_005 = np.append(p_005, X_MAX)
    p_3s = np.insert(p_3s, 0, X_MIN)
    p_3s = np.append(p_3s, X_MAX)
    p_5s = np.insert(p_5s, 0, X_MIN)
    p_5s = np.append(p_5s, X_MAX)

    CL_err = np.insert(CL_err, 0, X_MIN)
    CL_err= np.append(CL_err, X_MAX)
    p_005_err = np.insert(p_005_err, 0, 0)
    p_005_err = np.append(p_005_err, 0)
    p_3s_err = np.insert(p_3s_err, 0, 0)
    p_3s_err = np.append(p_3s_err, 0)
    p_5s_err = np.insert(p_5s_err, 0, 0)
    p_5s_err = np.append(p_5s_err, 0)
    
    fig, ax = plt.subplots(figsize=(14,6))
    ax.errorbar(CL, ERES, 0, CL_err, '--', color='black', )
    ax.errorbar(p_005, ERES, 0, p_005_err, '.-', color='black')
    ax.errorbar(p_3s, ERES, 0, p_3s_err, '.-', color='black')
    ax.errorbar(p_5s, ERES, 0, p_5s_err, '.-', color='black')
    
    ax.fill_between(p_005, 0, ERES, facecolor=col, alpha=0.1)
    ax.fill_between(p_3s, 0, ERES, facecolor=col, alpha=0.3)
    ax.fill_between(p_5s, 0, ERES, facecolor=col, alpha=0.8)

    plt.text(p_005[2]+(p_005[3]-p_005[2])/4.*3., 70, r'$CL$', fontsize=18, weight='normal', bbox=dict(facecolor='white', alpha = 0.5, edgecolor='black', boxstyle='round,pad=1'))
    plt.text(p_005[2]+(p_005[3]-p_005[2])/4.*3. + 100, 50, r'$p = 0.05$', fontsize=20, weight='normal', bbox=dict(facecolor='white', alpha = 0.5, edgecolor='black', boxstyle='round,pad=0.3'))
    plt.text(p_3s[2]+(p_3s[3]-p_3s[2])/4.*2. + 100, 30, r'$3\sigma$', fontsize=20, weight='normal', bbox=dict(facecolor='white', alpha = 0.5, edgecolor='black', boxstyle='round,pad=0.3'))
    plt.text(p_5s[2]+(p_5s[3]-p_5s[2])/4.*1 + 100, 10, r'$5\sigma$', fontsize=20, weight='normal', bbox=dict(facecolor='white', alpha = 0.5, edgecolor='black', boxstyle='round,pad=0.3'))

    ax.set_xlim(0, X_MAX)
    ax.set_xlabel('Number of detections '+r'between $5 \cdot 10^{20} eV\ <\ E\ <\ 5 \cdot 10^{22}$ eV', fontsize=20)
    ax.set_ylabel('Energy resolution (%)', fontsize=20)
    

    plt.tick_params(labelsize=20)
    if col == 'green':
        plt.title('Number of neutrino detections for significant discovery/rejection \nwith top down source', fontsize=28, y=1.03)
    elif col == 'red':
        plt.title('Number of neutrino detections for significant discovery/rejection \nwith bottom up source', fontsize=28, y=1.03)
    plt.tight_layout()
    plt.show()


def main():
    ERES = [0, 30, 80]

    TD_CL = [1584, 1684, 1978]
    TD_p_005 = [1598, 1709, 2047]
    TD_p_3s = [3890, 4140, 4940]
    TD_p_5s = [10536, 11191, 13330]

    TD_CL_err = [22, 24, 29]
    TD_p_005_err = [23, 25, 31]
    TD_p_3s_err = [53, 57, 73]
    TD_p_5s_err = [142, 153, 196]

    BU_CL = [5832, 6881, 4928]
    BU_p_005 = [5890, 6897, 5143]
    BU_p_3s = [13789, 16083, 12027]
    BU_p_5s = [36692, 42716, 31985]

    BU_CL_err = [150, 198, 110]
    BU_p_005_err = [153, 198, 120]
    BU_p_3s_err = [357, 460, 278]
    BU_p_5s_err = [947, 1221, 737]

    X_MAX = max(max(TD_CL[:], TD_p_005[:], TD_p_3s[:], TD_p_5s[:], BU_CL[:], BU_p_005[:], BU_p_3s[:], BU_p_5s[:]))

    plot_ERES_N('green', X_MAX, ERES, TD_CL, TD_p_005, TD_p_3s, TD_p_5s, TD_CL_err, TD_p_005_err, TD_p_3s_err, TD_p_5s_err)
    plot_ERES_N('red', X_MAX, ERES, BU_CL, BU_p_005, BU_p_3s, BU_p_5s, BU_CL_err, BU_p_005_err, BU_p_3s_err, BU_p_5s_err)

if __name__ == "__main__":
    sys.exit(main())
