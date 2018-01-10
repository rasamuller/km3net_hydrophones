# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 14:46:12 2017

@author: mullerrs
"""

import matplotlib.pyplot as plt


plt.plot([10, 20, 30, 50, 60],[0.12, 0.62, 0.99, 2.20, 4.41], 'b', label = '1Vpp')
plt.plot([10, 20, 30, 50, 60], [0.37, 1.87, 2.96, 6.63, 13.23], 'g', label = '3Vpp')
plt.plot([10, 20, 30, 50, 60], [0.62, 3.11, 4.94, 11.05, 22.05], 'r', label = '5Vpp')
plt.plot([10, 20, 30, 50, 60], [0.87, 4.36, 6.91, 15.47, 30.86], 'y', label = '7Vpp')
plt.plot([10, 20, 30, 50, 60], [1.24, 6.23, 9.87, 22.01, 44.09], 'k', label = '10Vpp')
plt.xlabel('frequency')
plt.legend()
plt.show()
plt.clf()