#!/usr/local/bin/python3

import math
import numpy as np
import matplotlib.pyplot as plt

ch1_1 = np.loadtxt("ch1-1.txt", delimiter=" ", dtype=[('x', 'f8'), ('num', 'f8')])
ch1_2 = np.loadtxt("ch1-2.txt", delimiter=" ", dtype=[('x', 'f8'), ('num', 'f8')])
ch2_2 = np.loadtxt("ch2-2.txt", delimiter=" ", dtype=[('x', 'f8'), ('num', 'f8')])
ch2_1 = np.loadtxt("ch2-1.txt", delimiter=" ", dtype=[('x', 'f8'), ('num', 'f8')])


sums = []
for plot in [ch1_1, ch1_2, ch2_2, ch2_1]:
	totalNum = 0
	for row in plot:
		totalNum += row['num']
	sums.append(totalNum)
	
for i, plot in enumerate([ch1_1, ch1_2, ch2_2, ch2_1]):
	for row in plot:
		row['num'] /= sums[i]

f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)

for axis, plot, title in zip([ax1, ax2, ax3, ax4], [ch1_1, ch1_2, ch2_2, ch2_1], ["PMT1; 1 Start", "PMT1; 2 Start", "PMT2; 2 Start", "PMT2; 1 Start"]):
	axis.plot([row['x'] for row in plot if row['x'] < 10000], [row['num'] for row in plot if row['x'] < 10000], linewidth=0.5, color='black')
	axis.set_yscale("log")
	axis.set_xlim((0, 10000))
	axis.tick_params(labelsize=12)
	axis.text(7500, 1e-3, title, fontdict={'fontsize':12})

ax1.set_title("Photon Arrival Time", fontsize=12)
ax1.set_ylabel("Counts [arb.]", fontsize=12)
ax4.set_xlabel("MCS Clock [x800 ps]", fontsize=12, x=0.85)

f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

#plt.show()
plt.savefig("Resample_Photon_Arrival.pdf")