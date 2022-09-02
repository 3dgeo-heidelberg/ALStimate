"""Plots and compares hists of distribution of point along x y and z axis for two point clouds"""
import matplotlib.pyplot as plt
from funcs import *

#original = r"E:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\parametrix\10_swaths.laz"
#simulated = r"E:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\parametrix\missisquoi_merged12.las"
original1 = r"E:\ws_21_22\bachelor_arbeit\helios-1.1.0\output\Survey Playback\missisquoi\2022-05-12_16-04-00\points\strip004_point.laz"
simulated1 = r"E:\ws_21_22\bachelor_arbeit\helios-1.1.0\output\Survey Playback\missisquoi\2022-05-12_16-04-00\points\strip003_point.laz"

original = r"F:\ss_22\bachelor_arbeit\data\otherOpenTopography\teton_conservation\13_swaths.laz"
simulated =r"F:\ss_22\bachelor_arbeit\data\otherOpenTopography\teton_conservation\teton_merged12.laz"

data1 = read_point_cloud_xyz(original, first10=True)
data2 = read_point_cloud_xyz(simulated, small_caps=True, scale_offset=False, first10=True)

#for teton
data1 = data1[:, data1[2, :] < 3000]

#for missisquoi
#data1 = data1[:, data1[2, :] < 120]

fig, axs = plt.subplots(2, 3, figsize=(15, 9), sharey=True)

for i in range(2):
    for k in range(3):
        axs[i][k].grid()
        axs[i][k].set_axisbelow(True)

axs[0][0].hist(data1[0,:], color='darkorange')
#axs[0][0].set_title('X Distribution')
axs[0][0].set_ylabel('Number of Points', fontsize=18)#, fontweight="bold")#\n(Input Point Cloud)')
#axs[0][0].set_xlabel('X')

axs[0][1].hist(data1[1,:], color='darkorange')
axs[0][1].set_title('Input Point Cloud', fontsize=18, fontweight="bold")
#axs[0][1].set_xlabel('Y')

axs[0][2].hist(data1[2,:], color='darkorange')
#axs[0][2].set_title('Z Distribution')
# For missisquoi:
#axs[0][2].set_xlim(10,130)
#axs[0][2].set_xlabel('Z')

axs[1][0].hist(data2[0,:])
#axs[1][0].set_title('X Distribution')
axs[1][0].set_ylabel('Number of Points', fontsize=18)#\n(Simulated Point Cloud)')
axs[1][0].set_xlabel('X', fontsize=18)

axs[1][1].hist(data2[1,:])
axs[1][1].set_title('Simulated Point Cloud', fontsize=18, fontweight="bold")
axs[1][1].set_xlabel('Y', fontsize=18)

axs[1][2].hist(data2[2,:])
#axs[1][2].set_title('Z Distribution')
axs[1][2].set_xlabel('Z', fontsize=18)
# For missisquoi:
#axs[1][2].set_xlim(10,130)

#fig.suptitle('Verification of Parameter Estimation Results for: Teton Conservation District', fontweight="bold", fontsize=18)
#fig.title('Verification of Parameter Estimation Results for: Teton Conservation District',fontsize=10)
#\nComparison of X, Y and Z Distributions between Input and Simulated Point Clouds')
for n, label in enumerate(axs[0][1].xaxis.get_ticklabels()):
    if n % 2 != 1:
        label.set_visible(False)
for n, label in enumerate(axs[1][1].xaxis.get_ticklabels()):
    if n % 2 != 1:
        label.set_visible(False)


#plt.setp(axs[0][1].get_xticklabels(), rotation=30, horizontalalignment='right')
#plt.setp(axs[1][1].get_xticklabels(), rotation=30, horizontalalignment='right')
#plt.locator_params(axis='x', nbins=5)
plt.tight_layout()
plt.show()
