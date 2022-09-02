from funcs import *
import statistics
import os
import datetime
import time

filenames = glob.glob('*.laz')
filenames = [r"G:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\teton_conservation\points.laz"]
#filenames = [r"G:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\parametrix\points_true.laz"]
results_file = os.path.join(os.getcwd(), 'results_{}.txt'.format(time.time()))
tab_results = os.path.join(os.getcwd(), 'tab_results_{}.txt'.format(time.time()))
functions_file = r'E:\ws_21_22\bachelor_arbeit\coding\funcs.py'
las_sort = r"E:\LAStools\bin\lassort.exe"
alt_results = []
v_results = []
alpha_rad = 200

'''with open(results_file, 'a') as f:
    f.write('Test of altitude and flight speed detection functions.\n')
    f.write('Time : {}\n'.format(datetime.datetime.now()))
    f.write('Most recent edit of funcs.py - {}\n'.format(datetime.datetime.fromtimestamp(os.path.getmtime(functions_file)).strftime('%c')))
    f.write("Alpharadius: {}\n".format(alpha_rad))
    f.write("Alphaval: {}\n\n".format(1/alpha_rad))'''

visualize = False
visualze_swaths = True
#sort_gpstime(las_sort, r"E:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\terrapoint\points5_copy.laz")
devs = [[], []]
for j, point_cloud in enumerate(filenames):
    #metadata_file = 'metadata_{}.txt'.format(point_cloud.split('merged_')[1].split('.laz')[0])
    metadata_file = 'metadata_16.txt'
    true_vals = read_true_values(metadata_file)

    data = read_point_cloud(point_cloud, subsample=1000)
    results = estimate_alt_v_exp(data, num_swaths=18, scan_angle=true_vals['scan_angle'], alpha_rad=alpha_rad,
                                 visualize=visualize, visualize_swaths=visualze_swaths, ft=False)
    #data = read_point_cloud(point_cloud, stop=1000)
    #prf = estimate_prf(data)
    #write_results(results_file, point_cloud, devs, true_vals, results, prf)
    #write_results_tab_no_sim(tab_results, point_cloud, true_vals, results, prf)


#write_overall_results(results_file, devs)