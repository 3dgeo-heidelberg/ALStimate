"""Check las point clouds for unique values in columns.
Reads each column present in data. If all values in column are the same sets value in variation flag column to 0. Else to 1."""
import laspy
import numpy as np
from tabulate import tabulate

def analyse_pc(file, sample_size):
    print('Analysing point cloud: ' + file)
    with laspy.open(file) as f:
        point_format = f.header.point_format
        point_count = f.header.point_count
        print("Number of points:   {}".format(point_count))
        #print(list(point_format.standard_dimension_names))

        variation_flags = []

        print('Reading file... Getting {} points. (Every {}th point)'.format(sample_size, int(point_count/sample_size)))
        las = laspy.read(file)[::int(point_count/sample_size)]
        '''for val in las.scan_angle_rank:
            print(val)'''
        for attr in list(point_format.standard_dimension_names):
            #print(attr)
            current_vals = np.array([las[attr]])
            #print(current_vals)
            different_vals = 0
            for i in range(len(current_vals[0])-1):
                if current_vals[0, i] != current_vals[0, i+1] and not different_vals:
                    different_vals = 1

            variation_flags.append([attr, different_vals])

        print(tabulate(variation_flags, headers=['Attribute', 'Variation Flag'], tablefmt="grid"))

#analyse_pc(r'D:\ws_21_22\bachelor_arbeit\data\espana_total\data\PNOA_2016_GAL_E_648-4776_ORT-CLA-CIR.laz', 2000)

datasets = [r'D:\ws_21_22\bachelor_arbeit\data\norway_small\data\merged.las',
            r'D:\ws_21_22\bachelor_arbeit\data\antarctica_small\points.laz',
            r'D:\ws_21_22\bachelor_arbeit\data\sanandreas_sanjacinto_first\points.laz',
            r'D:\ws_21_22\bachelor_arbeit\data\pais_vasco\LT20_08_au_RedNAP08.laz',
            r'D:\ws_21_22\bachelor_arbeit\data\syssifoss\ALS_BR01_2019-07-05_300m.laz',
            r'D:\ws_21_22\bachelor_arbeit\data\santorini\data\merged.las',
            r'D:\ws_21_22\bachelor_arbeit\data\espana_total\galicia_2016\PNOA_2016_GAL_E_648-4776_ORT-CLA-CIR.laz']

for dataset in datasets[4:]:
    analyse_pc(dataset, 2000)


#with laspy.open(datasets[0]) as f:
 #   point_format = f.header.point_format
  #  point_count = f.header.point_count
   # print("Number of points:   {}".format(point_count))
    #for points in f.chunk_iterator(100):
     #   print(type(points))
      #  for i in points:
       #     print(i.X)
        #break
