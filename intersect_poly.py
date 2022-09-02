"""Calculate intersections for corresponding swaths in two point clouds
Indice pairs were assessed visually.
"""
from funcs import *
from osgeo import ogr

pc_1 = r"E:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\parametrix\10_swaths_gps.laz"
pc_2 = r"E:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\parametrix\missisquoi_merged12.las"
outfile = r"E:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\parametrix\swath_comparison.txt"

pc_1 = r"F:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\teton_conservation\13_swaths_gps.laz"
pc_2 = r"F:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\teton_conservation\teton_merged12.laz"
outfile = r"F:\ws_21_22\bachelor_arbeit\data\otherOpenTopography\teton_conservation\swath_comparison_teton.txt"

pc_data1 = read_point_cloud(pc_1, subsample=1000, sort=True)
pc_data2 = read_point_cloud(pc_2, subsample=1000, sort=True)

alphas1, sizes1, widths1 = get_alphashapes(pc_data1, 13, 200)
alphas2, sizes2, widths2 = get_alphashapes(pc_data2, 13, 200)


indice_pairs = []
for shape in alphas1:
    alphavs1 = [xy.split(' ') for xy in shape.split('((')[1].split('))')[0].split(', ')]
    maxx1 = max([float(xy[0]) for xy in alphavs1])
    maxy1 = max([float(xy[1]) for xy in alphavs1])
    minx1 = min([float(xy[0]) for xy in alphavs1])
    print('x1: ' + str(maxx1))
    all_dists = []
    for shape2 in alphas2:
        alphavs2 = [xy.split(' ') for xy in shape2.split('((')[1].split('))')[0].split(', ')]
        maxx2 = max([float(xy[0]) for xy in alphavs2])
        minx2 = min([float(xy[0]) for xy in alphavs2])
        maxy2 = max([float(xy[1]) for xy in alphavs2])
        print('x2: ' + str(maxx2))
        dist = abs(maxx1-maxx2)+abs(maxy1-maxy2) + abs(minx1-minx2)
        print('distance: ' + str(dist))
        all_dists.append(dist)
    print(all_dists)
    index_min = all_dists.index(min(all_dists))
    indice_pairs.append(index_min)
print(indice_pairs)
indice_pairs = [9, 6, 3, 8, 11, 7, 1, 10, 2, 4, 12, 8, 0]
results = []

sizes1 = [6022450, 5581647, 5529826, 5486395, 5405276, 5145278, 3754686, 3470664, 1702487, 1425710]
sizes2 = [5321426, 5319194, 5318131, 5294309, 4949601, 4818860, 3506540, 3228378, 1682387, 1391386]

sizes1 = [4924585, 4848266, 4537861, 3867610, 3588820, 3465655, 3376999, 3263639, 2699668, 2567423, 2467096, 2257289, 2202213]
sizes2 = [4524965, 4460493, 4432226, 4431515, 4430786, 4418996, 4359112, 4350297, 4347566, 3874709, 3653184, 2628667, 2127450]

for i, index in enumerate(indice_pairs):
    intersect = intersect_poly(alphas1[i], alphas2[index])
    visualize_wkt([alphas1[i], alphas2[index]], i)

    size1 = ogr.CreateGeometryFromWkt(alphas1[i]).GetArea()
    size2 = ogr.CreateGeometryFromWkt(alphas2[index]).GetArea()
    size_intersect = ogr.CreateGeometryFromWkt(intersect).GetArea()
    print('\nSwath {}'.format(i))
    print('area1: {}\narea2: {}\nintersect: {}'.format(size1,
                                                       size2,
                                                       size_intersect))
    print('% of 1 covered by 2: {}'.format(size_intersect / size1 * 100))
    print('% of 2 covered by 1: {}'.format(size_intersect / size2 * 100))
    results.append([sizes1[i], size1, widths1[i], sizes2[index], size2, widths2[index], size_intersect, size_intersect / size1 * 100, size_intersect / size2 * 100])


with open(outfile, 'w') as f:
    f.write('Number of Points Original PC;Area Original PC;SW Original;Number of Points Sim;Area Sim;SW Sim;Area Intersection;% of original swath covered by simulated swath;'
            '% of simulated swath covered by original swath\n')
    for swath in results:
        f.write('{};{};{};{};{};{};{};{};{};\n'.format(swath[0], swath[1],swath[2], swath[3], swath[4], swath[5], swath[6], swath[7], swath[8]))
