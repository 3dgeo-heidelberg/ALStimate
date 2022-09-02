# Author: Mark Searle, 3DGeo Research Group Heidelberg
# Contact: mark.searle@stud.uni-heidelberg.de
# -*- coding: utf-8 -*-
"""Script used to generate HELIOS++ surveys for flights over teton conservation dataset
Flight lines were entered manually from estimation
"""
import random
import os


helios_path = ''
lastools = r"G:\LAStools\bin\lasmerge.exe"
outdir = r"G:/ws_21_22/bachelor_arbeit/data/generated_data/test_data3"
height_variation = 60


def world2Pixel(gt, x, y):
    ulX = gt[0]
    ulY = gt[3]
    pixel_width = gt[1]
    pixel_height = -gt[5]
    rtnX = gt[2]
    rtnY = gt[4]
    row = int((x - ulX) / pixel_width)
    col = int((ulY - y) / pixel_height)
    #print("x: {} - {} / {} = {}".format(x, ulX, pixel_width, row))
    #print("y: {} - {} / {} = {}".format(ulY, y, pixel_height, col))
    return (row, col)


for i in range(1):
    # Use flight_planner tool from helios to generate the legs
    # Use XML parse to select a scanner
    # Use random.randint for scan angle, altitude, velocity
    survey_file_name = 'data/surveys/teton_01.xml'
    scene = 'data/scenes/demo/teton.xml#teton'

    # Generate random survey parameters.
    platform = "data/platforms.xml#sr22"
    speed =  72 # in m/s - roughly estimated from point cloud
    altitude = 1393
    scanner = "data/scanners_als.xml#leica_als50"
    pulse_freq = 76269
    scan_freq = 36
    scan_angle = 20
    survey_name = 'teton'
    mean_height = -27
    trajectory_miss = [[[450356.6194915057, 262316.7544264039], [448498.14187853027, 258823.4744493477]],
                [[452112.32948792336, 262310.97639677546], [450255.11765544454, 258825.26707968133]],
                [[450961.06614172744, 262332.0863335057], [449085.53403192054, 258809.949702514]],
                [[452580.11567910167, 262099.6833243992], [450832.04363462504, 258818.7267768149]],
                [[451514.5746349774, 262307.47585786635], [449677.2627814148, 258816.2782849349]],
                [[449756.3840071246, 262326.5123338983], [448001.7502938127, 259000.96388646305]],
                [[452827.79454316257, 261496.841461925], [451434.2983415371, 258831.33419390657]],
                [[449169.5180374859, 262319.27098827856], [447664.22632246604, 259404.84337431696]],
                [[448628.05769253883, 262328.86696802406], [447667.0554880783, 260537.03797890784]],
                [[452825.90613454476, 260391.52972196692], [452004.6937500057, 258828.88220930172]]]
    trajectory = [[[507462.079328874, 4822383.532786906], [505659.7129898629, 4818435.345173506]],
                [[509547.2982111833, 4822362.578434537], [507564.36735469115, 4818454.004510125]],
                [[510349.54561082204, 4821418.816993396], [508965.6958481325, 4818517.063056077]],
                [[506834.4889249103, 4822413.194812896], [504853.97623804235, 4818493.354266439]],
                [[508810.0893951234, 4822370.461200926], [506941.4127961871, 4818451.764528735]],
                [[506118.6354790516, 4822367.725245757], [504330.7129018605, 4818715.385905422]],
                [[508797.59405889444, 4822296.578351715], [507021.74550570187, 4818530.917455384]],
                [[506878.4721886738, 4822392.657034541], [504786.8949214745, 4818558.178635273]],
                [[505442.33658443124, 4822389.474067379], [504244.29853682045, 4819753.253509631]],
                [[509479.4925847404, 4822287.127059998], [507669.6597556056, 4818533.776227919]],
                [[507517.53888838354, 4822424.008536536], [505459.95591092145, 4818564.62235054]],
                [[510383.78529091715, 4820171.887574716], [509671.99346106965, 4818502.99712624]],
                [[504903.42253106367, 4822340.435825322], [504270.88253106386, 4821116.3658253215]]]

    #from osgeo import ogr
    from osgeo import ogr, gdal, osr

    input_raster = r"data/sceneparts/teton_3m.tif"
    raster = gdal.Open(input_raster)
    proj = osr.SpatialReference(wkt=raster.GetProjection())
    print(proj.GetAttrValue('AUTHORITY',1))
    geotransform = raster.GetGeoTransform()
    raster_b1 = raster.GetRasterBand(1)
    print(geotransform)


    # Get trajectory from vector file.
    all_points = []
    for stripid, feature in enumerate(trajectory):
        for i, pt in enumerate(feature):
            if i==0:
                all_points.append([pt[0], pt[1], altitude, "true", stripid])
            else:
                all_points.append([pt[0], pt[1], altitude, "false", stripid])

    #all_points_corrected = all_points
    false_points = 0
    # Update trajectory values with values over raster
    all_points_corrected = []
    # Update altitude.
    for i in range(len(all_points)):

        # Only if leg is active..:
        if all_points[i][3] == "true":
            print('Trajectory active')
            #print('Moving between points {} and {}'.format(all_points[i], all_points[i+1]))#
            print('Accessing point: X={}, Y={}'.format(all_points[i][0], all_points[i][1]))
            try:
                # Get raster height
                pixel, line = world2Pixel(geotransform, all_points[i][0], all_points[i][1])
                print(pixel, line)
                value = raster_b1.ReadAsArray(pixel, line, 1, 1)[0, 0]
                # Add raster height to flight alt at beginning of trajectory.
                all_points[i][2] += value
                # Add point entry to list with new corrected altitude values.
                all_points_corrected.append([all_points[i][0], all_points[i][1], all_points[i][2], all_points[i][3], all_points[i][4]])
            except TypeError:
                print('Point is outside of raster.. Assigning standard height value of {}.'.format(mean_height))
                all_points_corrected.append([all_points[i][0], all_points[i][1], mean_height, all_points[i][3], all_points[i][4]])
                false_points+=1
            # For every leg, check for altitude 100 times.
            frequency = 100
            current_loc = all_points[i]
            current_alt = all_points[i][2]
            # Distance to be moved between each altitude check.
            move_per_step = [(all_points[i + 1][0] - all_points[i][0]) / frequency,
                             (all_points[i + 1][1] - all_points[i][1]) / frequency]
            print("Move per step is {}".format(move_per_step))
            # Perform altitude checks, update altitude value if necessary.
            current_gnd_z = 0
            for j in range(frequency):
                try:
                    pixel, line = world2Pixel(geotransform, current_loc[0], current_loc[1])
                    current_gnd_z = raster_b1.ReadAsArray(pixel, line, 1, 1)[0, 0]
                    distance_to_optimal = current_gnd_z + altitude - current_alt
                    if abs(distance_to_optimal) > height_variation:
                        # If altitude is out of range:
                        print('new alt: {}'.format(current_gnd_z + altitude))
                        # Update value in survey leg.
                        all_points_corrected.append(
                            [current_loc[0], current_loc[1], current_gnd_z + altitude, "true", all_points[i][4]])

                        try:
                            if all_points_corrected[-2][2] == mean_height:
                                all_points_corrected[-2][2] = all_points_corrected[-1][2]
                                for k in range(3, 30):
                                    if all_points_corrected[-k][2] == mean_height:
                                        all_points_corrected[-k][2] = all_points_corrected[-1][2]
                                    else:
                                        break
                        except IndexError:
                            pass

                        # Update value for upcoming checks.
                        current_alt = current_gnd_z + altitude

                except TypeError:
                    if current_gnd_z==0:
                        print('Point is outside of raster.. Assigning previous height value of {}.'.format(mean_height))
                        all_points_corrected.append(
                            [all_points[i][0], all_points[i][1], mean_height, all_points[i][3], all_points[i][4]])
                    else:
                        print('Point is outside of raster.. Assigning previous height value of {}.'.format(current_gnd_z + altitude))
                        all_points_corrected.append(
                            [all_points[i][0], all_points[i][1], current_gnd_z + altitude, all_points[i][3], all_points[i][4]])
                    false_points += 1
                # Move to next checkpoint.
                current_loc[0] += move_per_step[0]
                current_loc[1] += move_per_step[1]

        # If leg is not active, fill in current altitude and move on to next leg.
        if all_points[i][3] == "false":
            print('Trajectory not active')
            all_points_corrected.append([all_points[i][0], all_points[i][1], current_alt, "false", all_points[i][4]])

    print('False Points: {}'.format(false_points))

    survey_file = open(survey_file_name, "w")
    survey_file.write('<?xml version="1.0"?>\n')
    survey_file.write('<document>\n')
    survey_file.write('    <survey name="{survey}" platform="{platform}" scanner="{scanner}" scene="{scene}">\n'
                      .format(survey=survey_name, platform=platform, scanner=scanner, scene=scene))
    # Write legs to file.
    for point in all_points_corrected:
        survey_file.write('        <leg stripId="{}">\n'.format(point[4]))
        survey_file.write('            <platformSettings  x="{x}" y="{y}" z="{z}" onGround="false" '
                          'movePerSec_m="{v}"/>\n'.format(x=point[0], y=point[1],
                                                          z=point[2], v=speed))
        survey_file.write('            <scannerSettings  active="{active_flag}" pulseFreq_hz="{pulse_freq}" scanAngle_deg="{scan_angle}" '
                          'scanFreq_hz="{scan_freq}" headRotatePerSec_deg="0.0" headRotateStart_deg="0.0" '
                          'headRotateStop_deg="0.0" trajectoryTimeInterval_s="0.05"/>\n'
                          .format(pulse_freq=pulse_freq, scan_angle=scan_angle, active_flag=point[3], scan_freq=scan_freq))
        survey_file.write('        </leg>\n')
    survey_file.write('    </survey>\n</document>')
    survey_file.close()