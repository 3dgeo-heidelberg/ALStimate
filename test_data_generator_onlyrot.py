# Author: Mark Searle, 3DGeo Research Group Heidelberg
# Contact: mark.searle@stud.uni-heidelberg.de
# -*- coding: utf-8 -*-
"""HELIOS++ test data generation (only rotating polygon mirror scanners)
This script allows the construction and execution of a series of HELIOS++ ALS simulations with a random rotating polygon
mirror scanner and random settings for the scanner. The trajectories are defined by shapefiles located in the 'input_shp'
directory andthe altitude is set randomly, and held over the 'input_raster' DEM with a max deviation of
'height_variation'. The scenefile is defined by the variable 'scene'. A metadata file is automatically written to the
output location indicating all the parameters used in the survey.
"""
import random
import os

helios_path = ''
lastools = r"G:\LAStools\bin\lasmerge.exe"
outdir = r"G:/ws_21_22/bachelor_arbeit/data/generated_data/test_data3"
height_variation = 60
num_surveys = 20
#l = [1, 10, 15, 17, 19, 25, 9, 0, 11, 18, 20, 22, 23, 4]

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


for i in range(num_surveys):
    # Use flight_planner tool from helios to generate the legs
    # Use XML parse to select a scanner
    # Use random.randint for scan angle, altitude, velocity
    survey_file_name = 'data_gen_{}.xml'.format(i)#i
    scene = 'data/scenes/demo/hd_demo.xml#hd_demo'

    scanner_selector = random.randint(0,1)
    #pulse_freq_selector = random.randint(0,2)
    #scan_freq_selector = random.randint(0,2)
    #scan_angle_selector = random.randint(0,1)

    # List with available scanners. Format: [name, [min_PRF, max_PRF], [min_scan_freq, max_scan_freq], [min_scan_angle, max_scan_angle]]
    scanner_repo = [['riegl_lms-q780', [100000, 300000], [70, 200], [15, 30]],
                    ['riegl_lms-q560', [100000, 200000], [80, 160], [15, 22]]]

    # Generate random survey parameters.
    platform = "data/platforms.xml#sr22"
    speed = random.randint(45, 70)  # in m/s - roughly estimated from point cloud
    altitude = random.randint(1000, 2000)
    scanner = "data/scanners_als.xml#" + scanner_repo[scanner_selector][0]
    pulse_freq = random.randint(scanner_repo[scanner_selector][1][0], scanner_repo[scanner_selector][1][1])
    scan_freq = random.randint(scanner_repo[scanner_selector][2][0], scanner_repo[scanner_selector][2][0])
    scan_angle = random.randint(scanner_repo[scanner_selector][3][0], scanner_repo[scanner_selector][3][1])
    survey_name = survey_file_name.split('.')[0]

    survey_file = open(survey_file_name, "w")
    survey_file.write('<?xml version="1.0"?>\n')
    survey_file.write('<document>\n')
    survey_file.write('    <survey name="{survey}" platform="{platform}" scanner="{scanner}" scene="{scene}">\n'
                      .format(survey=survey_name, platform=platform, scanner=scanner, scene=scene))

    from osgeo import ogr, gdal, osr

    input_shp = r"flight_paths\tiny_flight_{}.shp".format(random.randint(1,5))
    print(input_shp)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(input_shp, 0)
    layer = dataSource.GetLayer()

    input_raster = r"data\sceneparts\tiff\dem_hd.tif"
    raster = gdal.Open(input_raster)
    proj = osr.SpatialReference(wkt=raster.GetProjection())
    print(proj.GetAttrValue('AUTHORITY',1))
    geotransform = raster.GetGeoTransform()
    raster_b1 = raster.GetRasterBand(1)


    # Get trajectory from vector file.
    all_points = []
    for feature in layer:
        geom = feature.GetGeometryRef()
        for i in range(0, geom.GetPointCount()):

            # GetPoint returns a tuple not a Geometry
            pt = geom.GetPoint(i)
            print("%d). POINT (%f %f)" % (i, pt[0], pt[1]))

            '''if isinstance(v, QgsPoint):
                v_xy = QgsPointXY(v)
            v_tranformed = crs_transform.transform(v_xy)'''

            # For last vertice, set 'active' flag to false.
            if i+1 != geom.GetPointCount():
                all_points.append([pt[0], pt[1], altitude, "true"])
                print("active")
            else:
                all_points.append([pt[0], pt[1], altitude, "false"])
                print("inactive")

    # Update trajectory values with values over raster
    all_points_corrected = []
    # Update altitude.
    for i in range(len(all_points)):

        # Only if leg is active..:
        if all_points[i][3] == "true":
            print('Trajectory active')
            #print('Moving between points {} and {}'.format(all_points[i], all_points[i+1]))
            # Get raster height
            pixel, line = world2Pixel(geotransform, all_points[i][0], all_points[i][1])
            print(pixel, line)
            value = raster_b1.ReadAsArray(pixel, line, 1, 1)[0, 0]
            # Add raster height to flight alt at beginning of trajectory.
            all_points[i][2] += value
            # Add point entry to list with new corrected altitude values.
            all_points_corrected.append([all_points[i][0], all_points[i][1], all_points[i][2], all_points[i][3]])
            # For every leg, check for altitude 100 times.
            frequency = 100
            current_loc = all_points[i]
            current_alt = all_points[i][2]
            # Distance to be moved between each altitude check.
            move_per_step = [(all_points[i + 1][0] - all_points[i][0]) / frequency,
                             (all_points[i + 1][1] - all_points[i][1]) / frequency]
            print("Move per step is {}".format(move_per_step))
            # Perform altitude checks, update altitude value if necessary.
            for j in range(frequency):
                pixel, line = world2Pixel(geotransform, current_loc[0], current_loc[1])
                current_gnd_z = raster_b1.ReadAsArray(pixel, line, 1, 1)[0, 0]
                distance_to_optimal = current_gnd_z + altitude - current_alt
                if abs(distance_to_optimal) > height_variation:
                    # If altitude is out of range:
                    print('new alt: {}'.format(current_gnd_z + altitude))
                    # Update value in survey leg.
                    all_points_corrected.append(
                        [current_loc[0], current_loc[1], current_gnd_z + altitude, "true"])
                    # Update value for upcoming checks.
                    current_alt = current_gnd_z + altitude
                # Move to next checkpoint.
                current_loc[0] += move_per_step[0]
                current_loc[1] += move_per_step[1]

        # If leg is not active, fill in current altitude and move on to next leg.
        if all_points[i][3] == "false":
            print('Trajectory not active')
            all_points_corrected.append([all_points[i][0], all_points[i][1], current_alt, "false"])


    # Write legs to file.
    for point in all_points_corrected:
        survey_file.write('        <leg>\n')
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

    import sys
    import pyhelios

    pyhelios.loggingVerbose2()
    pyhelios.setDefaultRandomnessGeneratorSeed("123")

    sim = pyhelios.Simulation(
        str(survey_file_name),
        'assets/',
        'output/',
        0,  # Num Threads
        False,  # LAS output
        False,  # LAS1.0 output
        False,  # ZIP output
    )

    # Enable final output.
    sim.finalOutput = True

    # Set sim frequency.
    sim.simFrequency = 0

    # Load survey file. Further configuration of survey possible.
    sim.loadSurvey(
        True,  # Leg Noise Disabled FLAG
        False,  # Rebuild Scene FLAG
        False,  # Write Wavef orm FLAG
        False,  # Calc Echowidth FLAG
        False,  # Full Wave Noise FLAG
        True  # Platform Noise Disabled FLAG
    )

    sim.start()
    print('Simulation has started!\nSurvey Name: {survey_name}\n{scanner_info}'.format(
        survey_name=sim.getSurvey().name,
        scanner_info=sim.getScanner().toString())
    )

    sim.join()

    # Find survey output directory (to inform the user).
    output_dir = os.path.join("output/", survey_name)
    # Find latest folder in survey output directory
    all_survey_outputs = [os.path.join(output_dir, d) for d in os.listdir(output_dir) if
                          os.path.isdir(os.path.join(output_dir, d))]
    latest_survey_output = max(all_survey_outputs, key=os.path.getmtime)

    # Write metadata file with all information from the survey
    metadata = os.path.join(latest_survey_output, 'metadata.txt')

    with open(metadata, 'w') as f:
        f.write('Survey no.: {}\n'.format(survey_file_name))
        f.write('trajectory: {}\n'.format(input_shp))
        f.write('scanner: {}\n'.format(scanner))
        f.write('pulse_freq: {}\n'.format(pulse_freq))
        f.write('scan_freq: {}\n'.format(scan_freq))
        f.write('scan_angle: {}\n'.format(scan_angle))
        f.write('platform_velocity: {}\n'.format(speed))
        f.write('platform_altitude(a.g.l.): {}\n'.format(altitude))