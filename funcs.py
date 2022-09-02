"""Library of functions for ALS point cloud parameter estimation
These functions can be combined in specific ways to estimate parameters from ALS point clouds. Some functions were
specifically designed to wok with the data generation script 'test_data_generator(_onlyrot).py' and the corresponding
script for merging the generated test data, 'merge_output.py'. See example workflow for parameter estimation in script
'param_estimate_multiple_strips.py'.
"""
import laspy
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import os
import pdal
import numpy.lib.recfunctions as rfn
import alphashape
import shapely
import open3d as o3d
import glob
import time
from scipy.spatial import ConvexHull


def minimum_bounding_rectangle(points):
    """
    Find the smallest bounding rectangle for a set of points.
    Returns a set of points representing the corners of the bounding box.

    :param points: an nx2 matrix of coordinates
    :rval: an nx2 matrix of coordinates
    """
    from scipy.ndimage.interpolation import rotate
    pi2 = np.pi/2.

    # get the convex hull for the points
    hull_points = points[ConvexHull(points).vertices]

    # calculate edge angles
    edges = np.zeros((len(hull_points)-1, 2))
    edges = hull_points[1:] - hull_points[:-1]

    angles = np.zeros((len(edges)))
    angles = np.arctan2(edges[:, 1], edges[:, 0])

    angles = np.abs(np.mod(angles, pi2))
    angles = np.unique(angles)

    # find rotation matrices
    # XXX both work
    rotations = np.vstack([
        np.cos(angles),
        np.cos(angles-pi2),
        np.cos(angles+pi2),
        np.cos(angles)]).T
#     rotations = np.vstack([
#         np.cos(angles),
#         -np.sin(angles),
#         np.sin(angles),
#         np.cos(angles)]).T
    rotations = rotations.reshape((-1, 2, 2))

    # apply rotations to the hull
    rot_points = np.dot(rotations, hull_points.T)

    # find the bounding points
    min_x = np.nanmin(rot_points[:, 0], axis=1)
    max_x = np.nanmax(rot_points[:, 0], axis=1)
    min_y = np.nanmin(rot_points[:, 1], axis=1)
    max_y = np.nanmax(rot_points[:, 1], axis=1)

    # find the box with the best area
    areas = (max_x - min_x) * (max_y - min_y)
    best_idx = np.argmin(areas)

    # return the best box
    x1 = max_x[best_idx]
    x2 = min_x[best_idx]
    y1 = max_y[best_idx]
    y2 = min_y[best_idx]
    r = rotations[best_idx]

    rval = np.zeros((4, 2))
    rval[0] = np.dot([x1, y2], r)
    rval[1] = np.dot([x2, y2], r)
    rval[2] = np.dot([x2, y1], r)
    rval[3] = np.dot([x1, y1], r)

    return rval


# Functions to calculate distance between to 3D points and flight alt from swath width
def calc_dist(p1, p2):
    """ Calculates euclidian distance between two 3D points
    :param p1: point 1 coords [x, y, z]
    :param p2: point 2 coords [x, y, z]
    :return: distance
    """
    xdist = p2[0] - p1[0]
    ydist = p2[1] - p1[1]
    zdist = p2[2] - p1[2]

    dist = pow((pow(xdist, 2) + pow(ydist, 2) + pow(zdist, 2)), 0.5)
    #print('Distance between {}, {}, {} and {}, {}, {} is {}'.format(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], dist))
    return dist


def calc_dist_2d(p1, p2):
    """ Calculates euclidian distance between two 2D points
        :param p1: point 1 coords [x, y]
        :param p2: point 2 coords [x, y]
        :return: distance
        """
    xdist = p2[0] - p1[0]
    ydist = p2[1] - p1[1]

    dist = pow((pow(xdist, 2) + pow(ydist, 2)), 0.5)
    return dist


def calc_centre(p1, p2):
    """Calculates the centre point between two 3D points
    :param p1: point 1 coords [x, y, z]
    :param p2: point 2 coords [x, y, z]
    :return: centre point coords [x, y, z]
    """
    x = p1[0]+(p2[0]-p1[0])/2
    y = p1[1]+(p2[1]-p1[1])/2
    z = p1[2]+(p2[2]-p1[2])/2

    return [x, y, z]


def calc_centre_2d(p1, p2):
    """Calculates the centre point between two 2D points
        :param p1: point 1 coords [x, y]
        :param p2: point 2 coords [x, y]
        :return: centre point coords [x, y]
        """
    x = p1[0]+(p2[0]-p1[0])/2
    y = p1[1]+(p2[1]-p1[1])/2

    return [x, y]


def calc_flight_alt(swath_width, scan_angle):
    """ Calculates the flight altitude of an airborne laser scanner based on the scan angle and the swath width
    :param swath_width: swath width in meters
    :param scan_angle: scan angle in degrees
    :return: altitude in meters
    """
    from math import tan, radians
    altitude = swath_width / 2 / tan(radians(scan_angle))
    return altitude


def read_point_cloud(filename, subsample=None, start=None, stop=None, num_returns=False, sort=True):
    """ Reads a las point cloud to a numpy array
    :param filename: path of the point cloud
    :param subsample: subsample step, every nth point will be stored in array
    :param start: starting point to read point cloud from, will read point cloud from nth point
    :param stop: read point cloud up until nth point
    :param num_returns: if true, omit reading return_number point record from point cloud
    :param sort: if true, sorts the array by gps time before returning, lowest to highest
    :return: numpy array with values read from point cloud
    """

    print('Reading file...: {}'.format(filename))
    try:
        with laspy.open(filename) as f:
            point_format = f.header.point_format
            scale = f.header.scale
            offset = f.header.offset
            point_num = f.header.point_count
            # print(f"Point format:       {point_format}")
            print(f"Number of points in file:   {point_num}")
            # print(f"Number of vlrs:     {len(f.header.vlrs)}")
            #print(list(point_format.standard_dimension_names))
            # print(f.header.offset)
    except Exception:
        print('The file "{}" could not be opened. This file will be skipped.'.format(filename))
        return

    if subsample != None:
        print("Subsampling to {} points.".format(int(point_num / subsample)))

    # Read file to numpy array
    las = laspy.read(filename)[start:stop:subsample]
    if num_returns:
        filtered_pc = np.array([las.X, las.Y, las.Z, las.gps_time])
    else:
        filtered_pc = np.array([las.X, las.Y, las.Z, las.gps_time, las.return_number])
    # Scale values in accordance with factors from header
    for i in range(3):
        filtered_pc[i, :] = filtered_pc[i, :] * scale[i] + offset[i]
    las = None

    if sort:
        # Sort points by gps time - lowest to highest
        filtered_pc = filtered_pc[:, filtered_pc[3, :].argsort()]

    return filtered_pc


def read_point_cloud_xyz(filename, subsample=None, start=None, stop=None, small_caps=False, first10=False, scale_offset=True):
    """Reads an xyz point cloud to a numpy array
    :param filename: path of the point cloud
    :param subsample: subsample step, every nth point will be stored in array
    :param start: starting point to read point cloud from, will read point cloud from nth point
    :param stop: read point cloud up until nth point
    :param small_caps: set to true if point records are stored as small caps. E.g. [las.x, las.y, las.z]
    :param first10: print first ten read points
    :param scale_offset: apply scale and offset calculation before returning
    :return: numpy array with values read from point cloud
    """

    print('Reading file...: {}'.format(filename))
    try:
        with laspy.open(filename) as f:
            point_format = f.header.point_format
            scale = f.header.scale
            offset = f.header.offset
            point_num = f.header.point_count
            # print(f"Point format:       {point_format}")
            print(f"Number of points in file:   {point_num}")
            # print(f"Number of vlrs:     {len(f.header.vlrs)}")
            # print(list(point_format.standard_dimension_names))
            print('Offset: {}'.format(f.header.offset))
            print('Scale: {}'.format(f.header.scale))
    except Exception:
        print('The file "{}" could not be opened. This file will be skipped.'.format(filename))
        return

    # Read file to numpy array
    las = laspy.read(filename)[start:stop:subsample]
    if small_caps:
        filtered_pc = np.array([las.x, las.y, las.z], dtype=np.float64)
    else:
        filtered_pc = np.array([las.X, las.Y, las.Z], dtype=np.float64)

    if first10:
        for p in range(10):
            print(filtered_pc[:, p])

    if scale_offset:
        # Scale values in accordance with factors from header
        for i in range(3):
            filtered_pc[i, :] = filtered_pc[i, :] * scale[i] + offset[i]
    las = None

    return filtered_pc


def estimate_alt_v_exp(pc_data, num_swaths, scan_angle, alpha_rad=200, visualize=False, visualize_swaths=False, ft=False):
    """Estimates flight altitude, flight speed and flight trajectory from point cloud
    :param pc_data: Point cloud as numpy array, may be subsampled. See function read_point_cloud for format.
    :param num_swaths: Number of full swaths present in the data #TODO: mention paper?
    :param scan_angle: Scan angle used by ALS scanner, in degrees
    :param alpha_rad: Alpharadius for calculation of alphashape
    :param visualize: Visualise histogram of points along GPStime flag
    :param visualize_swaths: Plot all swaths along x and y axis with swaths coloured individually
    :param ft: Distance unit is feet instead of meters flag
    :return: List of dicts with estimates for each swath
    """
    # Determine flight swaths
    sensitivity = 1

    # Add empty col for trajectory vals to pc array
    trajectory_col = np.empty((1, np.shape(pc_data)[1]), dtype=float)
    pc_data = np.vstack([pc_data, trajectory_col])

    # Iterate through gps time vals to determine gaps large enough to indicate new swath
    trajectory = 0
    swath_extents = [0]
    for i in range(np.shape(pc_data)[1] - 1):
        pc_data[4, i] = trajectory
        if (abs(pc_data[3, i] - pc_data[3, i + 1]) > sensitivity):
            swath_extents.append(i + 1)
            trajectory += 1
    swath_extents.append(np.shape(pc_data)[1] - 1)

    # Get sizes of swaths
    swath_sizes = [swath_extents[i + 1] - swath_extents[i] for i in range(len(swath_extents) - 1)]

    # Output info on swaths
    print('Number of swaths found: {}'.format(len(swath_extents) - 1))
    formatted_swaths = [[i + 1, swath_sizes[i]] for i in range(len(swath_sizes))]
    formatted_swaths.insert(0, ["swath nr.", "number of points"])
    print(tabulate(list(zip(*formatted_swaths))))

    if visualize:
        plt.figure(figsize=(18, 6))
        # Plot histogram of gps times, coloured by individual swaths
        for i in range(len(swath_extents) - 1):
            if i == swath_sizes.index(max(swath_sizes)):
                plt.hist(pc_data[3, swath_extents[i]:swath_extents[i + 1]], bins=10,
                         label='{} (largest)'.format(i + 1), color='red')
            else:
                plt.hist(pc_data[3, swath_extents[i]:swath_extents[i + 1]], bins=10, label=i + 1)  # , color='blue')

        #plt.title('{} - subsample\n{} Points'.format(os.path.basename(filename), np.shape(pc_data)[1]))
        plt.legend(title='Swath number')
        plt.ylabel('Amount of points')
        plt.xlabel('GPS time')
        plt.show()

    swath_sizes_sorted = swath_sizes.copy()
    swath_sizes_sorted.sort(reverse=True)
    print(swath_sizes_sorted)
    swath_nr=[]
    for i in range(num_swaths):
        swath_nr.append(swath_sizes.index(swath_sizes_sorted[i])+1)
        #swath_nr.append(swath_sizes.index(max(swath_sizes)) + 1)
        #swath_sizes.remove(max(swath_sizes))
    res = 1

    if visualize_swaths:
        fig1 = plt.figure(num=12)
        ax1 = fig1.add_subplot()

    results = []
    for counter, swth in enumerate(swath_nr):
        print('Analysing Swath Nr. {}'.format(swth))

        # Filter out swath
        swath = pc_data[:4, swath_extents[swth - 1] + 1:swath_extents[swth]].swapaxes(0, 1)
        print('Number of point in swath: {}'.format(np.shape(swath)))

        # Radius of circle with which alphashape is calculated
        alpha_val = 1 / alpha_rad
        print('Using alpha value of {}'.format(alpha_val))

        # Calculate alphashape
        alpha = alphashape.alphashape(swath[:, 0:2], alpha=alpha_val)

        # Calculate minimum bounding rectangle on alphashape
        alpha_vertices = np.asarray([[float(val) for val in xy.split(' ')] for xy in alpha.wkt.split('((')[1].split('))')[0].split(', ')])
        rect = minimum_bounding_rectangle(alpha_vertices)

        if visualize:
            # Add alphashape to graph above.
            fig = plt.figure(num=1)
            ax = fig.add_subplot()
            ax.scatter(swath[::res, 0], swath[::res, 1], c=swath[::res, 2], cmap='YlGnBu')
            ax.axis('equal')
            #plt.title('{} - Swath nr. {}\n{} Points'.format(os.path.basename(filename), swth, np.shape(swath)[0]))
            plt.xlabel('x')
            plt.ylabel('y')
            alpha_vertices = [xy.split(' ') for xy in alpha.wkt.split('((')[1].split('))')[0].split(', ')]
            # Empty arrays for bbox line. Will be updated after calculation of bbox.
            x = np.array([])
            y = np.array([])
            line1, = ax.plot(x, y, c='red')
            line1.set_xdata([float(xy[0]) for xy in alpha_vertices])
            line1.set_ydata([float(xy[1]) for xy in alpha_vertices])
            line1.set_label('Alphashape')
            print(rect)
            for number, rect_pt in enumerate(rect):
                print(rect[0], rect[1])
                ax.scatter(rect_pt[0], rect_pt[1], label='rect_pt_{}'.format(number))
            # fig.canvas.draw()
            plt.legend()
            plt.show()

        # Calculate lengths of minimum bounding rectangle on alphashape
        side1 = calc_dist_2d(rect[0], rect[3])
        side2 = calc_dist_2d(rect[0], rect[1])

        # Determine length of sides. Shorter side is set as width.
        if side1 < side2:
            print('Side 1 is width')
            width = side1
            length = side2
            trajectory_p1 = calc_centre_2d(rect[0], rect[3])
            trajectory_p2 = calc_centre_2d(rect[1], rect[2])
        else:
            print('Side 2 is width')
            width = side2
            length = side1
            trajectory_p1 = calc_centre_2d(rect[0], rect[1])
            trajectory_p2 = calc_centre_2d(rect[2], rect[3])

        if ft==True:
            width=width*0.3048
            length=length*0.3048

        if visualize_swaths:# and swth in [4, 3, 9, 8, 6, 5, 10, 11]:
            #colours = ['lightsalmon', 'dodgerblue', 'dodgerblue', 'lightsalmon', 'lightsalmon'
                #, 'dodgerblue', 'dodgerblue', 'lightsalmon', 'dodgerblue', 'lightsalmon', 'lightsalmon', 'dodgerblue']
            #colours = ['lightsalmon', 'dodgerblue', 'lightsalmon', 'dodgerblue', 'lightsalmon', 'dodgerblue', 'lightsalmon', 'dodgerblue', 'lightsalmon', 'dodgerblue',
                       #'lightsalmon', 'dodgerblue','lightsalmon', 'dodgerblue','lightsalmon', 'dodgerblue','lightsalmon', 'dodgerblue']
            colour = 'lightsalmon' if swth%2==0 else 'dodgerblue'
            print('Colour is {}'.format(colour))
            print(alpha.wkt)
            alpha_vertices = [xy.split(' ') for xy in alpha.wkt.split('((')[1].split('))')[0].split(', ')]
            ax1.scatter(swath[::res, 0], swath[::res, 1], alpha=0.5, c=colour, s=1)
            # Plot first and last points of each swath
            #ax1.scatter(swath[0, 0], swath[0, 1], s=5, label='First Point S{}'.format(swth))
            #ax1.scatter(swath[-1, 0], swath[-1, 1], s=5, label='Last Point S{}'.format(swth))
            #ax1.plot([trajectory_p1[0], trajectory_p2[0]], [trajectory_p1[1], trajectory_p2[1]], label='Swath {}'.format(swth))
            #ax1.plot([float(xy[0]) for xy in alpha_vertices], [float(xy[1]) for xy in alpha_vertices], alpha=0.2, c='black')

        print('Approximate swath width: {:0.2f}m'.format(width))
        print('Scan angle: {}'.format(scan_angle))

        alt = calc_flight_alt(width, scan_angle)
        print('Likely flight altitude: {:0.2f}m'.format(alt))

        # Swath length:
        print('Approximate swath length: {:0.2f}m'.format(length))

        # Calculate flight speed
        flight_time = swath[-1, 3] - swath[0, 3]
        print('Approximate flight time: {:0.2f}s'.format(flight_time))

        flight_speed = length / flight_time
        print('Approximate flight speed: {:0.2f}m/s\n\n'.format(flight_speed))

        results.append({'num_points' : np.shape(swath)[0], 'altitude' : alt, 'swath_width' : width, 'flight_time' : flight_time, 'flight_speed' : flight_speed, 'trajectory' : [trajectory_p1, trajectory_p2]})

    if visualize_swaths:
        # Turn off tick labels
        #ax1.set_yticklabels([])
        #ax1.set_xticklabels([])
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.legend()
        plt.show()

    return results


def filter_swaths(pc_data, output_pc, num_swaths):
    """Creates copy of point cloud with only n largest swaths
    :param pc_data: Input point cloud as numpy array (see read_point_cloud function)
    :param output_pc: Filename and path of output point cloud
    :param num_swaths: Number of swaths to be left in copy
    :return: nothing
    """
    # Determine flight swaths
    sensitivity = 1

    # Add empty col for trajectory vals to pc array
    trajectory_col = np.empty((1, np.shape(pc_data)[1]), dtype=float)
    pc_data = np.vstack([pc_data, trajectory_col])

    # Iterate through gps time vals to determine gaps large enough to indicate new swath
    trajectory = 0
    swath_extents = [0]
    for i in range(np.shape(pc_data)[1] - 1):
        pc_data[4, i] = trajectory
        if (abs(pc_data[3, i] - pc_data[3, i + 1]) > sensitivity):
            swath_extents.append(i + 1)
            trajectory += 1
    swath_extents.append(np.shape(pc_data)[1] - 1)

    # Get sizes of swaths
    swath_sizes = [swath_extents[i + 1] - swath_extents[i] for i in range(len(swath_extents) - 1)]

    # Output info on swaths
    print('Number of swaths found: {}'.format(len(swath_extents) - 1))
    formatted_swaths = [[i + 1, swath_sizes[i]] for i in range(len(swath_sizes))]
    formatted_swaths.insert(0, ["swath nr.", "number of points"])
    print(tabulate(list(zip(*formatted_swaths))))

    swath_sizes_sorted = swath_sizes.copy()
    swath_sizes_sorted.sort(reverse=True)
    print(swath_sizes_sorted)
    swath_nr = []
    for i in range(num_swaths):
        swath_nr.append(swath_sizes.index(swath_sizes_sorted[i]) + 1)
        # swath_nr.append(swath_sizes.index(max(swath_sizes)) + 1)
        # swath_sizes.remove(max(swath_sizes))
    res = 1

    outlas = np.array([[], [], [], [], []])

    results = []
    for counter, swth in enumerate(swath_nr):
        print('Analysing Swath Nr. {}'.format(swth))

        # Filter out swath
        swath = pc_data[:, swath_extents[swth - 1] + 1:swath_extents[swth]]
        print('Number of point in swath: {}'.format(np.shape(swath)))
        outlas = np.concatenate((outlas, swath), axis=1)

    print(np.shape(outlas))

    # 1. Create a new header
    header = laspy.LasData(header)
    header.scales = np.array([0.1, 0.1, 0.1])

    # 3. Create a LasWriter and a point record, then write it
    with laspy.open(output_pc, mode="w", header=header) as writer:
        point_record = laspy.ScaleAwarePointRecord.zeros(outlas.shape[0], header=header)
        point_record.x = outlas[:3, 0]
        point_record.y = outlas[:3, 1]
        point_record.z = outlas[:3, 2]

        writer.write_points(point_record)


def sort_gpstime(las_sort, las_file):
    """Sort a point cloud by GPS time using lassort
    :param las_sort: filepath of lassort.exe
    :param las_file: las file to sort
    :return: nothing
    """
    os.popen('{} {} -gps_time'.format(las_sort, las_file))


def estimate_prf(pc_data, show_values=False):
    """ Estimates the pulse repetition frequency from a set of points from a point cloud
    :param pc_data: consecutive points from a point cloud, with accurate GPS time values
    :param show_values: true to print all three estimates
    :return: estimate of pulse repetition frequency
    """
    import statistics

    print('Estimating PRF...')
    # Leave only first return points.
    pc_data = pc_data[:, pc_data[4, :] < 1.1]
    # First estimate from total duration and total number of pulses
    total_duration = pc_data[3, -1] - pc_data[3, 1]
    total_pulses = np.shape(pc_data)[1]
    pulse_freq = total_pulses / total_duration

    # Estimate from actual distances between points.
    # If several point in a row have the same value, the pulse count is increased
    # and the consecutive gps time gap is divided by number of same pulses
    num_pulses = 1
    freqs = []
    for i in range(np.shape(pc_data)[1] - 1):
        # print("{} {}".format(pc_data[3, i], pc_data[3, i + 1]))
        if pc_data[3, i] != pc_data[3, i + 1]:
            freq = 1 / (pc_data[3, i + 1] - pc_data[3, i]) * num_pulses
            num_pulses = 1
            freqs.append(freq)
        else:
            num_pulses += 1

    median_val = statistics.median(freqs)

    if show_values:
        print('Estimate from total flight duration:')
        print(pulse_freq)
        print('Estimate per point difference:')
        print(freqs)
        print('Median of estimates per point difference:')
        print(median_val)

    #mean_val = statistics.mean(freqs)
    #mode_val = statistics.mode(freqs)
    print('Estimated PRF - {}'.format(median_val))
    return median_val


def read_true_values(metadata_file):
    """ Reads the metadata values from a metadata file
    :param metadata_file: A metadata as written by
    :return:
    """
    true_vals = {}

    metadata_true = open(metadata_file, "r")
    lines = metadata_true.readlines()
    metadata_true.close()

    true_vals['scan_angle'] = float(lines[-3].split(": ")[1])
    true_vals['altitude'] = float(lines[-1].split(": ")[1])
    true_vals['speed'] = float(lines[-2].split(": ")[1])
    true_vals['PRF'] = float(lines[-5].split(": ")[1])
    true_vals['traj_num'] = int(lines[-7].split('tiny_flight_')[1].split('.shp')[0])
    true_vals['scanner'] = lines[-6].split('#')[1].split('\n')[0]

    return true_vals


def write_results_tab(results_file, filename, true_vals, results, prf):
    """Writes results of parameter estimation of a point cloud generated using merge_output.py and test_data_generator.py to a csv-type output file. Format:
    dataset_num;scanner;shp_num;true_scan_angle;true_altitude;true_speed;true_prf;[swath_width altitude flight_speed]x3;prf
    :param results_file: Location and name of results file
    :param filename: Location and name of point cloud, generated using merge_output.py and test_data_generator.py
    :param true_vals: True value object as generated by read_true_values function
    :param results: Results object as generated by estimate_alt_v function
    :param prf: Estimate of pulse repetition frequency, calculated using estimate_prf function
    :return: nothing
    """

    dataset_num = filename.split('merged_')[1].split('.la')[0]
    with open(results_file, 'a') as f:
        f.write('{};{};{};{};{};{};{};'.format(dataset_num, true_vals['scanner'], true_vals['traj_num'], true_vals['scan_angle'], true_vals['altitude'], true_vals['speed'], true_vals['PRF']))

        for i, swath in enumerate(results):
            f.write('{:.2f};{:.2f};{:.2f};{};'.format(swath['swath_width'], swath['altitude'], swath['flight_speed'], str(swath['trajectory'])))

        f.write('{:.2f}\n'.format(prf))


def write_results_tab_no_sim(results_file, true_vals, results, prf):
    """Writes results of parameter estimation to an output file. Format:
    dataset_num;scanner;shp_num;true_scan_angle;true_altitude;true_speed;true_prf;[swath_width altitude flight_speed]x3;prf
    :param results_file: Location and name of results file
    :param true_vals: True value object as generated by read_true_values function
    :param results: Results object as generated by estimate_alt_v function
    :param prf: Estimate of pulse repetition frequency, calculated using estimate_prf function
    :return: nothing
    """

    with open(results_file, 'a') as f:
        f.write('{};{};{};{};{};{};\n'.format(true_vals['scanner'], true_vals['traj_num'], true_vals['scan_angle'], true_vals['altitude'], true_vals['speed'], true_vals['PRF']))

        for i, swath in enumerate(results):
            f.write('{};{:.2f};{:.2f};{:.2f};{};\n'.format(swath['num_points'], swath['swath_width'], swath['altitude'], swath['flight_speed'], str(swath['trajectory'])))

        f.write('{:.2f}\n'.format(prf))


def run_sim(survey_name):
    """ Run a HELIOS++ simulation
    :param survey_name: absolute or relative path to a HELIOS++ survey file
    :return: output filepath of survey
    """
    import pyhelios

    pyhelios.loggingVerbose()
    pyhelios.setDefaultRandomnessGeneratorSeed("123")

    sim = pyhelios.Simulation(
        survey_name,
        'assets/',
        'output/',
        0,  # Num Threads
        True,  # LAS output
        False,  # LAS1.0 output
        False,  # ZIP output
    )

    # Enable final output.
    sim.finalOutput = True

    # Set sim frequency.
    sim.simFrequency = 1000

    # Load survey file. Further configuration of survey possible.
    sim.loadSurvey(
        True,  # Leg Noise Disabled FLAG
        False,  # Rebuild Scene FLAG
        False,  # Write Wavef orm FLAG
        False,  # Calc Echowidth FLAG
        False,  # Full Wave Noise FLAG
        True  # Platform Noise Disabled FLAG
    )

    print('Simulation has started!\nSurvey Name: {survey_name}\n{scanner_info}'.format(
        survey_name=sim.getSurvey().name,
        scanner_info=sim.getScanner().toString()))
    sim.start()

    output = sim.join()

    output_folder = output.filepath

    return output_folder


def intersect_poly(wkt1, wkt2):
    """ Calculates the intersection of two wkt polygons
    :param wkt1: a wkt polygon
    :param wkt2: a wkt polygon
    :return: wkt polygon of intersection
    """
    from osgeo import ogr

    poly1 = ogr.CreateGeometryFromWkt(wkt1)
    poly2 = ogr.CreateGeometryFromWkt(wkt2)

    intersection = poly1.Intersection(poly2)

    return intersection.ExportToWkt()


def get_alphashapes(pc_data, num_swaths, alpha_rad):
    """ Calculates alphashapes of swaths in ALS point cloud
    Usage example: intersect_poly.py
    :param pc_data: Point cloud data, as numpy array (see function read_point_cloud)
    :param num_swaths: Number of full swaths present in data #TODO Reference to thesis? Allowed?
    :param alpha_rad: Alpha radius with which to calculate alphashapes
    :return: [alphashape as wkt, number of points in swath, estimated swath width]
    """
    # Determine flight swaths
    sensitivity = 1

    # Add empty col for trajectory vals to pc array
    trajectory_col = np.empty((1, np.shape(pc_data)[1]), dtype=float)
    pc_data = np.vstack([pc_data, trajectory_col])

    # Iterate through gps time vals to determine gaps large enough to indicate new swath
    trajectory = 0
    swath_extents = [0]
    for i in range(np.shape(pc_data)[1] - 1):
        pc_data[4, i] = trajectory
        if (abs(pc_data[3, i] - pc_data[3, i + 1]) > sensitivity):
            swath_extents.append(i + 1)
            trajectory += 1
    swath_extents.append(np.shape(pc_data)[1] - 1)

    # Get sizes of swaths
    swath_sizes = [swath_extents[i + 1] - swath_extents[i] for i in range(len(swath_extents) - 1)]

    # Output info on swaths
    print('Number of swaths found: {}'.format(len(swath_extents) - 1))
    formatted_swaths = [[i + 1, swath_sizes[i]] for i in range(len(swath_sizes))]
    formatted_swaths.insert(0, ["swath nr.", "number of points"])
    print(tabulate(list(zip(*formatted_swaths))))

    swath_sizes_sorted = swath_sizes.copy()
    swath_sizes_sorted.sort(reverse=True)
    print(swath_sizes_sorted)
    swath_nr = []
    for i in range(num_swaths):
        swath_nr.append(swath_sizes.index(swath_sizes_sorted[i]) + 1)
        # swath_nr.append(swath_sizes.index(max(swath_sizes)) + 1)
        # swath_sizes.remove(max(swath_sizes))

    res = 1
    print(swath_nr)
    alphas = []
    sizes = []
    widths = []
    for swth in swath_nr:
        print('Analysing Swath Nr. {}'.format(swth))

        # Filter out swath
        swath = pc_data[:4, swath_extents[swth - 1] + 1:swath_extents[swth]].swapaxes(0, 1)

        # Radius of circle with which alphashape is calculated
        alpha_val = 1 / alpha_rad
        print('Using alpha value of {}'.format(alpha_val))

        # Calculate alphashape
        alpha = alphashape.alphashape(swath[:, 0:2], alpha=alpha_val)
        alphas.append(alpha.wkt)
        sizes.append(swath.shape[0])

        alpha_vertices = np.asarray(
            [[float(val) for val in xy.split(' ')] for xy in alpha.wkt.split('((')[1].split('))')[0].split(', ')])
        rect = minimum_bounding_rectangle(alpha_vertices)
        # Calculate lengths of minimum bounding rectangle on alphashape
        side1 = calc_dist_2d(rect[0], rect[3])
        side2 = calc_dist_2d(rect[0], rect[1])

        # Determine length of sides. Shorter side is set as width.
        if side1 < side2:
            print('Side 1 is width')
            width = side1
            length = side2
            tp1 = calc_centre_2d(rect[0], rect[3])
            tp2 = calc_centre_2d(rect[1], rect[2])
        else:
            print('Side 2 is width')
            width = side2
            length = side1
            tp1 = calc_centre_2d(rect[0], rect[1])
            tp2 = calc_centre_2d(rect[2], rect[3])
        widths.append(width)

    return alphas, sizes, widths


def visualize_wkt(wkts, labels, title, data=None):
    """Visualises multiple wkt polygons
    :param wkts: wkt polygons as list of wkt strings
    :param labels: list of strings containing labels for each polygon
    :param title: title of figure (str)
    :param data: point cloud to visualise on top of polygons (numpy array with x y z values [[x, y, z, ..., ...],[...],[...]])
    :return: nothing
    """
    from matplotlib.patches import Ellipse, Polygon

    # Scatterplot of swath
    fig = plt.figure(num=1)
    ax = fig.add_subplot()

    # Turn off tick labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    # plt.title('Swath nr. {}\n{} Points'.format(swath_nr, np.shape(swath)[0]))
    # plt.title('Alphashape of Single Flight Swath')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(title)
    plt.axis('equal')
    if data:
        ax.scatter(data[:,0], data[:, 1], c = data[:,2], cmap='YlGnBu', label='Flight Swath', s=2.5, alpha=0.5)
    #line3 = ax.plot([tp1[0], tp2[0]], [tp1[1], tp2[1]], '--', c='red', linewidth=2, label='Flight Trajectory')
    for i, wkt in enumerate(wkts):
        print('Wkt {}: '.format(i) + wkt)
        try:
            alpha_vertices = [xy.split(' ') for xy in wkt.split('((')[1].split('))')[0].split(', ')]
            if i==1:
                ax.plot([float(xy[0]) for xy in alpha_vertices], [float(xy[1]) for xy in alpha_vertices], label=labels[i])
            else:
                ax.plot([float(xy[0]) for xy in alpha_vertices], [float(xy[1]) for xy in alpha_vertices], '--', label=labels[i])
        except ValueError:
            alpha_vertices = [xy.split(' ') for xy in wkt.split('((')[1].split('))')[0].split(',')]
            #print(alpha_vertices)
            ax.add_patch(Polygon([(xy[0], xy[1]) for xy in alpha_vertices], facecolor='green', alpha=0.3, label=labels[i]))
            #ax.plot([float(xy[0]) for xy in alpha_vertices], [float(xy[1]) for xy in alpha_vertices], '--', label='Wkt {}'.format(i))

    plt.legend(loc='best')
    plt.show()


def write_overall_results(results_file, devs):
    """
    :param results_file:
    :param devs:
    :return:
    """
    with open(results_file, 'a') as f:
        f.write('\nOverall deviation for all files:\n')
        f.write(
            '----------------------------------------------------------------------------------------------------\n')
        f.write('Mean altitude deviation: {}\n'.format(statistics.mean(devs[0])))
        f.write('Mean speed deviation: {}\n\n\n'.format(statistics.mean(devs[1])))


def wkt_translate(wkt, x, y):
    """Add a translation to a wkt polygon
    :param wkt: wkt polygon string
    :param x: x offset
    :param y: y offset
    :return: modified wkt polygon string
    """
    verts = [xy.split(' ') for xy in wkt.split('((')[1].split('))')[0].split(', ')]
    new_verts = 'POLYGON (('
    for xy in verts:
        new_verts += '{} {}, '.format(float(xy[0]) + x, float(xy[1]) + y)
    new_verts = new_verts[:-2]
    new_verts += '))'
    return new_verts


def check_pointrecords(file, sample_size):
    """Checks which columns in point clouds are filled with unique values, prints result
    :param file: point cloud (las/laz)
    :param sample_size: how many point to read from file
    :return:
    """
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


def frtxt(in_file, find, replace):
    """Replace text in a file with some other text
    :param in_file: File to edit
    :param find: text to replace
    :param replace: replacement text
    :return: -
    """
    with open(in_file, 'r') as file:
        # Read in the file
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace(find, replace)

    # Write the file out again
    with open(in_file, 'w') as file:
        file.write(filedata)