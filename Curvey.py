"""
Created on Wed Oct 08 08:50:44 2014

@author: Nicky Jay feat Angus M.
Welcome to Curvey
"""
import numpy as np
from numpy.linalg import eig, inv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ConfigParser import SafeConfigParser
import ast
from scipy import optimize
from PIL import Image
from scipy.ndimage.measurements import label, center_of_mass
import scipy.ndimage as ndi
import os

def fit_ellipse(x,y):
    """Fancy stuff to fit an ellipse to selected points. Using eigenvectors.
       works out the mean of the ellipse offset and subtracts it because this 
       can only work out ellipse points within a certain range. This is added
       back latter"""
    x = x[:,np.newaxis].astype(float)
    y = y[:,np.newaxis].astype(float)
    xm = x.mean()
    ym = y.mean()
    x -= xm
    y -= ym
    
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = np.real_if_close(V[:,n])
    return a, xm, ym


def ellipse_centre(a,z,xm,ym):
    """Finds centre of the fitted ellipse. Can also show this on the fitted
       ellipse... like in my dissertation. This is where ellipse mean offset
       is added back..."""
    #plt.plot(xm, ym, 'rx', ms=5)
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0 = (c*d-b*f)/num + xm
    y0 = (a*f-b*d)/num + ym
    z0 = z
    return np.array([x0,y0,z0])

def show_fitted_ellipse(coord_list, xx, yy):
    """shows the fitted ellipse for a single slice through the IHCs... only 
       used once for my dissertation so it looks pretty ugly..."""
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.scatter(coord_list[:8,2], coord_list[:8,0], c ='r')
    
    wdir = os.path.normpath(r'C:\Users\DELL\Desktop')
    mtiff_file = "P6_1_myo.tif"

    print "Importing data"

    im = os.path.join(wdir, mtiff_file)
    nimg = 120
    img = Image.open(im)
    img.seek(0)
    nx, ny = img.size

    print "Allocating array"

    big_arr = np.zeros([nimg, nx, ny], dtype=np.uint8)

    for i in range(nimg):
        img.seek(i)
        big_arr[i] = np.array(img, dtype=np.uint8)
    
    plt.plot(yy,xx, color = 'red')
    ax.scatter(centre[1], centre[0], c = 'r', marker= 'x')
    plt.imshow(big_arr[:,200, 0:400])
    
def old_ellipse_centre(a):
    """Finds the centre of the ellipse without adding on the ellipse offsets"""
    #plt.plot(xm, ym, 'rx', ms=5)
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0 = (c*d-b*f)/num 
    y0 = (a*f-b*d)/num
    
    return np.array([x0,y0])

def ellipse_angle_of_rotation(a):
    """Self explanatory really..."""
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 
    5*np.arctan(2*b/(a-c))

def ellipse_axis_length(a):
    """Also self explanatory"""
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=(np.sqrt(up/down1)) * 0.2
    res2=(np.sqrt(up/down2)) * 0.0689988270199406
    
    return res1/ res2

def load_data(parser,config_file):
    """Loads the coordinates which were clicked on in Clicky 2.0 which are
       used to form the IHC central axes. Could make it read the number of cells
    """
    parser.read(config_file)
    num_cells = 4
    multi_cell_list = []
    for i in xrange(num_cells):
        top_list = parser.get('Cell %s' % int(i + 1), 'Top')
        top_list = ast.literal_eval(top_list)
        multi_cell_list.append(top_list)
    return np.array(multi_cell_list)
    
def fit_many_ellipses(coord_list,z_list):
    """Fits all the ellipses for all slices levels and number of IHCs cells.
       Finds all the ellipse centres too... Uses the functions mentioned above
    """
    big_fitted_list = []
    big_xm_list = []
    big_ym_list = []
    big_centre_list = []
    for i in xrange(len(coord_list)):
        fitted_list = []
        xm_list = []
        ym_list = []
        """This is why you have to click 8 times for each ellipse in clicky 2.0... 
           It's ugly but it works... Certainly has room for improvement.
           8 is a good number of clicks anyways"""
        xarr = np.array([coord_list[i,:8,0],coord_list[i,8:16,0],coord_list[i,16:24,0],
                         coord_list[i,24:32,0],coord_list[i,32:40,0]])
       
        yarr = np.array([coord_list[i,:8,2],coord_list[i,8:16,2],coord_list[i,16:24,2],
                         coord_list[i,24:32,2],coord_list[i,32:40,2]])
                         
        #plt.scatter(coord_list[2,16:24,0],coord_list[2,16:24,2])
        for ii in range(5):
            a, xm, ym = fit_ellipse(xarr[ii],yarr[ii])
            fitted_list.append(a)
            xm_list.append(xm)
            ym_list.append(ym)
        big_fitted_list.append(fitted_list)
    
        big_xm_list.append(xm_list)
        big_ym_list.append(ym_list)
    
    for i in xrange(len(big_fitted_list)):
        centre_list = []
        for ii in xrange(len(big_fitted_list[i])):
            centre_list.append(ellipse_centre(big_fitted_list[i][ii],z_list[ii],
                                              big_xm_list[i][ii],big_ym_list[i][ii]))
        big_centre_list.append(centre_list)
        
    return big_centre_list, big_fitted_list, big_xm_list, big_ym_list

def get_ellipse_axes(big_fitted_list):
    for i in xrange(len(big_fitted_list)):
        print " "
        for ii in xrange(len(big_fitted_list[0])):
            print ellipse_axis_length(big_fitted_list[i][ii])

def curve_func(data, a, b, c):
    """input function for the fit_curved_axis function below :) """
    return a*data**2 + b*data + c

def fit_curved_axis(data):
    """Fits polynomials in both the x and y dimensions to make some snazzy curved
       central axes for the IHCs.
    """
    guess = (1,1,1)
    x_params_list = []
    y_params_list = []
    for i in xrange(len(data)):
        x = data[i,:,0] * 0.2
        y = data[i,:,1] * 0.0689988270199406
        z = data[i,:,2]
        x_params, pcov = optimize.curve_fit(curve_func, z, x, guess)
        y_params, pcov = optimize.curve_fit(curve_func, z, y, guess)
        x_params_list.append(x_params)
        y_params_list.append(y_params)
    return x_params_list, y_params_list
    
    
R = np.arange(0,2*np.pi, 0.1)
    
def find_ellipse_params(fitted_list, xm_list, ym_list):
    """Finds ellipse paramaters, useful if you're plotting the ellipses...
       Not so useful for other things. Not really used but nice to have
    """
    centre = old_ellipse_centre(fitted_list)
    phi = ellipse_angle_of_rotation(fitted_list)
    axes = ellipse_axis_length(fitted_list)
                                            
    centre[0] = centre[0] + xm_list[0] 
    centre[1] = centre[1] + ym_list[0]

    a, b = axes
    xx = centre[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
    yy = centre[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)
    return xx, yy
    
def show_fitted_curve(line_x, line_y, line_z, data, pt_coords):
    """Shows the fitted central axes. Also shows the puncta locations"""
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    plt.xlim(0,35)
    plt.ylim(0,35)
    plt.ylabel("y")
    plt.xlabel("x")
    
    for i in xrange(len(data)):
        xd, yd, zd  = (data[i]).T
        ellipse_x = xd * 0.2
        ellipse_y = yd * 0.0689988270199406
        ellipse_z = zd
        
        ax.scatter(ellipse_x, ellipse_y, ellipse_z)
        ax.plot(line_x[i][0], line_y[i][0], line_z)
        
    ax.scatter(pt_coords[:,0], pt_coords[:,1], pt_coords[:,2], c = 'r', s = 50)
    return ax
    
def load_image(config_file, im_file):
    """Loads puncta image into an array. Copy and pasted this function from 
       Puncty"""
    parser = SafeConfigParser()
    parser.read(config_file)
    im = (im_file)
    nimg = parser.getint('section1', 'number_of_slices')
    img = Image.open(im)
    img.seek(0)
    nx, ny = img.size
    print "Allocating array"
    big_arr = np.zeros([nimg, nx, ny], dtype=np.uint8)
    for i in range(nimg):
        img.seek(i)
        big_arr[i] = np.array(img, dtype=np.uint8)
    return big_arr, parser
    
def extract_puncta(image ,threshold, x_offset, y_offset):
    """Filters and thresholds puncta image. Then finds location of puncta in 3D
       space with centre_of_mass. Also figures out the number of voxels which 
       make up the puncta #copyandpaste"""
       
    print "Extracting puncta"
    #filtered_arr = ndi.median_filter(image, 2)
    thrshd = image > threshold
    
    #opened = ndi.binary_opening(thrshd, iterations = 2)
    labels, num = label(thrshd)
    count, val = np.histogram(labels, bins=labels.max())
    val = val.astype(int)
    
    for lbl in xrange(len(count)):
        if count[lbl] > 10000:
            count[lbl] = 0 
    
    too_small_labels = val[count < 4]
    
    for lbl in too_small_labels:
        labels[labels == lbl] = 0
                    
    puncta_coords = center_of_mass(thrshd, labels=labels, 
                                   index=val[count > 4])
                                    
    for i in xrange(len(puncta_coords)):
        puncta_coords[i] = list(puncta_coords[i])
        puncta_coords[i][0] =  puncta_coords[i][0]
        puncta_coords[i][0] *= 0.2
    
    """ I think it's upside down? whhaaaat?!? Figured it out... it was 
    upside down. Highest z-level for puncta is pixel value 0 not 1024
    """                          
    for i in xrange(len(puncta_coords)):
        puncta_coords[i] = list(puncta_coords[i])
        puncta_coords[i][1] += y_offset 
        puncta_coords[i][1] = 1024 - puncta_coords[i][1]
        puncta_coords[i][1] *= 0.0689988270199406
        
    for i in xrange(len(puncta_coords)):
        puncta_coords[i] = list(puncta_coords[i])
        puncta_coords[i][2] += x_offset
        puncta_coords[i][2] *= 0.0689988270199406

    voxel = 0.0689988270199406 * 0.0689988270199406 * 0.2
    
    new_puncta_coords = []
    
    """rearranging the puncta coords"""
    for i in xrange(len(puncta_coords)):
        new_puncta_coords.append((puncta_coords[i][0], puncta_coords[i][2], puncta_coords[i][1] ))
    
    new_count = []
    for i in xrange(len(count)):
        if count[i] > 4:
            new_count.append(count[i] * voxel)
            
    #new_puncta_coords.append((35,15,35))
    return new_puncta_coords, thrshd, labels, new_count
    
def dist_point_from_curve(x_params, y_params, pt_coords):
    """First need to work out point on line which a punctum is closest too, using
       differentiation :) Then can work out distance between the two points
    
    dist_equation    = (a*z**2 + b*z + c - xp)**2 
                     + (e*z**2 + f*z + g - yp)**2 
                     + (z - zp)**2
    diff (=min_dist) = 2*(a*z**2 + b*z + c - xp)*(2*a*z + b) 
                     + 2*(e*z**2 + f*z + g - yp)*(2*e*z + f) 
                     + 2*(z - zp)
    where a, b, c  = x plane polynomial coeffs. e,f,g = y plane polynomial coeffs
          xp,yp,zp = x,y,z coordinates for each puncta
    """
    print "Differentiating axes"
    big_min_dist_list = []
    big_line_points_list = []
    num_cells = len(x_params)
    
    """for each of the central axes"""
    for i in xrange(num_cells):
        a,b,c = x_params[i][0], x_params[i][1], x_params[i][2]
        e,f,g = y_params[i][0], y_params[i][1], y_params[i][2]
        line_points_list = []
        min_dist_list = []
        aa = 4*a**2 + 4*e**2
        bb = 6*a*b + 6*e*f
        """for each of the puncta"""
        for ii in xrange(len(pt_coords)):
            xp,yp,zp = pt_coords[ii][0], pt_coords[ii][1], pt_coords[ii][2]
            cc = 2 + 2*b**2 + 4*a*c + 2*f**2 + 4*e*g - 4*e*yp - 4*a*xp
            dd = 2*c*b + 2*f*g - (2*b)*xp - (2*f)*yp - 2*zp
            """simplified = aa*z**3 + bb*z**2 + cc*z + dd"""
            coeffs = [aa,bb,cc,dd]
            zed = np.roots(coeffs)
            #print "z =", zed
            """Should produce 3 roots (3rd order polynomial), 
               one being real and two being complex (having imaginary components).
               Finding only real answers"""
            for iii in xrange(len(zed)):
                if zed[iii].imag == 0:
                    zl = zed[iii].real
            
            xl = a*zl**2 + b*zl + c
            yl = e*zl**2 + f*zl + g
                #print "y solution:", yl
                #print "x solution:", xl
            """3D version of pythagoras equation. Distance: point to point"""
            min_dist = np.sqrt((xl - xp)**2 + (yl - yp)**2 + (zl - zp)**2)
            #print "minimum distance:", min_dist
            line_points_list.append([xl,yl,zl])
            min_dist_list.append(min_dist)
        """all the distances for puncta to axes"""
        big_min_dist_list.append(min_dist_list)
        """all the nearest points on the axes"""
        big_line_points_list.append(line_points_list)
    return big_line_points_list, big_min_dist_list

def define_plane(line_points_list, x_params, y_params, pt_coords, ax):
    """Name says it all really. Defines the plane which encapsualtes both a
       punctum and its respective closest point to the central axes. This plane
       is normal to the tangent of the closest point on an axis. This then uses
       trigonometry to determine the angle of a punctum around the plane.
    """
    print "Defining planes"
    big_angle_list = []
    
    for i in xrange(len(x_params)):
        angle_list = []
        a,b = x_params[i][0], x_params[i][1]
        e,f = y_params[i][0], y_params[i][1]
        for ii in xrange(len(line_points_list[i])):
            xl, yl, zl = line_points_list[i][ii][0], line_points_list[i][ii][1], line_points_list[i][ii][2]
            x_grad = (2*a*zl + b)
            y_grad = (2*e*zl + f)
            #tangent gradient (x & y) = normal to plane
            n = np.array([x_grad, y_grad, 1])
    
            """ now can define a plane pl using pt_coords and n    
                find point pl_pt in pl where x = line_point_x,
                i.e. need to find pl_pt_y, pl_pt_z"""
            pl_pt = np.array([xl, yl + 5, zl - 5 * n[1]])
            #ax.scatter(pl_pt[0], pl_pt[1], pl_pt[2], c='c', s=50)
            """find angle between vectors (pt_coords - line_point) and (pl_pt - line_point) 
               use dot product for this since dot(a, b) / |a||b| = cos theta_ab """
            la = pt_coords[ii] - line_points_list[i][ii]
            lb = pl_pt - line_points_list[i][ii]
            angle_rads = np.arccos(np.dot(la, lb) / \
                        (np.linalg.norm(la) * np.linalg.norm(lb)))
            if pt_coords[ii][0] - xl < 0:
                angle_rads= np.pi - angle_rads + np.pi
            #print "Angle (deg):", np.degrees(angle_rads)
            plane_x, plane_y = np.meshgrid(range(-8, 8), range(-8, 8))
            plane_x += xl
            plane_y += yl
            plane_z = zl - n[0] * (plane_x) - n[1] * (plane_y)
            #ax.plot_surface(plane_x, plane_y, plane_z, alpha=0.4)
            #ax.scatter(xl,yl,zl, c='g', s=50)
            angle_list.append(angle_rads)
        big_angle_list.append(angle_list)
    return big_angle_list
        
def create_scatter_plot(angle_list, dist_list, puncta_coords):
    """creates a polar scatter plot using the calculated distances and angles
       Also uses the y-level of the puncta for colour coding because colours
       are nice"""
    print "Making pretty graphs"

    dist_list = np.array(dist_list).T
    angle_list = np.array(angle_list).T
    distance = []
    angle = []
    cell1 = []
    cell2 = []
    cell3 = []
    cell4 = []
    for i in xrange(len(dist_list)):
        distance.append(min(dist_list[i]))
        angle.append(angle_list[i][np.argmin(dist_list[i])])
        if (np.argmin(dist_list[i])) == 0:
            cell1.append(1)
        if (np.argmin(dist_list[i])) == 1:
            cell2.append(1)
        if (np.argmin(dist_list[i])) == 2:
            cell3.append(1)
        if (np.argmin(dist_list[i])) == 3:
            cell4.append(1)
    print "Number of Puncta..."               
    print "Cell 1:", len(cell1)
    print "Cell 2:", len(cell2)
    print "Cell 3:", len(cell3)
    print "Cell 4:", len(cell4)
    print "Total:", len(cell1 + cell2 + cell3 + cell4)
    zlevel = []
    for ipuncta in xrange(len(puncta_coords)):
        zlevel.append((puncta_coords[ipuncta][2]))
    fig = plt.figure(1)
    ax = fig.add_subplot(111, polar=True)
    ax.scatter(angle, distance, c = zlevel)
    plt.show()
    return angle
    
def calculate_density(angle_list):
    """determines the number of puncta within a 2pi/8th of a wedge"""
    density_list = []
    for x in range(8):
        for y in angle_list:
            if y > (np.pi/4 * x) and y < (np.pi/4 * x + np.pi/4):
                density_list.append(x)
                
    return density_list

def make_density_plot(density_list):
    
    """Uses the number of puncta within each 2pi/8 wedge to obtain a relative
    density comparing the density in wedge '0 - pi/4' with all other wedges.
    This relative density is used as the bar distance for the graph"""
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

    theta = []

    for x in range(4):
        eighth = 2 * x * np.pi/4 - np.pi/4 
        theta.append(eighth)
        
    radii = []
    
    for x in range(8):
        radii.append(density_list.count(x))
     
    new_radii = []
    
    new_radii.append(radii[0] + radii[7])
    new_radii.append(radii[1] + radii[2])
    new_radii.append(radii[3] + radii[4])
    new_radii.append(radii[5] + radii[6])
    
    total_puncta = sum(new_radii)
    
    perc_1 = float(new_radii[0])/ total_puncta *100
    perc_2 = float(new_radii[1])/ total_puncta *100
    perc_3 = float(new_radii[2])/ total_puncta *100
    perc_4 = float(new_radii[3])/ total_puncta *100
        
    relative_radii = []
    for x in new_radii:
        relative_radii.append(float(x)/float(new_radii[0]))
    
        
    width = np.pi/2
    ax.bar(theta, relative_radii, width=width, bottom=0.0)
    
    plt.show()
    return perc_1, perc_2, perc_3, perc_4
    
def puncta_size(new_count, angle_list):
    median = np.median(new_count)
    relative_puncta_size = []
    for i in xrange(len(new_count)):
        relative_puncta_size.append(new_count[i]/median)
    segment1 = []
    segment2 = []
    segment3 = []
    segment4 = []
    
    for i in xrange(len(angle_list)):
        if angle_list[i] > 0 and angle_list[i] < np.pi/4:
            segment1.append(relative_puncta_size[i])
        if angle_list[i] > np.pi * 1.75 and angle_list[i] < np.pi *2:
            segment1.append(relative_puncta_size[i])
        if angle_list[i] > np.pi/4 and angle_list[i] < np.pi/2:
            segment2.append(relative_puncta_size[i])
        if angle_list[i] > np.pi/2 and angle_list[i] < np.pi/4 *3:
            segment2.append(relative_puncta_size[i])
        if angle_list[i] > np.pi/4 *3 and angle_list[i] < np.pi:
            segment3.append(relative_puncta_size[i])
        if angle_list[i] > np.pi and angle_list[i] < np.pi + np.pi/4:
            segment3.append(relative_puncta_size[i])
        if angle_list[i] > np.pi + np.pi/4 and angle_list[i] < np.pi *1.5:
            segment4.append(relative_puncta_size[i])
        if angle_list[i] > np.pi *1.5 and angle_list[i] < np.pi *1.75:
            segment4.append(relative_puncta_size[i])
        
    right_size = sum(segment1)/len(segment1)
    neural_size = sum(segment2)/len(segment2)
    left_size = sum(segment3)/len(segment3)
    abneural_size = sum(segment4)/len(segment4)
    
    return neural_size, right_size, abneural_size, left_size

    
def misc_arrange(config_file, wdir, z_list, y_min, y_max, x_min, x_max, threshold):
    z_list = np.array(z_list) 
    z_list = (1024 - z_list)  * 0.0689988270199406
    parser = SafeConfigParser()     
    coord_list = load_data(parser, config_file)
    centre_list, fitted_list, xm_list, ym_list = fit_many_ellipses(coord_list,z_list)
    data = np.array(centre_list)
    fitted_list =  np.array(fitted_list)
    x_params_list, y_params_list = fit_curved_axis(data)
    line_z = np.linspace(z_list[0], z_list[-1] - 2, 1000)
    big_line_x = []
    big_line_y = []
    for i in xrange(len(x_params_list)):
        line_x = []
        line_y =[]
        x_params = x_params_list[i]
        y_params = y_params_list[i]
        line_x.append(curve_func(line_z, *x_params))
        line_y.append(curve_func(line_z, *y_params))
        big_line_x.append(line_x)
        big_line_y.append(line_y)
    parser.read(config_file)
    mtiff_file = parser.get('section1','puncta_file')
    im_file = os.path.join(wdir, mtiff_file)
    big_arr, config_parser = load_image(config_file, im_file)    
    smaller_arr = big_arr[:,y_min:y_max,x_min:x_max]    
    puncta_coords, thrshd, labels, new_count = extract_puncta(smaller_arr, threshold, 
                                                              x_min, y_min)
    puncta_coords = np.array(puncta_coords)
    get_ellipse_axes(fitted_list)
    ax = show_fitted_curve(big_line_x, big_line_y, line_z, data,  puncta_coords)
    line_points_list, min_dist_list = dist_point_from_curve(x_params_list, y_params_list, puncta_coords)

    angle_list = define_plane(line_points_list, x_params_list, y_params_list,  puncta_coords, ax)
    angle = create_scatter_plot(angle_list, min_dist_list, puncta_coords)
    density_list = calculate_density(angle)
    perc_1, perc_2, perc_3, perc_4 = make_density_plot(density_list)
    """
    print "Basal", perc_1
    print "Neural", perc_2
    print "Apical", perc_3
    print "Abneural", perc_4
    neural_size, right_size, abneural_size, left_size = puncta_size(new_count, angle)
    print "Basal size:", right_size
    print "Neural size:", neural_size
    print "Apical size:", left_size
    print "Abneural size:", abneural_size
    """
    
if __name__=='__main__':
    """Input section"""
    config_file = 'C:\Users\DELL\Desktop\curved_config\P6_1_ribbon.ini'
    wdir = os.path.normpath(r'C:/Users/Dell/Desktop/confocal images/glyofixx/p6')
    """y slice levels from Clicky 2.0""" 
    z_list = [260,300,350,400,450]
    """The array slicing dimensions used for making the smaller 3D array"""
    array_y_min = 325
    array_y_max = 514
    array_x_min = 65
    array_x_max = 540
    threshold = 40

    misc_arrange(config_file, wdir, z_list, array_y_min, array_y_max,
                 array_x_min, array_x_max, threshold)