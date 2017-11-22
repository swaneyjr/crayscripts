#!/usr/bin/env python

import random
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import brentq
from collections import namedtuple
from operator import itemgetter


Pixel = namedtuple('Pixel',['x','y','val'])

# experiment options
phone = "S4"
ROTATED = False
n_days = 30
res = "720p"
pix_depth = 15 # in microns

# phone settings
n_min = 5
n_max = 10000
l1thresh = 10
l2thresh = 8
average_pix = False

# convolving algorithm
n_dim = 7
k = 3
CONVOLUTION = np.zeros((n_dim,n_dim))
for i,j in np.ndindex(CONVOLUTION.shape):
    CONVOLUTION[i][j] = math.exp(-k*math.sqrt((i-n_dim//2)**2+(j-n_dim//2)**2))

# plot options
plotScatter = True
plotLenStats = True
plotHoughStats = True
plotRhoVals = True

# efficiency as a function of angle of incidence and energy
def random_val(energy, isotropic):
    if not isotropic:
        return max(0,random.expovariate(.01))
    else:
        return 0

# generates an energy distribution
def random_energy():

    return 0



muon_flux = 1/60
isotropic_flux = 0
    
phone_dict = {"S4": [(4128,3096), .000112], "S3": [(3264,2448), .000147]}
res_dict = {"1080p": (1920,1080), "720p": (1280,720), "480p": (640,480)}

resx,resy = res_dict[res]
max_res, phys_pix = phone_dict[phone]
max_resx, max_resy = max_res

phys_per_vis = min(max_resx//resx, max_resy//resy)
vis_pix = phys_pix*phys_per_vis
pix_depth /= 10000 # microns -> cm

cmos_x, cmos_y = vis_pix*resx, vis_pix*resy

if ROTATED:
    muon_flux /= 2
t_tot = n_days*24*3600


# calculates the angle of a track via hough transform
def hough_theta(pixels):

    weight = False

    if weight:
        intersections = []
    

        def theta_weight(dx, dy):
            return dx**2 + dy**2
        
        for i,p1 in enumerate(pixels):
            for p2 in pixels[i+1:]:
                dx = p1.x-p2.x
                dy = p1.y-p2.y
                weight = theta_weight(dx, dy)         
                if dy == 0:
                    intersections.append((math.pi/2, weight))
                else:
                    intersections.append((math.atan(-dx/dy), weight))

        # take weighted average
        theta_weighted_sum = 0
        for th, weight in intersections:
            theta_weighted_sum += cmath.rect(weight, 2*th)
        theta = cmath.phase(theta_weighted_sum)/2
        rho = 0
        for p in pixels:
            rho += p.x*math.cos(theta) + p.y*math.sin(theta)
        rho /= len(pixels)

    else:
        # orthogonal least squares fit
        sum_x = sum_y = sum_x_sq = sum_y_sq = sum_xy = 0
        n = len(pixels)
        for p in pixels:
            sum_x += p.x
            sum_y += p.y
            sum_x_sq += p.x**2
            sum_y_sq += p.y**2
            sum_xy += p.x*p.y

        num = n*sum_xy-sum_x*sum_y
        den = n*sum_x_sq-n*sum_y_sq-sum_x**2+sum_y**2
        if den == 0:
            theta = math.copysign(math.pi/4, -num)
        else:
            theta = 1/2 * math.atan2(-2*num, -den)
        rho = (sum_x*math.cos(theta)+sum_y*math.sin(theta))/n

    # calculate other stats based on theta fit

    rho_vals = [p.x*math.cos(theta)+p.y*math.sin(theta) for p in pixels]
    pix_vals = [p.val for p in pixels]
    sigma_vals = [p.x*math.sin(theta)-p.y*math.cos(theta) for p in pixels]

    drho_val_pairs = [(abs(r-rho),v) for r,v in zip(rho_vals,pix_vals)]
    
    if len(pixels) > 2:
        rho_std = np.std(rho_vals, ddof = 2)
    else:
        rho_std = None
    sigma_vals, rho_vals = map(np.array, zip(*sorted(list(zip(sigma_vals,rho_vals)), key=itemgetter(0))))
    sigma_mid = np.add(sigma_vals[1:],sigma_vals[:-1])/2
    rho_diff = np.diff(rho_vals)
    if len(pixels) < 3:
        resid_slope = None
    else:
        resid_slope = np.polyfit(sigma_mid, rho_diff, deg=1)[1]
    
    return theta, drho_val_pairs, rho_std, resid_slope


# given a flux and source distribution, appends tracks onto the CMOS
def generate_tracks(flux, isotropic):

    class Point(object):
        def __init__(self, x, y):
            self.x = x
            self.y = y

        def add(self, p):
            return Point(self.x + p.x, self.y + p.y)

        def pixel(self):
            return Point(self.x/phys_pix, self.y/phys_pix)

        def pixel_xy(self):
            return self.x/phys_pix, self.y/phys_pix

    class LineSegment(object):
        def __init__(self, start, end):
            self.start = start
            self.end = end

        def y(self, x):
            slope = (self.end.y-self.start.y)/(self.end.x-self.start.x)
            return self.start.y + slope*(x-self.start.x)

    # calculates vector of travel through CMOS
    def find_dist(theta, phi):
        if ROTATED:
            return Point(pix_depth/math.tan(phi), pix_depth/math.tan(theta)/math.sin(phi))
        else:
            length = pix_depth * math.tan(theta)
            return Point(length*math.cos(phi), length*math.sin(phi))

    # generates a cos^2 distribution of theta values
    def random_theta(isotropic):

        if isotropic:
            return random.uniform(0, math.pi)
        else:
            rand = random.random()
            cdf_minus_rand = lambda theta: 2.0/math.pi*(theta + .5 * \
                math.sin(2.0 * theta)) - rand
            return brentq(cdf_minus_rand, 0., math.pi/2.0)

    
    tracks = []
    hough_stats = []
    angle_diffs = []
    drho_val_pairs = []
    
    t = 0
    if flux > 0:
        percent = 10
        while t < t_tot:
            t += random.expovariate(flux * cmos_x * cmos_y)
            if (t*100)/t_tot >= percent:
                print(' ',percent,'%', sep='')
                percent += 10
            
            # generate trajectory of muon
            theta = random_theta(isotropic)
            phi = random.uniform(-math.pi, math.pi)
            energy = random_energy()
            
            # compute muon path through CMOS
            dist = find_dist(theta, phi)
            start = Point(random.uniform(0,cmos_x), random.uniform(0,cmos_y))
            end = start.add(dist)
            start_pix = start.pixel()
            end_pix = end.pixel()
            path = LineSegment(start_pix, end_pix)
            
            # determine affected pixels
            pixels = []
            start_x, start_y = start.pixel_xy()
            end_x, end_y = end.pixel_xy()
            min_x = min(start_x, end_x)
            max_x = max(start_x, end_x)
            x_vals = [min_x] + list(range(int(min_x)+1, int(max_x)+1)) + [max_x]
            x_vals = [x for x in x_vals if x >= 0 and x < max_resx]
                    
            for i in range(len(x_vals)-1):
                y1 = path.y(x_vals[i])
                y2 = path.y(x_vals[i+1])
                y_vals = list(range(int(min(y1, y2)), int(max(y1, y2))+1))
                y_vals = [y for y in y_vals if y >= 0 and y < max_resy]
                for y in y_vals:
                    pixels.append(Pixel(int(x_vals[i]),y, random_val(energy, isotropic)))
            pixels = process(pixels)
            if len(pixels) >= n_min and len(pixels) <= n_max and pixels != []:
                tracks.append(pixels)
                if len(pixels) > 1:
                    h_theta, drho_val_track, rho_std, resid_slope = hough_theta(pixels)
                    hough_stats.append((h_theta, rho_std, resid_slope))
                    drho_val_pairs += drho_val_track
                    theta_diff = h_theta - phi - math.pi/2
                    while theta_diff > math.pi/2:
                        theta_diff -= math.pi
                    while theta_diff < -math.pi/2:
                        theta_diff += math.pi
                    angle_diffs.append(theta_diff)

    return tracks, hough_stats, angle_diffs, drho_val_pairs

# convolutes and samples pix vals
def process(raw_pix):

    # simulate energy transfer
    if raw_pix == []: return []
    pix_x, pix_y, pix_val = zip(*raw_pix)
    
    min_x = min(pix_x)
    x_range = max(pix_x)-min_x+n_dim
    min_y = min(pix_y)
    y_range = max(pix_y)-min_y+n_dim

    vals_array = np.zeros((x_range, y_range))
    
    for p in raw_pix:
        for i,j in np.ndindex(CONVOLUTION.shape):
            vals_array[p.x-min_x+i][p.y-min_y+j] += p.val * CONVOLUTION[i][j]

    # apply sampling algorithm
    vals_array = sample_vals(vals_array, min_x-n_dim//2, min_y-n_dim//2)
    min_samp_x = 2*(min_x - n_dim//2)//(2*phys_per_vis)
    min_samp_y = 2*(min_y - n_dim//2)//(2*phys_per_vis)

    # demosaic and convert to luminance values
    vals_array = luminance_vals(vals_array, min_samp_x, min_samp_y)

    if np.amax(vals_array) < l1thresh: return []

    
    new_pix = []

    for i,j in np.ndindex(vals_array.shape):
        if vals_array[i][j] >= l2thresh and i >= -min_samp_x and j >= -min_samp_y \
           and i < resx - min_samp_x and j < resy - min_samp_y:
            new_pix.append(Pixel(i+min_samp_x, j+min_samp_y, min(vals_array[i][j],256)))
 
    return new_pix

# from array of physical pixel values, returns visual pixel values
# according to Bayer filter pattern
def sample_vals(vals_array, min_x_phys, min_y_phys):
    
    vals_size = vals_array.shape
    
    min_x_vis = 2*min_x_phys//(2*phys_per_vis)
    min_y_vis = 2*min_y_phys//(2*phys_per_vis)
    max_x_vis = 2*(min_x_phys+vals_size[0]-1)//(2*phys_per_vis)+1
    max_y_vis = 2*(min_y_phys+vals_size[1]-1)//(2*phys_per_vis)+1
    
    bayer_vals_array = np.zeros((phys_per_vis**2, max_x_vis-min_x_vis+1, max_y_vis-min_y_vis+1))

    # assign to spots in Bayer pattern
    for iphys,jphys in np.ndindex(vals_size):
        ivis = 2*(iphys + min_x_phys)//(2*phys_per_vis) - min_x_vis + (iphys % 2)
        jvis = 2*(jphys + min_y_phys)//(2*phys_per_vis) - min_y_vis + (jphys % 2)
        position = (iphys%(2*phys_per_vis))//2 + phys_per_vis*((jphys%(2*phys_per_vis))//2)
        bayer_vals_array[position][ivis][jvis] = vals_array[iphys][jphys]

    if average_pix:
        return np.mean(bayer_vals_array, axis=0)
    else:
        np.random.shuffle(bayer_vals_array)
        return bayer_vals_array[0][:][:]

# simulate demosaicing and processing RGB -> BW
def luminance_vals(vals_array, min_x, min_y):

    vals_size = vals_array.shape

    def avg_partners(vals_array,i,j,partners):
        val = 0
        valx,valy = vals_array.shape
        for pi,pj in partners:
            if i >= -pi and j >= -pj and i < valx-pi and j < valy-pj:
                val += vals_array[i+pi,j+pj]
        return val/len(partners)
    
    adj = [(1,0),(0,1),(-1,0),(0,-1)]
    diag = [(1,1),(1,-1),(-1,-1),(-1,1)]
    lr = [(1,0),(-1,0)]
    ud = [(0,1),(0,-1)]
    
    # obtain rgb values
    rgb_array = np.zeros((vals_size[0], 3, vals_size[1]))
    for i,j in np.ndindex(vals_size):
        if (i+min_x)%2 == 0 and (j+min_x)%2 == 0: #red
            rgb_array[i,0,j] = vals_array[i,j]
            rgb_array[i,1,j] = avg_partners(vals_array, i, j, adj)
            rgb_array[i,2,j] = avg_partners(vals_array, i, j, diag)
        elif (i+min_x)%2 == 1 and (j+min_x)%2 == 1: #blue
            rgb_array[i,0,j] = avg_partners(vals_array, i, j, diag)
            rgb_array[i,1,j] = avg_partners(vals_array, i, j, adj)
            rgb_array[i,2,j] = vals_array[i,j]
        elif (i+min_x)%2 == 1 and (j+min_x)%2 == 0: #green
            rgb_array[i,0,j] = avg_partners(vals_array, i, j, lr)
            rgb_array[i,1,j] = vals_array[i,j]
            rgb_array[i,2,j] = avg_partners(vals_array, i, j, ud)
        if (i+min_x)%2 == 0 and (j+min_x)%2 == 1: #green
            rgb_array[i,0,j] = avg_partners(vals_array, i, j, ud)
            rgb_array[i,1,j] = vals_array[i,j]
            rgb_array[i,2,j] = avg_partners(vals_array, i, j, lr)
            

    # find sRGB Y value
    def find_luminance(rgb_array):
        rgb_transform = np.array([.2126, .7152, .0722])
        return np.dot(rgb_transform, rgb_array)

    luminance_array = find_luminance(rgb_array)
        
            
    return luminance_array

    
print('Generating muon tracks...')
muon_tracks, muon_hough, muon_angle_diffs, muon_rho_val = generate_tracks(muon_flux, isotropic=False)
print('Generating isotropic background...')
iso_tracks, iso_hough, iso_angle_diffs, iso_rho_val = generate_tracks(isotropic_flux, isotropic = True)

tracks = muon_tracks + iso_tracks

hough_stats = muon_hough + iso_hough
theta, rho_std, resid_slope = zip(*hough_stats)

rho_std = [r for r in rho_std if r != None]
resid_slope = [r for r in resid_slope if r!= None]

angle_diffs = muon_angle_diffs + iso_angle_diffs
rho_val_pairs = muon_rho_val + iso_rho_val

pix_n = [len(tr) for tr in tracks]

discr_len = [] # distance of tracks in discrete metric
tr_len = [] # Euclidean metric


#
# make plots
#

print('Plotting stats...')
plotSum = int(plotScatter)
if plotScatter:
    plt.figure(plotSum)
    #xy scatter
    Colors = plt.get_cmap('gist_rainbow')
    for i in range(len(tracks)):
        pix_x, pix_y, pix_val = zip(*tracks[i])
        discr_len.append(max(pix_x)-min(pix_x)+max(pix_y)-min(pix_y)+1)
        tr_len.append(math.sqrt((max(pix_x)-min(pix_x))**2+(max(pix_y)-min(pix_y))**2))

        colz = Colors(i/len(tracks))
        
        plt.scatter(pix_x, pix_y, s=1, marker = ',',c=colz, edgecolors='face')
        plt.title('Tracks on CMOS')
        plt.xlabel('pix_x')
        plt.ylabel('pix_y')

    max_pix = max(discr_len)
    pix_n_cutoff = [min(n, 35) for n in pix_n]
    discr_len_cutoff = [min(d, 50) for d in discr_len]
    tr_len_cutoff = [min(l, 50) for l in tr_len]

if plotLenStats:
    plotSum += plotLenStats
    plt.figure(plotSum)
    plt.subplot(221)
    plt.hist(pix_n_cutoff, bins=min(max_pix,35)-n_min+1, range=(n_min, min(max_pix, 50)))
    plt.xlabel("pix_n")
    plt.ylabel("Frequency")
    plt.xscale("log")
    plt.yscale("log")

    plt.subplot(222)
    plt.hist(discr_len_cutoff, bins=min(max_pix,50)-n_min+1, range=(n_min, min(max_pix,50)))
    plt.xlabel("Discrete Length (Pixels)")
    plt.ylabel("Frequency")
    plt.xscale("log")
    plt.yscale("log")

    plt.subplot(223)
    plt.hist(tr_len_cutoff, bins=min(max_pix,50)-n_min+1, range=(n_min, min(max_pix,50)))
    plt.xlabel("Track Length (Pixels)")
    plt.ylabel("Frequency")
    plt.xscale("log")
    plt.yscale("log")

    if max_pix > 2:
        plt.subplot(224)
        efficiency = [(pix_n[i]-2)/(discr_len[i]-2) for i in range(len(tracks)) if discr_len[i]>2]
        plt.hist(efficiency, bins = 50, weights = [d-2 for d in discr_len if d>2])
        plt.xlabel("Efficiency")
        plt.ylabel("Frequency")

if plotHoughStats:
    plotSum += plotHoughStats
    plt.figure(plotSum)
    plt.subplot(221)
    plt.hist(theta, bins = 50, range=(-math.pi/2, math.pi/2))
    plt.xlabel(r'$\theta$')
    plt.ylabel("Frequency")

    if not ROTATED:
        plt.subplot(222)
        plt.hist(angle_diffs, bins = 50)
        plt.xlabel(r'$\theta - \phi$')
        plt.ylabel("Frequency")

    plt.subplot(223)
    plt.hist(rho_std, bins = 50)
    plt.xlabel(r'$\sigma_\rho$')
    plt.ylabel('Frequency')

    plt.subplot(224)
    plt.hist(resid_slope, bins = 50)
    plt.xlabel(r'Residual slope')
    plt.ylabel('Frequency')

if plotRhoVals:
    plotSum += plotRhoVals
    plt.figure(plotSum)
    drho_vals, val_vals = zip(*rho_val_pairs)
    max_val = max(val_vals)
    plt.hist2d(drho_vals, val_vals, bins = [30,(max_val-l2thresh)//2+1], \
               range=[(0,max(drho_vals)),(l2thresh,max_val)], norm=LogNorm())
    plt.xlabel(r'$\rho-\bar{\rho}$')
    plt.ylabel('pix_val')
    plt.colorbar()


plt.show()






