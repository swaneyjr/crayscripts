#!/usr/bin/env python

import random
import math
import scipy.optimize as optimization
import numpy as np
import pylab as py
import matplotlib.pyplot as plt

t_max = 1.
n_frames = 15
n_bins = 150
minus_bg = True
eff_fit = False

frame_rate = 15.
bins_per_frame = 20

if n_frames != None:
    n_bins = int(bins_per_frame*(n_frames+0.5))
    t_max = (n_frames + 0.5)/frame_rate

led = []
muon = []

t_tot = 14*24*3600
paddle_side_x = math.sqrt(10)
paddle_side_y = math.sqrt(10)
A_paddle = paddle_side_x * paddle_side_y
cmos_side_x = .454
cmos_side_y = .352
A_cmos = cmos_side_x * cmos_side_y
cmos_gap = 2.5
telescope_gap = 5
muon_flux = 1/60
cmos_noise = 8000/24/3600
paddle_noise = .0001
coinc_efficiency = 0.9
expected_delay = 0.0 #offset for data analysis
real_delay = 0.21
led_delay_sigma = 0.06
p_shutter = 0.0 #LED-LED shutter effect
p_shutter_muon = 0.0 #LED-muon shutter effect
processing_rate = 120.

def frame_sort(diffs, n_frames):
    frames = np.zeros(2*n_frames+1)
    for i in range(len(diffs)):
        frames[i//bins_per_frame] += diffs[i]
    return frames

# generates a cos^2 distribution of theta values
def random_theta():
    
    rand = random.random()
    cdf_minus_rand = lambda theta: 2.0/math.pi*(theta + .5 * \
        math.sin(2.0 * theta)) - rand
    return optimization.brentq(cdf_minus_rand, 0., math.pi/2.0)

# find the value of f = P(tel | CMOS)
f = 0
trials = 50000
for i in range(trials):
    theta = random_theta()
    phi = random.uniform(0, 2.0 * math.pi)
    x = random.uniform(-cmos_side_x/2.0, cmos_side_x/2.0)
    y = random.uniform(-cmos_side_y/2.0, cmos_side_y/2.0)
    x_first = x - cmos_gap * math.tan(theta) * math.cos(phi)
    y_first = y - cmos_gap * math.tan(theta) * math.sin(phi)
    x_third = x_first + telescope_gap * math.tan(theta) * math.cos(phi)
    y_third = y_first + telescope_gap * math.tan(theta) * math.sin(phi)
    if max(abs(x_first), abs(x_third)) < paddle_side_x/2.0 \
        and max(abs(y_first), abs(y_third)) < paddle_side_y/2.0:
        f += 1.0
f /= float(trials)

extra_muon_rate = muon_flux * A_cmos * (1-f)

# returns number of coincidences with an associated muon-like hit
def find_matches(muon, led, t_coinc, led_delay):
    muon_array = np.array(sorted(muon))
    led_array = np.array(sorted(led))
    matches = 0
    m = 0
    l = 0
    while m < len(muon_array) and l < len(led_array):
        gap = led_array[l] - muon_array[m] - led_delay 
        if abs(gap) < abs(t_coinc) and np.sign(gap) == np.sign(t_coinc):
            l += 1
            m += 1
            matches += 1
        elif gap < 0:
            l += 1
        else:
            m += 1
    return matches
        

# runs a simulation of the muon telescope
def run_telescope(efficiency):

    print("Running simulation...")

    # append muon hits through the scintillator and CMOS chip
    #f1 = f2 = g1 = g2 = h1 = h2 = 0
    t = 0
    last_led = -real_delay
    while t < t_tot:
        delta_t = random.expovariate(muon_flux * A_paddle)
        theta = random_theta()
        phi = random.uniform(0, 2.0 * math.pi)
        #cmos = False
        if t > 0:
            x = random.uniform(-paddle_side_x/2.0, paddle_side_x/2.0)
            y = random.uniform(-paddle_side_y/2.0, paddle_side_y/2.0)
            #g1 += 1
            x += cmos_gap * math.tan(theta) * math.cos(phi)
            y += cmos_gap * math.tan(theta) * math.sin(phi)  
            if abs(x) < cmos_side_x/2.0 and abs(y) < cmos_side_y/2.0 \
               and random.random() < efficiency:
                muon.append(t)
                #cmos = True
                #f1 += 1
            x += (telescope_gap-cmos_gap) * math.tan(theta) * math.cos(phi)
            y += (telescope_gap-cmos_gap) * math.tan(theta) * math.sin(phi) 
            if abs(x) < paddle_side_x/2.0 and abs(y) < paddle_side_y/2.0 \
               and random.random() < coinc_efficiency \
               and t - last_led > real_delay: #circuit effects
                last_led = t + random.gauss(real_delay, led_delay_sigma)
                led.append(last_led)
                #g2 += 1
                #h1 += 1
                #if cmos:
                    #f2 += 1
                    #h2 += 1
        t += delta_t
    #f = f2/float(f1)
    #g = g2/float(g1)
    #h = h2/float(h1)
    #print("f = ", f)
    #print("g = ", g)
    #print("h = ", h)

    # add muons that missed first paddle
    t = 0
    while t < t_tot:
        if t > 0:
            muon.append(t)
        t += random.expovariate(extra_muon_rate)

    # append additional noise on the CMOS chip
    t = 0
    while t < t_tot:
        if t > 0:
            muon.append(t)
        t += random.expovariate(cmos_noise)
    # and the scintillators
    t = 0
    while t < t_tot:
        if t > 0:
            led.append(t)
        t += random.expovariate(paddle_noise)

    return


# find efficiency of simulation

efficiency = float(input("Efficiency = "))
run_telescope(efficiency)

# get lists of frame numbers
print("Processing times...")
muon_frames = [(int)(m*frame_rate) for m in muon]
led_frames = [(int)(l*frame_rate) for l in led]

# rolling shutter effect
led_frames_plus_shutter = []
for l in led_frames:
    r = 0.5 - random.random()
    if abs(r) < p_shutter_muon/2.0:
        muon_frames.append(l + np.sign(r))
    elif abs(r) < p_shutter/2.0:
        led_frames_plus_shutter.append(l + np.sign(r))
    led_frames_plus_shutter.append(l)

# remove copies
muon_frames = list(set(muon_frames))
led_frames = list(set(led_frames_plus_shutter))

print("Muon count = ", len(muon_frames))
print("LED count = ", len(led_frames))

# generate processed timestamps
muon = [m/frame_rate + random.expovariate(processing_rate) for m in muon_frames]
led = [l/frame_rate + random.expovariate(processing_rate) for l in led_frames]



print("Counting matches...")

n_pts = 2*n_bins + 1 
matches_array = np.zeros(n_pts)
t_step = 2 * t_max / float(n_pts - 1)

for i in range(n_pts):
    matches_array[i] = find_matches(muon, led, -t_max + i*t_step, expected_delay)
    if i % 50 == 0:
        print(" {0:.1f}% ({1}/{2})".format(100.*i/n_pts ,i, n_pts))
diffs = abs(np.diff(matches_array))

# sort into relevant bins
if n_frames != None:
    print("Sorting into frames...")
    data_x = np.arange(-n_frames, n_frames+1)
    data_y = frame_sort(diffs, n_frames)
else:
    data_x = np.linspace(-t_max + 0.5*t_step, t_max - 0.5*t_step, n_pts - 1) + expected_delay
    data_y = diffs
            
if minus_bg:
    data_x = -data_x[0:(len(data_x)+1)//2]
    for j in range((len(data_y)+1)//2):
        data_y[j] = abs(data_y[len(data_y)-j-1]) - abs(data_y[j])
    data_y = data_y[0:(len(data_y)+1)//2]

plt.hist(data_x, weights = data_y, bins = data_x.size)
plt.xlabel(r'$t_{LED}-t_{muon}$')
plt.ylabel("Matches")
if eff_fit:
    print("Calculating efficiency...")

    guess = np.array([0.3, expected_delay, .05])
    sigma = np.sqrt(diffs)

    # model
    def n(t, eff, mu_t, sigma_t):
        n_real = eff * coinc_efficiency * muon_flux * t_tot * A_cmos * f
        return t_step*n_real / math.sqrt(2. * math.pi * sigma_t**2) * \
            np.exp(-(t-mu_t)**2/(2.*sigma_t**2)) + \
             (len(led) - n_real) * (len(muon) - n_real)*t_step/t_tot * \
             np.exp((-len(muon)-len(led)+2*n_real)*t_step/t_tot)


    x = np.linspace(-t_max + expected_delay, t_max + expected_delay, 3*n_pts)
    y = n(x, efficiency, expected_delay + led_delay_mu, led_delay_sigma)
    py.plot(x,y, '-r', label = 'Predicted curve')
    popt, pcov = optimization.curve_fit(n, t_coinc, diffs, guess, sigma)
    zopt = n(x, popt[0], popt[1], popt[2])
    py.plot(x,zopt, '-b', label = 'Best fit')
    print("Efficiency ~ ", popt[0])
    
plt.show()





        

                      

