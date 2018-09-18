#!/usr/bin/env python3

import numpy as np
import random
from scipy.optimize import brentq

# everything in mm
cmos_side_x = 5312 * 1.12e-3
cmos_side_y = 2988 * 1.12e-3
paddle_side_x = 13
paddle_side_y = 15
cmos_gap = 7.5 #distance from first paddle to CMOS
tel_gap = 15 #distance from first paddle to last paddle
tel_eff = 0.91
muon_flux = 1/(60*100)

n = 100000


# generates a cos^2 distribution of theta values in [0,2*pi)
def random_theta(n):
    rand = np.random.uniform(0, np.pi / 2, n)
    cdf = lambda theta, r: theta + np.sin(2.0 * theta) / 2 - r
    thetas = [brentq(cdf, 0., np.pi/2.0, args=rand[i]) for i in range(n)]
    return np.array(thetas)

# geometrical acceptance of telescope
tan_theta = np.tan(random_theta(n))
phi = np.random.uniform(0, 2.0 * np.pi, n)
cos_phi = np.cos(phi)
sin_phi = np.sin(phi)

x0 = np.random.uniform(-paddle_side_x/2.0, paddle_side_x/2.0, n)
y0 = np.random.uniform(-paddle_side_y/2.0, paddle_side_y/2.0, n)
xcmos = x0 + cmos_gap * tan_theta * cos_phi
ycmos = y0 + cmos_gap * tan_theta * sin_phi
cmos = (np.abs(xcmos) < cmos_side_x/2.0) & (np.abs(ycmos) < cmos_side_y/2.0)

xtel = x0 + tel_gap * tan_theta * cos_phi
ytel = y0 + tel_gap * tan_theta * sin_phi
tel = (np.abs(xtel) < paddle_side_x/2.0) & (np.abs(ytel) < paddle_side_y/2.0) \
    & (np.random.rand(n) < tel_eff)

g = tel_eff * tel.sum() / n
print("P(tel  | 1st) = ", g)
print("P(cmos | tel) = ", (tel & cmos).sum() / tel.sum())

# fraction of muons through CMOS also passing through paddles
tan_theta = np.tan(random_theta(n))
phi = np.random.uniform(0, 2.0 * np.pi, n)
cos_phi = np.cos(phi)
sin_phi = np.sin(phi)

xcmos = np.random.uniform(-cmos_side_x/2.0, cmos_side_x/2.0, n)
ycmos = np.random.uniform(-cmos_side_y/2.0, cmos_side_y/2.0, n)
x0 = xcmos - cmos_gap * tan_theta * cos_phi
y0 = ycmos - cmos_gap * tan_theta * sin_phi
xtel = x0 + tel_gap * tan_theta * cos_phi
ytel = y0 + tel_gap * tan_theta * sin_phi
p0 = (np.abs(x0) < paddle_side_x/2.0) & (np.abs(x0) < paddle_side_y/2.0) \
        & (np.random.rand(n) < tel_eff)
ptel = (np.abs(xtel) < paddle_side_x/2.0) & (np.abs(ytel) < paddle_side_y/2.0) \
        & (np.random.rand(n) < tel_eff)

f = (p0 & ptel).sum() / n

print("P(tel  | cmos) = ", (p0 & ptel).sum() / n)
print()
print("Muon rate:      %.3E Hz" % (muon_flux * cmos_side_x * cmos_side_y))
print("Telescope rate: %.3E Hz" % (tel_eff * g * muon_flux * paddle_side_x * paddle_side_y))
print("Match rate:     %.3E Hz" % (f * muon_flux * cmos_side_x * cmos_side_y))

