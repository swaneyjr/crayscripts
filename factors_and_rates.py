import random
import math
from scipy.optimize import brentq

cmos_side_x = .454
cmos_side_y = .352
paddle_side_x = math.sqrt(10.)
paddle_side_y = math.sqrt(10.)
cmos_gap = 1. #distance from first paddle to CMOS
telescope_gap = 2. #distance from first paddle to last paddle
muon_flux = 1./60.

trials = 1000000


# generates a cos^2 distribution of theta values in [0,2*pi)
def random_theta():
    
    rand = random.random()
    cdf_minus_rand = lambda theta: 2.0/math.pi*(theta + .5 * \
        math.sin(2.0 * theta)) - rand
    return brentq(cdf_minus_rand, 0., math.pi/2.0)

# geometrical acceptance of telescope
g = 0
h1 = h2 = 0

for i in range(trials):
    cmos = False
    theta = random_theta()
    phi = random.uniform(0, 2.0 * math.pi)
    x1 = random.uniform(-paddle_side_x/2.0, paddle_side_x/2.0)
    y1 = random.uniform(-paddle_side_y/2.0, paddle_side_y/2.0)
    x2 = x1 + cmos_gap * math.tan(theta) * math.cos(phi)
    y2 = y1 + cmos_gap * math.tan(theta) * math.sin(phi)
    if abs(x2) < cmos_side_x/2.0 and abs(y2) < cmos_side_y/2.0:
        cmos = True
    x3 = x1 + telescope_gap * math.tan(theta) * math.cos(phi)
    y3 = y1 + telescope_gap * math.tan(theta) * math.sin(phi)
    if abs(x3) < paddle_side_x/2.0 and abs(y3) < paddle_side_y/2.0:
        g += 1.0/float(trials)
        h1 += 1
        if cmos:
            h2 += 1
print("g = ", g)
print("h = ", h2/float(h1))

# fraction of muons through CMOS also passing through paddles
f = 0
for i in range(trials):
    theta = random_theta()
    phi = random.uniform(0, 2.0 * math.pi)
    x2 = random.uniform(-cmos_side_x/2.0, cmos_side_x/2.0)
    y2 = random.uniform(-cmos_side_y/2.0, cmos_side_y/2.0)
    x1 = x2 - cmos_gap * math.tan(theta) * math.cos(phi)
    y1 = y2 - cmos_gap * math.tan(theta) * math.sin(phi)
    x3 = x1 + telescope_gap * math.tan(theta) * math.cos(phi)
    y3 = y1 + telescope_gap * math.tan(theta) * math.sin(phi)
    if max(abs(x1), abs(x3)) < paddle_side_x/2.0 \
        and max(abs(y1), abs(y3)) < paddle_side_y/2.0:
        f += 1.0
f /= float(trials)
print("f = ", f)

print("Expected muon candidate rate: ", muon_flux * cmos_side_x \
      * cmos_side_y, "candidates/s")
print("Expected LED rate: ", g * muon_flux * paddle_side_x \
      * paddle_side_y, " LED/s")
print("Expected rate of muon-LED matches: ", f * muon_flux \
      * cmos_side_x * cmos_side_y, " matches/s")

print("g/f^2 = ", g/f**2)

