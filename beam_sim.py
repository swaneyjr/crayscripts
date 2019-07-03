#!/usr/bin/env python3

import os
import time
import uuid

import numpy as np
from scipy.sparse import coo_matrix, save_npz
from scipy.stats import multivariate_normal

BEAM_TIME = 4.2 # seconds
BEAM_PERIOD = 60 # seconds
PIX_SIZE = .00112 # mm

RES_X = 4656
RES_Y = 3492

def map_coords(xp, yp, phi_p, theta_p, xs, ys, zs, phi_s, theta_s, psi_s):
    """
    Parameters:

    xp, yp:                 initial position of the particle in the beam frame
    phi_p, theta_p:         trajectory of the particle with respect to the 
                            beam axis
    xs, ys:                 Coordinates of the center of the sensor in the
                            lab frame, relative to its intersection with
                            the beam axis
    zs:                     distance from particle's starting coordinates to 
                            the plane of the beam, as measured along the beam 
                            axis
    phi_s, theta_s, psi_s:  Euler angles of the sensor plane
   

    Returns:
    
    x_hit, y_hit:           coordinates of the particle hit relative to the
                            center of the sensor
    """

    # convert particle positions to polar coords
    rp = np.sqrt(xp**2 + yp**2)
    xi_p = np.arctan2(yp, xp)

    # convert sensor positions to polar coords
    rs = np.sqrt(xs**2 + ys**2).reshape(-1,1)
    xi_s = np.arctan2(ys, xs).reshape(-1,1)

    # calculate trig functions
    cos_xi_p = np.cos(xi_p - phi_s.reshape(-1,1))
    sin_xi_p = np.sin(xi_p - phi_s.reshape(-1,1))
    
    cos_xi_s = np.cos(xi_s - phi_s.reshape(-1,1))
    sin_xi_s = np.sin(xi_s - phi_s.reshape(-1,1))

    cos_phi = np.cos(phi_p - phi_s.reshape(-1,1))
    sin_phi = np.sin(phi_p - phi_s.reshape(-1,1))

    cos_psi = np.cos(psi_s).reshape(-1,1)
    sin_psi = np.sin(psi_s).reshape(-1,1)

    sec_theta = 1 / np.cos(theta_s).reshape(-1,1)
    tan_theta = np.tan(theta_s).reshape(-1,1)
    
    denom = 1 / np.tan(theta_p) - tan_theta * sin_phi

    zt = zs.reshape(-1,1) - rp*cos_xi_p*tan_theta + rs*cos_xi_s*tan_theta

    # result
    x_hit = rp * (cos_xi_p*cos_psi*sec_theta + sin_xi_p*sin_psi) \
          - rs * (cos_xi_s*cos_psi*sec_theta + sin_xi_s*sin_psi) \
          + zt * (sin_phi *sin_psi*sec_theta + cos_phi *cos_psi) / denom
    y_hit = rp * (sin_xi_p *cos_psi - cos_xi_p *sec_theta*sin_psi) \
          - rs * (sin_xi_s *cos_psi - cos_xi_s *sec_theta*sin_psi) \
          + zt * (sin_phi*cos_psi*sec_theta - cos_phi*sin_psi) / denom

    return x_hit, y_hit


def visualize(spot_size, n_particles, x_s, y_s, phi_s, theta_s, psi_s):

    global RES_X, RES_Y, PIX_SIZE

    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    fig = plt.figure()
    ax = fig.gca()

    # first draw beam profile
    beam_xy = multivariate_normal(mean=[0,0], cov=args.spot_size/2)
    sz = 1.5*spot_size
    x = np.linspace(-sz/2, sz/2, 100)
    y = np.linspace(-sz/2, sz/2, 100)
    pos = np.dstack(np.meshgrid(x,y))
    intensity = n_particles*beam_xy.pdf(pos)
    plt.imshow(intensity, extent=2*[-sz/2, sz/2], cmap='plasma')
    plt.colorbar()

    # next add rectangles for each phone
    n_phones = x_s.size
    
    lim_x = RES_X * PIX_SIZE / 2
    lim_y = RES_Y * PIX_SIZE / 2
    rects = np.repeat(np.array([[
            [-lim_x, -lim_y], 
            [lim_x, -lim_y], 
            [lim_x, lim_y],
            [-lim_x, lim_y]]]), n_phones, axis=0)

    # do rotations

    cos_phi = np.cos(phi_s)
    sin_phi = np.sin(phi_s)
    phi_mat = np.array([[cos_phi, -sin_phi], [sin_phi, cos_phi]])

    theta_mat = np.array([
        [np.cos(theta_s), np.zeros(n_phones)],
        [np.zeros(n_phones), np.ones(n_phones)]
        ])

    cos_psi = np.cos(psi_s)
    sin_psi = np.sin(psi_s)
    psi_mat = np.array([[cos_psi, -sin_psi], [sin_psi, cos_psi]])


    # FIXME: there's a more elegant way to do this with einsum

    for iphone in range(args.n_phones):
         
        for corner in range(4):
            psi_rot = np.matmul(psi_mat[:,:,iphone], rects[iphone][corner])
            theta_rot = np.matmul(theta_mat[:,:,iphone], psi_rot)
            phi_rot = np.matmul(phi_mat[:,:,iphone], theta_rot)
            rects[iphone][corner] = phi_rot + np.array([x_s, y_s])[:,iphone]


    polygons = [Polygon(rects[iphone]) for iphone in range(n_phones)]
    p = PatchCollection(polygons, edgecolors='r', facecolors='none', linewidths=1)
    ax.add_collection(p)

    plt.show()


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='')
    parser.add_argument('--out', default='beam_sim.npz', help='Base file name, with a .npz extension')

    parser.add_argument('--n_phones', default=3, type=int, help='Number of phones')
    parser.add_argument('--fps', default=1., type=float, help='Frame rate')
    parser.add_argument('--drift', default=0, type=float, help='Standard deviation of drift (dt / t)')
    parser.add_argument('--phone_gap', default=9.2, type=float, help='Average distance between planes of phone sensors')
    parser.add_argument('--eff', default=1., type=float, help='Efficiency of phones')
    parser.add_argument('--upload_frac', default=1., type=float, help='Likelihood that a frame will successfully upload')

    parser.add_argument('--gap_stdev', default=.05, type=float, help='Standard deviation applied to phone_gap')
    parser.add_argument('--theta_dev', default=.06, type=float, help='Standard deviation in Eulerian theta for phone rotations')
    parser.add_argument('--psi_dev', default=.1, type=float, help='Standard deviation of Eulerian psi (about -phi) for phone rotations')
    parser.add_argument('--xy_dev', default=0.5, type=float, help='Standard deviation of phone positions in rotated plane with respect to apparatus')
    parser.add_argument('--xy_offset', default=0.8, type=float, help='Standard deviation of the phone holder position with respect to the beam')
    parser.add_argument('--theta_offset', default=.1, type=float, help='Standard deviation of the phone holder rotation')

    parser.add_argument('--n_spills', default=1, type=int, help='Number of beam spills')
    parser.add_argument('--spill_size', default=5e5, type=int, help='Number of particles in each beam spill')
    parser.add_argument('--spot_size', default=6., type=float, help='Diameter of beam in mm')
    parser.add_argument('--collimation', default=5e-5, type=float, help='Angular variation of beam particles in radians')

    parser.add_argument('--visualize', action='store_true', help='Create a visualization of the phone positions')
    args = parser.parse_args()

    basename = os.path.basename(args.out)
    dirname = os.path.dirname(args.out)
    
    os.makedirs(os.path.join(dirname, 'raw'), exist_ok=True)
    f_truth_name = basename.replace('.npz', '_truth.npz')
    spill_template = basename.replace('.npz', '_p{}_t{}.npz')
    millis = int(time.time() * 1000)

    # create the phone geometries
    n = args.n_phones

    phone_hwid = [uuid.uuid4().hex[:16] for _ in range(n)]

    phone_x0 = np.random.normal(0, args.xy_dev, n)
    phone_y0 = np.random.normal(0, args.xy_dev, n)
    phone_z0 = np.arange(n) * args.phone_gap \
            + np.random.normal(0, args.gap_stdev, n)
    
    # now move the holder about its center
    holder_xy = np.random.normal(0, args.xy_offset, 2)
    holder_theta = np.random.normal(0, args.theta_offset, 1)

    z_av = (n-1)/2 * args.phone_gap

    # N.B. while we are rotating the phone (and thus passively rotating
    # the coordinate system) the shift in phone positions looks like
    # an active rotation
    phone_x = phone_x0 * np.cos(holder_theta) - (phone_z0 - z_av) * np.sin(holder_theta) + holder_xy[0]
    phone_y = phone_y0 + holder_xy[1]
    phone_z = z_av + phone_x0 * np.sin(holder_theta) + (phone_z0 - z_av) * np.cos(holder_theta) 

    # create the angles within the holder
    phone_theta0 = np.abs(np.random.normal(0, args.theta_dev, n))
    phone_phi0 = np.random.uniform(-np.pi, np.pi, n)
    phone_psi0 = -phone_phi0 + np.random.normal(0, args.psi_dev, n) 
    
    # now find the angles after the holder rotation
    sth = np.sin(phone_theta0)
    cth = np.cos(phone_theta0)
    sps = np.sin(phone_psi0)
    cps = np.cos(phone_psi0)
    sph = np.sin(phone_phi0)
    cph = np.cos(phone_phi0)
    sho = np.sin(holder_theta)
    cho = np.cos(holder_theta)
    
    phone_theta = np.arccos(sth*sph*sho + cth*cho)
    phone_phi = np.where(np.logical_or(sth, sho), \
            np.arctan2(sth*sph*cho - cth*sho, sth*cph), phone_phi0)
    phone_psi = np.where(np.logical_or(sth, sho), \
            np.arctan2(cps*cph*sho - cth*sph*sps*sho + sps*sth*cho, \
            cps*sth*cho - sps*cph*sho - cth*sph*cps*sho), \
            phone_psi0)

    drifts = np.random.normal(loc=0, scale=args.drift, size=n)
    offsets = np.random.uniform(0, 1/args.fps, n)

    np.savez(os.path.join(dirname, f_truth_name), hwid=phone_hwid, x=phone_x, \
            y=phone_y, phi=phone_phi, theta=phone_theta, psi=phone_psi, \
            t_drifts=drifts, t_offsets=offsets, eff=args.eff)

    if args.visualize:
        visualize(args.spot_size, args.spill_size, phone_x, phone_y, \
                phone_phi, phone_theta, phone_psi)
        
    # now create the beam spill

    # create the beam times with an offset relative to the phones
    dt_spl = np.repeat(BEAM_PERIOD, args.n_spills - 1) + np.random.normal(0, 5)
    spl_times = np.insert(np.cumsum(dt_spl), 0, 0) \
            + np.random.uniform(0, BEAM_PERIOD)
 
    phone_times = offsets.copy()
    for ispill in range(args.n_spills):

        print('Spill {}/{}'.format(ispill+1, args.n_spills))

        n_particles = np.random.poisson(args.spill_size)
        dt_avg = BEAM_TIME / n_particles 
        dt = np.random.exponential(dt_avg, n_particles - 1)
        particle_times = np.insert(np.cumsum(dt), 0, 0) + spl_times[ispill]

        x0 = np.random.normal(0, args.spot_size/2, n_particles)
        y0 = np.random.normal(0, args.spot_size/2, n_particles)
        phi = np.random.uniform(0, 2*np.pi, n_particles)
        theta = np.abs(np.random.normal(0, args.collimation, n_particles))

        xs, ys = map_coords(x0, y0, phi, theta, phone_x, phone_y, phone_z, \
                phone_phi, phone_theta, phone_psi)

        pix_x = xs / PIX_SIZE
        pix_y = ys / PIX_SIZE

        # determine the phone response  

        for iphone in range(args.n_phones):
            
            tmin = phone_times[iphone] - 1/args.fps
            while phone_times[iphone] < particle_times[0]:
                tmin = phone_times[iphone]
                phone_times[iphone] += 1/args.fps

            while tmin < particle_times[-1]: 
                frame_particles = (particle_times > tmin) & (particle_times < phone_times[iphone])

                xframe = pix_x[iphone, frame_particles]
                yframe = pix_y[iphone, frame_particles]

                hits = (np.abs(xframe) < RES_X/2) \
                        & (np.abs(yframe) < RES_Y/2) \
                        & (np.random.random(xframe.shape) < args.eff)

                xhits = (xframe[hits] + RES_X/2).astype(int)
                yhits = (yframe[hits] + RES_Y/2).astype(int)

                # save output 
                if np.random.random(0) > args.upload_frac: continue
                t_out = millis + int(1000 * (phone_times[iphone] + \
                        (phone_times[iphone]-offsets[iphone]) * drifts[iphone]))
                f_spill_name = spill_template.format(phone_hwid[iphone], t_out)
                fname = os.path.join(dirname, 'raw', f_spill_name)
                np.savez(fname, x=xhits, y=yhits, t=t_out, shape=np.array([RES_X, RES_Y]))

                tmin = phone_times[iphone]
                phone_times[iphone] += 1 / args.fps

