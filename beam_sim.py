#!/usr/bin/env python3

import time

import numpy as np
from scipy.sparse import coo_matrix, save_npz
from scipy.stats import multivariate_normal

BEAM_TIME = 4.2 # seconds
BEAM_PERIOD = 60 # seconds
PIX_SIZE = .00112 # mm

RES_X = 4656
RES_Y = 3492

def map_coords(xp, yp, phi_p, theta_p, zs, phi_s, theta_s, psi_s, xs0, ys0):
    """
    Parameters:

    xp, yp:                 initial position of the particle in the beam frame
    phi_p, theta_p:         trajectory of the particle with respect to the 
                            beam axis
    zs:                     distance from particle's starting coordinates to 
                            the plane of the beam, as measured along the beam 
                            axis
    phi_s, theta_s, psi_s:  Euler angles of the sensor plane
    xs0, ys0:                 Coordinates of the center of the sensor in the
                            sensor plane, relative to its intersection with
                            the beam axis

    Returns:
    
    xs, ys:                 coordinates of the particle hit relative to the
                            center of the sensor
    """

    # convert to polar coords
    rp = np.sqrt(xp**2 + yp**2)
    xi_p = np.arctan2(yp, xp)

    # calculate trig functions
    cos_xi = np.cos(xi_p - phi_s)
    sin_xi = np.sin(xi_p -phi_s)

    cos_phi = np.cos(phi_p - phi_s)
    sin_phi = np.sin(phi_p - phi_s)

    cos_psi = np.cos(psi_s)
    sin_psi = np.sin(psi_s)

    sec_theta = 1 / np.cos(theta_s)
    
    denom = 1 / np.tan(theta_p) + np.tan(theta_s) * cos_phi

    # result
    xs = -xs0 + rp * (cos_xi *cos_psi*sec_theta + sin_xi *sin_psi) \
              + zs * (cos_phi*cos_psi*sec_theta + sin_phi*sin_psi) / denom
    ys = -ys0 + rp * (sin_xi *cos_psi - cos_xi *sec_theta*sin_psi) \
              + zs * (sin_phi*cos_psi - cos_phi*sec_theta*sin_psi) / denom

    return xs, ys


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='')
    parser.add_argument('--out', default='beam_sim.npz')

    parser.add_argument('--n_phones', default=3, type=int)
    parser.add_argument('--fps', default=1., type=float)
    parser.add_argument('--phone_gap', default=10., type=float)
    parser.add_argument('--eff', default=1., type=float)

    parser.add_argument('--gap_stdev', default=1., type=float)
    parser.add_argument('--theta_dev', default=.05, type=float)
    parser.add_argument('--psi_dev', default=.05, type=float)
    parser.add_argument('--xy_dev', default=1, type=float)

    parser.add_argument('--n_spills', default=1, type=int)
    parser.add_argument('--spill_size', default=5e5, type=int)
    parser.add_argument('--spot_size', default=6., type=float)
    parser.add_argument('--collimation', default=.001, type=float)

    parser.add_argument('--visualization', action='store_true')
    args = parser.parse_args()

    fname = args.out.replace('.npz', '_p{}_t{}.npz')
    millis = int(time.time() * 1000)

    # create the phone geometries
    n = (args.n_phones, 1)

    phone_z = np.arange(args.n_phones).reshape(n) * args.phone_gap \
            + np.random.normal(0, args.gap_stdev, n)
    phone_phi = np.random.uniform(0, 2*np.pi, n)
    phone_psi = -phone_phi + np.random.normal(0, args.psi_dev, n)
    phone_theta = np.abs(np.random.normal(0, args.theta_dev, n))

    phone_x = np.random.normal(0, args.xy_dev, n)
    phone_y = np.random.normal(0, args.xy_dev, n)

    beam_xy = multivariate_normal(mean=[0,0], cov=args.spot_size/2)

    if args.visualization:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection

        fig = plt.figure()
        ax = fig.gca()

        # first draw beam profile
        sz = 1.5*args.spot_size
        x = np.linspace(-sz/2, sz/2, 100)
        y = np.linspace(-sz/2, sz/2, 100)
        pos = np.dstack(np.meshgrid(x,y))
        intensity = args.spill_size*beam_xy.pdf(pos)
        plt.imshow(intensity, extent=2*[-sz/2, sz/2], cmap='plasma')
        plt.colorbar()

        # next add rectangles for each phone
        lim_x = RES_X * PIX_SIZE / 2
        lim_y = RES_Y * PIX_SIZE / 2
        rects = np.repeat(np.array([[
                [-lim_x, -lim_y], 
                [lim_x, -lim_y], 
                [lim_x, lim_y],
                [-lim_x, lim_y]]]), args.n_phones, axis=0)

        # apply translation
        rects += np.dstack([phone_x, phone_y]).reshape(args.n_phones,1,2)

        # do rotations
        # FIXME: there's a more elegant way to do this
        for iphone in range(args.n_phones):
            cos_psi = np.cos(phone_psi[iphone,0])
            sin_psi = np.sin(phone_psi[iphone,0])
            psi_mat = np.array([[cos_psi, -sin_psi], [sin_psi, cos_psi]])

            theta_mat = np.array([[np.cos(phone_theta[iphone,0]), 0],[0,1]])

            cos_phi = np.cos(phone_phi[iphone,0])
            sin_phi = np.sin(phone_phi[iphone,0])
            phi_mat = np.array([[cos_phi, -sin_phi], [sin_phi, cos_phi]])

            for corner in range(4):
                psi_rot = np.matmul(psi_mat, rects[iphone][corner])
                theta_rot = np.matmul(theta_mat, psi_rot)
                phi_rot = np.matmul(phi_mat, theta_rot)
                rects[iphone][corner] = phi_rot


        polygons = [Polygon(rects[iphone]) for iphone in range(args.n_phones)]
        p = PatchCollection(polygons, edgecolors='r', facecolors='none', linewidths=1)
        ax.add_collection(p)

        plt.show()

    # now create the beam spill
    for ispill in range(args.n_spills):

        print('Spill {}/{}'.format(ispill+1, args.n_spills))

        n_particles = np.random.poisson(args.spill_size)
        dt_avg = BEAM_TIME / n_particles 
        dt = np.random.exponential(dt_avg, n_particles - 1)
        times = np.insert(np.cumsum(dt), 0, 0)

        x0 = np.random.normal(0, args.spot_size/2, n_particles)
        y0 = np.random.normal(0, args.spot_size/2, n_particles)
        phi = np.random.uniform(0, 2*np.pi, n_particles)
        theta = np.abs(np.random.normal(0, args.collimation, n_particles))

        xs, ys = map_coords(x0, y0, phi, theta, phone_z, phone_phi, \
                phone_theta, phone_psi, phone_x, phone_y)

        pix_x = xs / PIX_SIZE
        pix_y = ys / PIX_SIZE

        # determine the phone response
        offsets = np.random.uniform(0, 1/args.fps, args.n_phones)

        for iphone in range(args.n_phones):
            
            tmin = 0
            tmax = offsets[iphone]

            while tmin < times[-1]: 
                frame_particles = (times > tmin) & (times < tmax)

                xframe = pix_x[iphone, frame_particles]
                yframe = pix_y[iphone, frame_particles]

                hits = (np.abs(xframe) < RES_X/2) \
                        & (np.abs(yframe) < RES_Y/2) \
                        & (np.random.random(xframe.shape) < args.eff)

                xhits = (xframe[hits] + RES_X/2).astype(int)
                yhits = (yframe[hits] + RES_Y/2).astype(int)
                data = np.ones(xhits.shape, dtype=bool)


                sparse_hits = coo_matrix((data, (yhits, xhits)), \
                        shape=(RES_Y, RES_X))

                # save output 
                t_abs = millis + int(1000*(BEAM_PERIOD*ispill + tmin))
                save_npz(fname.format(iphone, t_abs), sparse_hits)

                tmin = tmax
                tmax += 1 / args.fps

