#!/usr/bin/env python3

import os
import time
import uuid

import numpy as np
from scipy.stats import multivariate_normal
from sklearn.cluster import DBSCAN

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

def gen_noise(n_avg):
    global RES_X, RES_Y
    
    # add noise in a radial lens-shade pattern
    n_noise = np.random.poisson(n_avg)
    x_noise = np.array([])
    y_noise = np.array([])

    diag = np.sqrt(RES_X**2 + RES_Y**2) / 2
                
    while x_noise.size < n_noise:
        r_n = diag * np.random.rand(n_noise)**(1/3)
        phi_n = np.random.uniform(size=n_noise, high=2*np.pi)

        x_new = r_n * np.cos(phi_n)
        y_new = r_n * np.sin(phi_n)

        cut = (np.abs(x_new) < RES_X/2) & (np.abs(y_new) < RES_Y/2)

        x_new = (x_new[cut] + RES_X/2).astype(int)
        y_new = (y_new[cut] + RES_Y/2).astype(int)

        x_noise = np.append(x_new, x_noise)
        y_noise = np.append(y_new, y_noise)

    return x_noise[:n_noise], y_noise[:n_noise]


def visualize(spot_type, spot_size, n_particles, x_s, y_s, phi_s, theta_s, psi_s):

    global RES_X, RES_Y, PIX_SIZE

    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    fig = plt.figure()
    ax = fig.gca()

    # first draw beam profile 
    sz = 1.5*spot_size
    x = np.linspace(-sz/2, sz/2, 100)
    y = np.linspace(-sz/2, sz/2, 100)
    if spot_type == 'gaussian':
        pos = np.dstack(np.meshgrid(x,y))
        beam_xy = multivariate_normal(mean=[0,0], cov=spot_size/2)
        intensity = n_particles*beam_xy.pdf(pos)
    elif spot_type == 'uniform':
        xx, yy = np.meshgrid(x,y)
        intensity = n_particles / (np.pi * spot_size**2 / 4) *  (xx**2 + yy**2 < spot_size**2 / 4)
    elif spot_type == 'square':
        xx, yy = np.meshgrid(x,y)
        intensity = 4 * n_particles / spot_size**2 * (np.maximum(np.abs(xx), np.abs(yy)) < spot_size / 2)
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
    parser.add_argument('--out', required=True, help='Base file name, with a .npz extension')

    parser.add_argument('--n_phones', default=3, type=int, help='Number of phones')
    parser.add_argument('--fps', default=1., type=float, help='Frame rate')
    parser.add_argument('--drift', default=0, type=float, help='Standard deviation of drift (dt / t)')
    parser.add_argument('--noise', default=0, type=float, help='Average number of background pixels per image')
    parser.add_argument('--bgfiles', default=250, type=int, help='Number of files per phone generated with background noise only')
    parser.add_argument('--phone_gap', default=9.7, type=float, help='Average distance between planes of phone sensors')
    parser.add_argument('--particles', default=[1.], type=float, nargs='+', help='Abundances of different particle species in the beam spill.  Arguments should add to 1')
    parser.add_argument('--eff', default=[1.], type=float, nargs='+', help='Efficiency of phones to different particle species.  Number of arguments should be the same as n_particles.')
    parser.add_argument('--upload_frac', default=1., type=float, help='Likelihood that a frame will successfully upload')



    parser.add_argument('--gap_stdev', default=0.5, type=float, help='Standard deviation applied to phone_gap')
    parser.add_argument('--theta_dev', default=.05, type=float, help='Standard deviation in Eulerian theta for phone rotations')
    parser.add_argument('--psi_dev', default=.05, type=float, help='Standard deviation of Eulerian psi (about -phi) for phone rotations')
    parser.add_argument('--xy_dev', default=0.5, type=float, help='Standard deviation of phone positions (in mm) in rotated plane with respect to apparatus')
    parser.add_argument('--xy_offset', default=0.8, type=float, help='Standard deviation of the phone holder position (in mm) with respect to the beam')
    parser.add_argument('--theta_offset', default=.05, type=float, help='Standard deviation of the phone holder rotation')

    parser.add_argument('--n_spills', default=1, type=int, help='Number of beam spills')
    parser.add_argument('--spill_size', default=1e6, type=int, help='Number of particles in each beam spill')
    parser.add_argument('--spot_size', default=10., type=float, help='Diameter of beam in mm')
    parser.add_argument('--spot_type', choices=('gaussian', 'uniform', 'square'), default='gaussian', help='Distribution of beam profile')
    parser.add_argument('--collimation', default=5e-5, type=float, help='Angular variation of beam particles in radians')

    parser.add_argument('--visualize', action='store_true', help='Create a visualization of the phone positions')
    parser.add_argument('--clustering', type=float, default=-1, help='Maximum pixel separation to be considered a cluster.  Default of -1, in which case, pixel multiplicities are included.')
    
    args = parser.parse_args()

    if np.abs(sum(args.particles) - 1) > 1e-5:
        print('ERROR: relative particle abundances not normalized.')
        resp = input('Fix normalization? (y/n) ')
        if resp.strip().lower() == 'y':
            args.particles /= sum(args.particles)
        else:
            print('Exiting')
            exit(0)

    if not len(args.eff) == len(args.particles):
        print('ERROR: Invalid number of entries for --eff: must be equal to n_particles')
        exit(1)

    basename = os.path.basename(args.out)
    dirname = os.path.dirname(args.out)
    
    os.makedirs(os.path.join(dirname, 'cluster'), exist_ok=True)
    os.makedirs(os.path.join(dirname, 'bg'), exist_ok=True)

    f_truth_name = basename.replace('.npz', '_truth.npz')
    f_config_name = basename.replace('.npz', '_config.npz')
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


    if args.visualize:
        visualize(args.spot_type, args.spot_size, args.spill_size, \
                phone_x, phone_y, phone_phi, phone_theta, phone_psi)
        if not input('Continue? y/n: ') == 'y': exit(0)

    
    drifts = np.random.normal(loc=0, scale=args.drift, size=n)
    offsets = np.random.uniform(0, 1/args.fps, n)

    # handle particle species and efficiencies
    particle_cumsum = np.cumsum(args.particles).reshape(-1,1)
    efficiencies = np.array(args.eff)

    # save the config and truth files
    np.savez(os.path.join(dirname, f_config_name), fps=args.fps, \
            res=np.array([RES_X, RES_Y]), particles=args.particles)

    np.savez(os.path.join(dirname, f_truth_name), hwid=phone_hwid, x=phone_x, \
            y=phone_y, phi=phone_phi, theta=phone_theta, psi=phone_psi, \
            t_drifts=drifts, t_offsets=offsets, eff=efficiencies)
    
    # generate noise-only files
    print('Generating background files')
    for hwid in phone_hwid:
        for i in range(args.bgfiles):
            x_noise, y_noise = gen_noise(args.noise)
            np.savez(os.path.join(dirname, 'bg', '{}_{}'.format(hwid, i)), \
                    x=x_noise, y=y_noise)

    # now create the beam spills

    # create the beam times with an offset relative to the phones
    dt_spl = np.repeat(BEAM_PERIOD, args.n_spills - 1) + np.random.normal(0, 5)
    spl_times = np.insert(np.cumsum(dt_spl), 0, 0) \
            + np.random.uniform(0, BEAM_PERIOD)

    phone_times = offsets.copy()
    for ispill in range(args.n_spills):

        print('Spill {}/{}'.format(ispill+1, args.n_spills), end="\r")

        n_particles = np.random.poisson(args.spill_size)
        dt_avg = BEAM_TIME / n_particles 
        dt = np.random.exponential(dt_avg, n_particles - 1)
        particle_times = np.insert(np.cumsum(dt), 0, 0) + spl_times[ispill]

        if args.spot_type == 'gaussian':
            x0 = np.random.normal(0, args.spot_size/2, n_particles)
            y0 = np.random.normal(0, args.spot_size/2, n_particles)
        elif args.spot_type == 'uniform':
            r = np.random.triangular(0, args.spot_size/2, args.spot_size/2, n_particles)
            phi = np.random.uniform(0, 2 * np.pi, n_particles)
            x0 = r * np.cos(phi)
            y0 = r * np.sin(phi)
        elif args.spot_type == 'square':
            x0 = np.random.uniform(-args.spot_size/2, args.spot_size/2, n_particles)
            y0 = np.random.uniform(-args.spot_size/2, args.spot_size/2, n_particles)
        else:
            print('Unrecognized beam profile. Exiting.')
            exit()

        phi = np.random.uniform(0, 2*np.pi, n_particles)
        theta = np.abs(np.random.normal(0, args.collimation, n_particles))

        xs, ys = map_coords(x0, y0, phi, theta, phone_x, phone_y, phone_z, \
                phone_phi, phone_theta, phone_psi)

        pix_x = xs / PIX_SIZE
        pix_y = ys / PIX_SIZE

        particle_eff = efficiencies[np.argmin(particle_cumsum < np.random.random(n_particles), axis=0)]

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
                        & (np.random.random(xframe.shape) < particle_eff[frame_particles])

                x_hits = (xframe[hits] + RES_X/2).astype(int)
                y_hits = (yframe[hits] + RES_Y/2).astype(int)
                
                # add noise
                x_noise, y_noise = gen_noise(args.noise)

                xtot = np.hstack([x_hits, x_noise])
                ytot = np.hstack([y_hits, y_noise])

                if args.clustering == 0:
                    # just get rid of duplicates
                    itot = np.unique(xtot + RES_X * ytot)
                
                    xtot = itot % RES_X
                    ytot = itot // RES_X
                
                if args.clustering > 0:
                    # use an actual clustering algorithm
                    xy = np.column_stack((xtot,ytot))
                    clustering = DBSCAN(eps=args.clustering, min_samples=1)
                    clustering.fit(xy)
                    
                    n_groups = clustering.labels_.max() + 1
                    group_indices = np.array([np.argmax(clustering.labels_ == i) for i in range(n_groups)])
                    
                    xtot = xtot[group_indices]
                    ytot = ytot[group_indices]
                    

                # save output 
                if np.random.random() > args.upload_frac: continue
                t_out = millis + 1000 * (phone_times[iphone] + \
                        (phone_times[iphone]-offsets[iphone]) * drifts[iphone])
                f_spill_name = spill_template.format(phone_hwid[iphone], int(t_out))
                fname = os.path.join(dirname, 'cluster', f_spill_name)
                np.savez(fname, 
                        x=xtot, 
                        y=ytot, 
                        t=t_out, 
                        # truth information
                        particles=np.arange(n_particles)[frame_particles][hits],
                        #ptimes=1000*particle_times[frame_particles][hits] + millis,
                        )

                tmin = phone_times[iphone]
                phone_times[iphone] += 1 / args.fps

    print('Spill {0}/{0}'.format(args.n_spills))
