import numpy as np
import show_precal

def get_correlation(infiles):

    # initialize arrays
    
    precal = show_precal.get_precal(infiles[0])
    weights = show_precal.get_weight_array(precal)
    temp = precal.battery_temp

    exp_w = weights
    exp_w_sq = weights**2
    exp_t = temp
    exp_t_sq = temp**2
    exp_wt = weights*temp

    # loop over rest of files
    
    for pathname in infiles[1:]:

        precal = show_precal.get_precal(pathname)
        weights = show_precal.get_weight_array(precal)
        temp = precal.battery_temp
            
        exp_w += weights
        exp_w_sq += weights**2
        exp_t += temp
        exp_t_sq += temp**2
        exp_wt += weights*temp

    # calculate correlation coefficients

    n_frames = len(infiles)
    exp_w /= n_frames
    exp_w_sq /= n_frames
    exp_t /= n_frames
    exp_t_sq /= n_frames
    exp_wt /= n_frames

    cov = exp_wt - exp_w * exp_t
    sig_w = np.sqrt(exp_w_sq - exp_w**2)
    sig_t = (exp_t_sq - exp_t**2)**0.5

    return cov/(sig_w * sig_t)
        

if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    from argparse import ArgumentParser

    parser = ArgumentParser(description='')
    parser.add_argument('--in', nargs='+', help='')
    args = parser.parse_args()

    matplotlib.use('tkagg')
    plt.figure(1)
    plt.imshow(get_correlation(args.in), cmap='coolwarm', interpolation='nearest')
    plt.colorbar()
    plt.show()
