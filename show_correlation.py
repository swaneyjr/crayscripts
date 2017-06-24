import numpy as np
import show_precal

def get_stats(infiles):

    # initialize arrays
    
    precal = show_precal.get_precal(infiles[0])
    weights = show_precal.get_weight_array(precal)
    temp = precal.battery_temp

    sum_w = weights
    sum_w_sq = weights**2
    sum_t = temp
    sum_t_sq = temp**2
    sum_wt = weights*temp

    # loop over rest of files
    
    for pathname in infiles[1:]:

        precal = show_precal.get_precal(pathname)
        weights = show_precal.get_weight_array(precal)
        temp = precal.battery_temp
            
        sum_w += weights
        sum_w_sq += weights**2
        sum_t += temp
        sum_t_sq += temp**2
        sum_wt += weights*temp

    # calculate correlation coefficients

    n_frames = len(infiles)
    exp_w = sum_w/n_frames
    exp_w_sq = sum_w_sq/n_frames

    sig_w = np.sqrt(exp_w_sq - exp_w**2)

    rho = (n_frames*sum_wt - sum_w*sum_t) / (n_frames*sig_w*(n_frames*sum_t_sq - sum_t**2)**0.5)
    return exp_w, sig_w, rho
        

if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Converts PrecalibrationResult.bin files into array of correlation coefficients with temperature and weights')
    parser.add_argument('--in', dest='infiles', nargs='+', help='.bin files to use')
    args = parser.parse_args()

    mean, std, corr = get_stats(args.infiles)

    matplotlib.use('tkagg')
    plt.figure(1)
    plt.title('Mean Weights')
    plt.imshow(mean, cmap='plasma', interpolation='nearest')
    plt.colorbar()

    plt.figure(2)
    plt.title('Weight Standard Deviations')
    plt.imshow(std, cmap='viridis', interpolation='nearest')
    plt.colorbar()

    plt.figure(3)
    plt.title('Weight Correlation Coefficients')
    plt.imshow(corr, cmap='coolwarm', interpolation='nearest')
    plt.colorbar()


    plt.show()
