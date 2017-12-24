from show_precal import get_precal, get_weight_array
import numpy as np

def get_thresh_array(fname, thresh):
    precal = get_precal(fname)
    res_x = precal.res_x
    res_y = res_x * precal.sample_res_y / precal.sample_res_x
    weights = get_weight_array(precal, size=(res_x, res_y))

    return np.floor((thresh+0.5)/weights)



if __name__ == '__main__':
    from argparse import ArgumentParser
    import matplotlib
    import matplotlib.pyplot as plt

    parser = ArgumentParser(description='Find the spatially dependent thresholds equivalent to the calculated weights')
    parser.add_argument('precal', help='PreCalibration.bin file to process')
    parser.add_argument('--thresh', required=True, type=int, help='L1/L2 threshold used on adjusted pix_vals')

    args = parser.parse_args()

    plt.figure(1)
    plt.imshow(get_thresh_array(args.precal, args.thresh), cmap='jet', interpolation='nearest')
    plt.colorbar()

    plt.title('Spatially dependent thresholds')
    plt.show()    

