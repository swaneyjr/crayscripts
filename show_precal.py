import crayfis_data_pb2 as pb
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cv2

def get_precal(fname):
    dc = pb.DataChunk()
    f = open(fname)
    dc.ParseFromString(f.read())
    f.close()

    return dc.precalibration_results[0]

def get_weight_array(precal, size=None)
   
    res_x = precal.sample_res_x
    res_y = precal.sample_res_y
    downsample = np.array(precal.weights).reshape(res_y,res_x)
    if size:
        resize = cv2.INTER_CUBIC
        if precal.interpolation:
            resize = precal.interpolation
	return cv2.resize(downsample, size, interpolation=resize)
    return downsample


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('precal', help='protobuf .bin file')
    parser.add_argument('--resample', help='Comma-delimited size to rescale')
    args = parser.parse_args()

    if args.resample:
	args.resample = tuple(map(int, args.resample.split(',')))

    matplotlib.use('tkagg')
    plt.figure(1)
    plt.imshow(get_precal_array(args.precal, args.resample), cmap='RdYlBu', interpolation='nearest')
    plt.colorbar()

    plt.show()
    
