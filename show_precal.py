import crayfis_data_pb2 as pb
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cv2

def get_precal_array(fname, size=None):
    dc = pb.DataChunk()
    f = open(fname)
    dc.ParseFromString(f.read())
    f.close()
   
    precal = dc.precalibration_results[0]
    res_x = precal.sample_res_x
    res_y = precal.sample_res_y
    downsample = np.array(precal.weights).reshape(res_y,res_x)
    if size:
	return cv2.resize(downsample, size, interpolation=cv2.INTER_CUBIC)
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
    plt.imshow(get_precal_array(args.precal, args.resample), cmap='plasma', interpolation='nearest')
    plt.colorbar()

    plt.show()
    
