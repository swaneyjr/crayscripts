#!/usr/bin/env python

import crayfis_data_pb2 as pb
import numpy as np
import cv2

def get_precal(fname):
    dc = pb.DataChunk()
    f = open(fname)
    dc.ParseFromString(f.read())
    f.close()

    return dc.precalibration_results[0]

def get_weight_array(precal, size=None, compressed=True):
   
    res_x = precal.sample_res_x
    res_y = precal.sample_res_y
    
    if compressed:
	compressed_array = np.array(bytearray(precal.compressed_weights))

        downsample = cv2.imdecode(compressed_array, 0)/255.
    else:
        downsample = np.array(precal.weights).reshape(res_y,res_x)
    
    if size:
        if precal.interpolation:
            return cv2.resize(downsample, size, interpolation=precal.interpolation)
        return cv2.resize(downsample, size)
    return downsample


if __name__ == '__main__':
    from argparse import ArgumentParser
    import matplotlib
    import matplotlib.pyplot as plt
    
    parser = ArgumentParser()
    parser.add_argument('precal', help='protobuf .bin file')
    parser.add_argument('--resample', help='Comma-delimited size to rescale')
    compression = parser.add_mutually_exclusive_group(required=True)
    compression.add_argument('-u','--uncompressed', action='store_true', help='Use uncompressed weights')
    compression.add_argument('-c','--compressed', action='store_true', help='Use compressed weights')
    compression.add_argument('-x','--compare', action='store_true', help='Find ratio of compressed/uncompressed weights')
    args = parser.parse_args()

    if args.resample:
	args.resample = tuple(map(int, args.resample.split(',')))

    matplotlib.use('tkagg')
    
    if not args.compare:
        plt.figure(1)
        plt.imshow(get_weight_array(get_precal(args.precal), args.resample, args.compressed), cmap='plasma', interpolation='nearest')
        plt.colorbar()
    else:
        precal = get_precal(args.precal)

        compressed = get_weight_array(precal, args.resample, True)
        uncompressed = get_weight_array(precal, args.resample, False)
        
        plt.figure(1)
        plt.title('Uncompressed')
        plt.imshow(uncompressed, cmap='plasma', interpolation='nearest')
        plt.colorbar()

        plt.figure(2)
        plt.title('Compressed')
        plt.imshow(uncompressed, cmap='plasma', interpolation='nearest')
        plt.colorbar()

        plt.figure(3)
        plt.title('Compressed/Uncompressed')
        plt.imshow(compressed/uncompressed, cmap='coolwarm', interpolation='nearest')
        plt.colorbar()
        

    plt.show()
    
