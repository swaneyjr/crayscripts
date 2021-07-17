#!/usr/bin/env python

import crayfis_data_pb2 as pb
import numpy as np
import cv2

def get_precal(fname):
    dc = pb.DataChunk()
    f = open(fname, 'rb')
    dc.ParseFromString(f.read())
    f.close()

    return dc.precalibration_results[0]

def get_weight_array(precal, size=None, compressed=True):
    
    if compressed:
        compressed_array = np.array(bytearray(precal.compressed_weights))

        downsample = cv2.imdecode(compressed_array, 0)/255.
    else:
        downsample = np.array(precal.weights)
    
    if size:
        if precal.interpolation:
            return cv2.resize(downsample, size, interpolation=precal.interpolation).T
        return cv2.resize(downsample, size).T
    return downsample.T


if __name__ == '__main__':
    from argparse import ArgumentParser
    import matplotlib.pyplot as plt
    
    parser = ArgumentParser()
    parser.add_argument('precal', help='protobuf .bin file')
    parser.add_argument('--resample', type=int, nargs=2, help='Comma-delimited size to rescale')
    compression = parser.add_mutually_exclusive_group(required=True)
    compression.add_argument('-u','--uncompressed', action='store_true', help='Use uncompressed weights')
    compression.add_argument('-c','--compressed', action='store_true', help='Use compressed weights')
    compression.add_argument('-x','--compare', action='store_true', help='Find ratio of compressed/uncompressed weights')
    args = parser.parse_args()
    

    resample = (args.resample[1], args.resample[0]) if args.resample else None
    
    if not args.compare:
        plt.figure(figsize=(5,2.5), tight_layout=True)
        precal = get_precal(args.precal)
        weights = get_weight_array(precal, resample, args.compressed)
        plt.imshow(1/weights, interpolation='nearest', origin='lower', 
                cmap='jet', vmin=1, extent=[0,precal.res_x,0,precal.res_y])
        cbar = plt.colorbar()
        cbar.set_label('Lens-shading gains', rotation=270, labelpad=16)
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
    else:
        precal = get_precal(args.precal) 

        compressed = get_weight_array(precal, resample, True)
        uncompressed = get_weight_array(precal, resample, False)
        
        plt.figure(1)
        plt.title('Uncompressed')
        plt.imshow(uncompressed, cmap='plasma', interpolation='nearest', origin='lower')
        plt.colorbar()

        plt.figure(2)
        plt.title('Compressed')
        plt.imshow(uncompressed, cmap='plasma', interpolation='nearest', origin='lower')
        plt.colorbar()

        plt.figure(3)
        plt.title('Compressed/Uncompressed')
        plt.imshow(compressed/uncompressed, cmap='coolwarm', interpolation='nearest', origin='lower')
        plt.colorbar()
        

    plt.show()
    
