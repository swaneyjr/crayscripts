import crayfis_data_pb2 as pb
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def get_precal_array(fname):
    dc = pb.DataChunk()
    f = open(fname)
    dc.ParseFromString(f.read())
    f.close()
   
    precal = dc.precalibration_results[0]
    res_x = precal.sample_res_x
    res_y = precal.sample_res_y
    return np.array(precal.weights).reshape(res_y,res_x) 


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('precal', help='protobuf .bin file')
    args = parser.parse_args()

    matplotlib.use('tkagg')
    plt.figure(1)
    plt.imshow(get_precal_array(args.precal), cmap='plasma', interpolation='nearest')
    plt.colorbar()

    plt.show()
    
