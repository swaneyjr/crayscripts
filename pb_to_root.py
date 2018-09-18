#!/usr/bin/env python2

import crayfis_data_pb2 as pb
import ROOT as r
import numpy as np
from collections import namedtuple
from os import remove

type_set = namedtuple('type_set', ['full','short','numpy'])

def merge_pb(infiles):
    dc = pb.DataChunk()
    for fname in infiles:
        f = open(fname)
        dc.MergeFromString(f.read())
        f.close()
        remove(fname)
    merge_name = infiles[0].split('_')[-2] + '.bin'
    print "Compiling messages into %s" % merge_name
    merge_file = open(merge_name, "wb")
    merge_file.write(dc.SerializeToString())
    merge_file.close()
    print "Messages saved."
    return dc

def pb_to_hists(dc):
    xb_hist_tot = np.zeros(256)
    evt_hist_tot = np.zeros(256)

    for xb in dc.exposure_blocks:
        xb_hist = np.array(xb.hist)
        xb_hist.resize(256)
        xb_hist_tot += xb_hist

        for evt in xb.events:
            evt_hist = np.array(evt.hist)
            evt_hist.resize(256)
            evt_hist_tot += evt_hist

            for pix in evt.pixels:
                if pix.val > xb.L1_thresh:
                    xb_hist_tot[pix.adjusted_val] += 1
                    evt_hist_tot[pix.adjusted_val] += 1

    L1_pass_hist = r.TH1F('L1pass','Pass', 256, 0, 256)
    L1_skip_hist = r.TH1F('L1skip','Skip', 256, 0, 256)
    
    for val in xrange(256):
        L1_pass_hist.SetBinContent(val, evt_hist_tot[val])
        L1_pass_hist.SetBinError(val, evt_hist_tot[val]**0.5)

        L1_skip_hist.SetBinContent(val, xb_hist_tot[val] - evt_hist_tot[val])
        L1_skip_hist.SetBinError(val, (xb_hist_tot[val] + evt_hist_tot[val])**0.5)
    
    return L1_pass_hist, L1_skip_hist
        
def pb_to_precal(dc):

    for precal in dc.precalibration_results:

        weight_hist = r.TH2C(precal.run_id, 'Precalibration weights', precal.sample_res_x, 0, precal_sample_res_x, precal.sample_res_y, 0, precal_sample_res_y)
        
        compressed_weights = np.array(bytearray(precal.compressed_weights))
        uncompressed = cv2.imdecode(compressed_weights, 0)

        for ixy, weight in np.nditer(uncompressed):
            weight_hist.SetBinAt(ixy[0], ixy[1], weight)

def pb_to_trees(dc):
    
    # write to ROOT
    exposure = r.TTree('exposure', 'Exposure from protobuf data')
    events = r.TTree('events', 'Events from protobuf data')

    
    type_name = {int: type_set('UInt_t', '/i', np.uint32,), \
	long: type_set('ULong64_t', '/l', np.uint64), \
	float: type_set('Double_t', '/D', np.float64), \
	bool: type_set('Bool_t', '/O', np.bool_), \
	unicode: type_set('Char_t', '/C', unicode) }
    

    # assign appropriate branches to trees

    xb_containers = {}
    evt_containers = {}
    pix_containers = {}
    pix_n = np.zeros(1, dtype=int)
    events.Branch('pix_n', pix_n, 'pix_n/I')    
    
    # fill tree
    print "Filling tree"
    
    total_xb = len(dc.exposure_blocks)
    for ixb,xb in enumerate(dc.exposure_blocks):

        if ixb % (total_xb/10) == 0:
            print " %d/%d (%0.2f)" % (ixb, total_xb, 100.*ixb/total_xb)
        
        for xb_field,val in xb.ListFields():
            if xb_field.label == 3: continue
            if not xb_field.name in xb_containers.iterkeys(): 
                xb_containers[xb_field.name] = np.zeros(1,dtype=type_name[type(val)].numpy)
                exposure.Branch(xb_field.name, xb_containers[xb_field.name], \
                                xb_field.name+type_name[type(val)].short)
            xb_containers[xb_field.name][0] = val
	    if type(val) == unicode:
                xb_containers[xb_field.name][0] += '\0'
        exposure.Fill()
        
        for evt in xb.events:
            pix_n[0] = 0
            for evt_field,val in evt.ListFields():
                if evt_field.label == 3 or evt_field.type == 11: continue
                if not evt_field.name in evt_containers.iterkeys():
                    evt_containers[evt_field.name] = np.zeros(1,dtype=type_name[type(val)].numpy)
                    events.Branch(evt_field.name, evt_containers[evt_field.name], \
                                  evt_field.name + type_name[type(val)].short)
                evt_containers[evt_field.name][0] = val
                if type(val) == unicode:
                    evt_containers[evt_field.name][0] += '\0'
            
            evt_containers['timestamp'][0] = xb.start_time
            # reset pixel vectors
            for v in pix_containers.itervalues():
                v.clear()
            
            if evt.zero_bias.val:
                if not 'zb_x' in pix_containers.iterkeys():
                    pix_containers['zb_x'] = r.vector('UInt_t')()
                    pix_containers['zb_y'] = r.vector('UInt_t')()
                    pix_containers['zb_val'] = r.vector('UInt_t')()
                    events.Branch('zb_x', pix_containers['zb_x'])
                    events.Branch('zb_y', pix_containers['zb_y'])
                    events.Branch('zb_val', pix_containers['zb_val'])


                s = int(len(evt.zero_bias.val)**0.5)

                for dy in xrange(s):
                    for dx in xrange(s):
                        pix_containers['zb_x'].push_back(evt.zero_bias.x_min + dx)
                        pix_containers['zb_y'].push_back(evt.zero_bias.y_min + dy)
                        pix_containers['zb_val'].push_back(evt.zero_bias.val[dx + s * dy])

            if evt.byte_block.val:
                pix_n[0] += len(evt.byte_block.x)
                offset = (int)((evt.byte_block.side_length-1)/2)
                ival = 0
                xy_set = set([])
                for ixy, x in enumerate(evt.byte_block.x):
                    y = evt.byte_block.y[ixy]
                    for dy in xrange(evt.byte_block.side_length):
                        iy = y + dy - offset
                        if iy < 0 or iy >= xb.res_y: continue
                        for dx in xrange(evt.byte_block.side_length):
                            ix = x + dx - offset
                            if ix < 0 or ix >= xb.res_x: continue

                            # make sure we are not overriding existing pixel
                            if (ix, iy) in xy_set: continue
                            xy_set.add((ix, iy))

                            if not 'bb_x' in pix_containers.iterkeys():
                                pix_containers['bb_x'] = r.vector('UInt_t')()
                                pix_containers['bb_y'] = r.vector('UInt_t')()
                                pix_containers['bb_val'] = r.vector('UInt_t')()
                                events.Branch('bb_x', pix_containers['bb_x'])
                                events.Branch('bb_y', pix_containers['bb_y'])
                                events.Branch('bb_val', pix_containers['bb_val'])
                                
                            pix_containers['bb_x'].push_back(ix)
                            pix_containers['bb_y'].push_back(iy)
                            pix_containers['bb_val'].push_back(evt.byte_block.val[ival])

                            ival += 1

                assert ival == len(evt.byte_block.val)

            for pix in evt.pixels:
                pix_n[0] += 1
                for pix_field,val in pix.ListFields():
                    if not pix_field.name in pix_containers.iterkeys():
                        pix_containers[pix_field.name] = r.vector(type_name[type(val)].full)()
                        events.Branch('pix_'+pix_field.name, pix_containers[pix_field.name])
                    pix_containers[pix_field.name].push_back(val)
                
                
            events.Fill()

    return exposure, events
        

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--in", required=True, dest="infiles", nargs='+', help='.bin files to be converted')
    parser.add_argument("--out", default="pbfiles.root", help='Output file name')
    args = parser.parse_args()
    
    if len(args.infiles) > 1:
        dc = merge_pb(args.infiles)
    else:
        f = open(args.infiles[0])
        dc = pb.DataChunk()
        dc.ParseFromString(f.read())
        f.close()
    outfile = r.TFile(args.out, "recreate")
    exposure, events = pb_to_trees(dc)
    L1pass, L1skip = pb_to_hists(dc)
    outfile.Write()
    raw_input('Press any key to continue')
    outfile.Close()
