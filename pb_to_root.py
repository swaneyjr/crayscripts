import bgThresh_pb2 as pb
import ROOT as r
import numpy as np

def merge_pb(infiles):
    dc = pb.DataChunk()
    for fname in infiles:
        f = open(fname)
        dc.MergeFromString(f.read())
        f.close()
    merge_name = fname[0].split('_')[0] + '.bin'
    merge_file = open(merge_name, "wb")
    merge_file.write(dc.SerializeToString())
    merge_file.close()
    return merge_name


def pb_to_tree(fname):
    dc = pb.DataChunk()
    f = open(fname)
    dc.ParseFromString(fname)
    f.close()
    
    # write to ROOT
    exposure = r.TTree('exposure', 'Exposure from protobuf data')
    events = r.TTree('events', 'Events from protobuf data')

    str_type = {int:'/I', long:'/L', float:'/D', bool:'/O', str:'/C'}
    full_type = {int:'int', long:'long', float:'double', bool:'bool',str:'char'}

    # assign appropriate branches to trees

    xb_containers = {}
    evt_containers = {}
    pix_containers = {}
    pix_n = np.zeros(1, dtype=int)

    for xb_field,val in dc.exposure[0].ListFields():
        if xb_field != 'events' and xb_field != 'daq_state':
            xb_containers[xb_field]= np.zeros(1,dtype=type(val))
            exposure.TBranch(xb_field, xb_containers[xb_field], xb_field+str_type[val])
        

    for evt_field,val in dc.exposure[0].events[0].ListFields():
        if evt_field != 'pixels':
            evt_containers[evt_field] = np.zeros(1,dtype=type(val))
            events.TBranch(evt_field, evt_containers[evt_field], evt_field+str_type[val])

    for pix_field,val in dc.exposure[0].events[0].pixels[0].ListFields():
        pix_containers[pix_field] = r.vector(full_type[type(val)])()
        events.TBranch('pix_'+pix_field, pix_containers[pix_field])

    events.TBranch('pix_n', pix_n, 'pix_n/I')
    
    
    # fill tree
    for xb in dc.exposure:
        
        for xb_field,val in xb.ListFields():
            if xb_field == 'daq_state' or xb_field == 'events': continue
            xb_containers[xb_field][0] = val
        exposure.Fill()
        
        for evt in xb.events:
            pix_n[0] = len(evt.pixels)
            for evt_field,val in evt.ListFields()):
                if evt_field == 'pixels': continue
                evt_containers[evt_field][0] = val
                
            # reset pixel vectors
            for v in pix_containers.itervalues():
                v.clear()
                
            for pix in evt.pixels:
                for pix_field,val in pix.ListFields():
                    evt_containers[pix_field].push_back(val)
                
                
            events.Fill()

    return exposure, events
        

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--in", required=True, dest="infiles", nargs='+', help='.bin files to be converted')
    parser.add_argument("--out", default="pbfiles.root", help='Output file name')
    
    if len(args.infiles) > 1:
        fname = merge_pb(args.infiles)
    else:
        fname = args.infiles[0]
    outfile = r.TFile(args.out, "recreate")
    exposure, events = pb_to_events(fname)
    outfile.Write()
    outfile.Close()
