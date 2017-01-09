import bgThresh_pb2 as pb
import ROOT as r
import numpy as np

def merge_pb(infiles):
    dc = pb.DataChunk()
    for fname in infiles:
        f = open(fname)
        dc.MergeFromString(f.read())
        f.close()
    merge_name = infiles[0].split('_')[-2] + '.bin'
    print "Compiling messages into %s" % merge_name
    merge_file = open(merge_name, "wb")
    merge_file.write(dc.SerializeToString())
    merge_file.close()
    print "Messages saved."
    return merge_name


def pb_to_trees(fname):
    dc = pb.DataChunk()
    f = open(fname)
    dc.ParseFromString(f.read())
    f.close()
    
    # write to ROOT
    exposure = r.TTree('exposure', 'Exposure from protobuf data')
    events = r.TTree('events', 'Events from protobuf data')

    short_type = {int:'/i', long:'/l',float:'/D',bool:'/O',unicode:'/C'}
    full_type = {int:'int', long:'long',float:'double',bool:'bool',unicode:'char'}

    # assign appropriate branches to trees

    xb_containers = {}
    evt_containers = {}
    pix_containers = {}
    pix_n = np.zeros(1, dtype=int)

    for xb_field,val in dc.exposure_blocks[0].ListFields():
        if xb_field.label != 3: # save repeated fields for events
            xb_containers[xb_field.name]= np.zeros(1,dtype=object)
            exposure.Branch(xb_field.name, xb_containers[xb_field.name], \
                            xb_field.name+short_type[type(val)])
        

    for evt_field,val in dc.exposure_blocks[0].events[0].ListFields():
        if evt_field.label != 3: # save repeated fields for pixels
            evt_containers[evt_field.name] = np.zeros(1,dtype=object)
            events.Branch(evt_field.name, evt_containers[evt_field.name], \
                          evt_field.name+short_type[type(val)])

    for pix_field,val in dc.exposure_blocks[0].events[0].pixels[0].ListFields():
        pix_containers[pix_field.name] = r.vector(full_type[type(val)])()
        events.Branch('pix_'+pix_field.name, pix_containers[pix_field.name])

    events.Branch('pix_n', pix_n, 'pix_n/I')
    
    
    # fill tree
    print "Filling tree"
    
    total_xb = len(dc.exposure_blocks)
    for ixb,xb in enumerate(dc.exposure_blocks):

        if ixb % (total_xb/10) == 0:
            print " %d/%d (%0.2f)" % (ixb, total_xb, 100.*ixb/total_xb)
        
        for xb_field,val in xb.ListFields():
            if xb_field.label == 3: continue
            xb_containers[xb_field.name][0] = val
            if type(val) == unicode:
                xb_containers[xb_field.name][0] += '\0'
        exposure.Fill()
        
        for evt in xb.events:
            pix_n[0] = len(evt.pixels)
            for evt_field,val in evt.ListFields():
                if evt_field.label == 3: continue
                evt_containers[evt_field.name][0] = val
                
            # reset pixel vectors
            for v in pix_containers.itervalues():
                v.clear()
                
            for pix in evt.pixels:
                for pix_field,val in pix.ListFields():
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
        fname = merge_pb(args.infiles)
    else:
        fname = args.infiles[0]
    outfile = r.TFile(args.out, "recreate")
    exposure, events = pb_to_trees(fname)
    outfile.Write()
    outfile.Close()
