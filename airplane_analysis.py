#!/usr/bin/env python

import ROOT as r
import numpy as np

def find_runtime(times, cutoff=300000):
    time_diffs = np.diff(np.sort(times))
    return np.sum(time_diffs[time_diffs < cutoff])/60000. # in minutes

def get_hist(values, bins=None):
    if not bins:
        bins = np.arange(min(values), max(values)+1)
    hist,bins = np.histogram(values, bins)
    return hist

def get_pix_val_rate(t):
    runtime = find_runtime([evt.timestamp for evt in t])
    hist = get_hist([max(evt.pix_val) for evt in t])
    return hist/runtime, runtime

def get_tree(filename):
    f = r.TFile(filename)
    return f.Get('events')


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = '')
    parser.add_argument('--air', required='true', help='Name of TFile for airplane events')
    parser.add_argument('--ground', help='Name of TFile for sea level events')
    args = parser.parse_args()

    air = get_tree(args.air)
    ground = get_tree(args.ground)

    air_rate, airtime = get_pix_val_rate(air)
    ground_rate, groundtime = get_pix_val_rate(ground)

    air_val = r.TH1F("air", "Airplane", 256, 0, 256)
    ground_val = r.TH1F("ground", "Sea level", 256, 0, 256)

    air_val.GetXaxis().SetTitle('pix_val')
    ground_val.getXaxis().SetTitle('pix_val')

    for evt in air:
        air_val.Fill(max(air.pix_val), 1./airtime)

    for evt in ground:
        ground_val.Fill(max(ground.pix_val), 1./groundtime)

    c = r.TCanvas()

    air_val.Draw()
    ground_val.Draw('same')

    raw_input('Press enter to continue')


    if ground:
        import matplotlib.pyplot as plt

        plt.figure(1)
        plt.title('Ratio of air to ground rates')
        plt.xlabel('pix_val')

        bins = np.arange(256)[ground_rate != 0]
        air_nz = air_rate[ground_rate != 0]
        ground_nz = ground_rate[ground_rate != 0]

        plt.errorbar(bins, air_nz/ground_nz, yerr=air_nz/ground_nz*np.sqrt(1./np.sqrt(air_nz) + 1./np.sqrt(ground_nz)))
    
        cutoff = raw_input('Select a signal cutoff')

        print "Sea level rate: {}".format(ground_rate[cutoff:])
        print "Airplane rate: {}".format(air_rate[cutoff:])

