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

def get_pix_val_rate(t, bins=None):
    runtime = find_runtime([evt.timestamp for evt in t])
    print "Runtime: %f min" % runtime
    hist = get_hist([max(evt.pix_val) for evt in t], bins)
    hist.resize(256)
    return hist/runtime, runtime

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = '')
    parser.add_argument('--air', required='true', help='Name of TFile for airplane events')
    parser.add_argument('--ground', required='true', help='Name of TFile for sea level events')
    args = parser.parse_args()

    fair = r.TFile(args.air)
    fground = r.TFile(args.ground)

    air = fair.Get('events')
    ground = fground.Get('events')

    print "ROOT files successfully loaded"

    bins = range(60)+range(60,90,2)+range(90,150,3)+range(150,256,5)

    air_rate, airtime = get_pix_val_rate(air, bins)
    ground_rate, groundtime = get_pix_val_rate(ground, bins)

    air_val = r.TH1F("airhist", "Airplane", 256, 0, 256)
    ground_val = r.TH1F("groundhist", "Sea level", 256, 0, 256)

    air_val.GetXaxis().SetTitle('pix_val')
    ground_val.GetXaxis().SetTitle('pix_val')

    air_val.SetStats(False)
    ground_val.SetStats(False)
    
    air_val.SetLineColor(2)
    ground_val.SetLineColor(4)
    print "Histograms created"

    for evt in air:
        for val in evt.pix_val:
            air_val.Fill(val, 1./airtime)

    for evt in ground:
        for val in evt.pix_val:
            ground_val.Fill(val, 1./groundtime)

    ground_val.SetMaximum(10)
    ground_val.Draw()
    air_val.Draw('same')

    r.gPad.SetLogy()
    r.gPad.BuildLegend()

    import matplotlib.pyplot as plt

    fig = plt.figure()
    plt.ylim(1,100)
    plt.yscale('log')
    plt.title('Ratio of air to ground rates')
    plt.xlabel('pix_val')

    nonzero = np.logical_and(ground_rate != 0, air_rate != 0)

    xbins = np.arange(256)[nonzero]
    air_nz = air_rate[nonzero]
    ground_nz = ground_rate[nonzero]

    ratios = air_nz/ground_nz
    errorbars = ratios*np.sqrt(1./(air_nz*airtime) + 1./(ground_nz*groundtime))

    plt.errorbar(xbins, air_nz/ground_nz, yerr=errorbars)
    plt.show()

    cutoff = int(raw_input('Select a signal cutoff: '))

    sea_rate_tot = np.sum(ground_rate[cutoff:])
    sea_rate_err = (sea_rate_tot/groundtime)**0.5
    air_rate_tot = np.sum(air_rate[cutoff:])
    air_rate_err = (air_rate_tot/airtime)**0.5

    print "Sea level rate: {} +/- {} hit/min".format(sea_rate_tot, sea_rate_err)
    print "Airplane rate: {} +/- {} hit/min".format(air_rate_tot, air_rate_err)
    print "Average multiplier: {} +/- {}".format(air_rate_tot/sea_rate_tot, air_rate_tot/sea_rate_tot*((sea_rate_err/sea_rate_tot)**2 + (air_rate_err/air_rate_tot)**2)**0.5)

    fair.Close()
    fground.Close()

