#!/usr/bin/env python2

import ROOT as r
import numpy as np

def find_runtime(times, cutoff=300000):
    time_diffs = np.diff(np.sort(times))
    return np.sum(time_diffs[time_diffs < cutoff])/60000. # in minutes

def get_hist(values, bins=None):
    if not bins:
        bins = np.arange(256)
    hist, bins = np.histogram(values, bins)
    return hist

def get_pix_val_rate(t_tuple, bins=None, gps=False):
    runtime = 0
    hist = np.zeros(len(bins)-1) if bins else np.zeros(255)
    for t in t_tuple:
        runtime += find_runtime([evt.gps_fixtime if gps else evt.timestamp for evt in t])
        hist += get_hist([max(evt.pix_val) for evt in t], bins)
    print "Runtime: %f min" % runtime
    return hist/runtime, runtime

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = '')
    parser.add_argument('--air', required='true', nargs='+', help='Name of TFile for airplane events')
    parser.add_argument('--ground', required='true', nargs='+', help='Name of TFile for sea level events')
    parser.add_argument('--expected', type=float, help='Expected ratio of air/ground cosmic rates')
    args = parser.parse_args()

    fair = map(r.TFile, args.air)
    fground = map(r.TFile, args.ground)

    air = map(lambda f: f.Get('events'), fair)
    ground = map(lambda f: f.Get('events'), fground)

    print "ROOT files successfully loaded"

    bins = range(14,40)+range(40,60,2)+range(60,90,3)+range(90,120,5)+range(120,160,20)+[170,256]

    air_rate, airtime = get_pix_val_rate(air, bins)
    ground_rate, groundtime = get_pix_val_rate(ground, bins, gps=True)

    air_val = r.TH1F("airhist", "Airplane", 256, 0, 256)
    ground_val = r.TH1F("groundhist", "Sea level", 256, 0, 256)

    air_val.GetXaxis().SetTitle('pix_val')
    ground_val.GetXaxis().SetTitle('pix_val')

    air_val.SetStats(False)
    ground_val.SetStats(False)
    
    air_val.SetLineColor(2)
    ground_val.SetLineColor(4)
    print "Histograms created"

    for ta in air:
        for evt in ta:
            for val in evt.pix_val:
                air_val.Fill(val, 1./airtime)

    for tg in ground:
        for evt in tg:
            for val in evt.pix_val:
                ground_val.Fill(val, 1./groundtime)

    ground_val.SetMaximum(10)
    ground_val.Draw()
    air_val.Draw('same')

    r.gPad.SetLogy()
    r.gPad.BuildLegend()

    import matplotlib.pyplot as plt

    fig = plt.figure()
    plt.ylim(0,40)
    #plt.yscale('log')
    plt.title('Ratio of air to ground rates')
    plt.xlabel('pix_val')
    plt.xscale('log')

    nonzero = np.logical_and(ground_rate != 0, air_rate != 0)

    xbins = np.array(bins)[:-1][nonzero]
    air_nz = air_rate[nonzero]
    ground_nz = ground_rate[nonzero]

    ratios = air_nz/ground_nz
    errorbars = ratios*np.sqrt(1./(air_nz*airtime) + 1./(ground_nz*groundtime))
    plt.figure(1)
    plt.errorbar(xbins, air_nz/ground_nz, yerr=errorbars, fmt='o')
    plt.xlabel('pix_val')
    plt.title('Ratio of air/ground rates by pix_val')

    if args.expected:
        signal_frac = (ratios - 1) / (args.expected - 1)
        exp_err = 5
        frac_errors = signal_frac * np.sqrt(exp_err**2 / args.expected**2 + errorbars**2 / ratios**2)
        plt.figure(2)
        plt.errorbar(xbins, signal_frac, yerr=frac_errors, fmt='o')
        plt.xlabel('pix_val')
        plt.ylabel('signal / total')
        plt.semilogx()
        plt.title('Signal rate as a fraction of total rate')

    plt.show()

    cutoff = int(raw_input('Select a signal cutoff: '))

    sea_rate_tot = np.sum(ground_rate[cutoff:])
    sea_rate_err = (sea_rate_tot/groundtime)**0.5
    air_rate_tot = np.sum(air_rate[cutoff:])
    air_rate_err = (air_rate_tot/airtime)**0.5

    print "Sea level rate: {} +/- {} hit/min".format(sea_rate_tot, sea_rate_err)
    print "Airplane rate: {} +/- {} hit/min".format(air_rate_tot, air_rate_err)
    print "Average multiplier: {} +/- {}".format(air_rate_tot/sea_rate_tot, air_rate_tot/sea_rate_tot*((sea_rate_err/sea_rate_tot)**2 + (air_rate_err/air_rate_tot)**2)**0.5)

    map(lambda f: f.Close(), fair+fground)

