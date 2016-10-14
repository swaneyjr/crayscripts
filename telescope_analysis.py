#!/usr/bin/env python

import ROOT as r
import numpy as np

FRAME_RATE = 15.
BINS_PER_FRAME = 20

def cleaned_diffs(times):
    diffs = sorted(list(np.diff(np.array(sorted(times)))))
    i = len(diffs) - 1
    while diffs[i] > 1.5 * diffs[i-9]:
        diffs.remove(diffs[i])
        i -= 1
    return diffs

def find_rate(times, max_time = None):
    diffs = cleaned_diffs(times)
    if max_time == None:
        max_time = max(diffs)
    h = r.TH1F("h","diffs", 200, 0, max_time)
    for d in diffs:
        h.Fill(d)
    h.Fit("expo")
    r.gPad.SetLogy()
    return h

def frame_sort(diffs, n_frames):
    frames = np.zeros(2*n_frames+1)
    for i in xrange(len(diffs)):
        frames[i/BINS_PER_FRAME] += diffs[i]
    return frames
    


class TelData(object):
    def __init__(self, m_times = [], l_times = []):
        self.muon_times = m_times
        self.led_times = l_times
    
    # removes muon data when LED is turned off
    def remove_off_time(self):
        t_cutoff = 1000 #gap above which LED is considered off
        led_array = np.array(sorted(self.led_times))
        muon_list = [m for m in self.muon_times if m > led_array[0] \
                      and m < led_array[len(led_array)-1]]
        m = 0
        l = 0
        while m < len(muon_list) and l < len(led_array)-1:
            if muon_list[m] > led_array[l] and muon_list[m] < led_array[l+1]:
                if led_array[l+1] - led_array[l] > t_cutoff:
                    muon_list.remove(muon_list[m])
                else:
                    m += 1
            else:
                l += 1
        self.muon_times = muon_list
        self.led_times = list(led_array)
        return self

    # gets rid hits less than t_shutter after a previous one
    def remove_shutter_effect(self, t_shutter = 0.15):
        # clears closeby muon hits
        muon_array = np.array(sorted(self.muon_times))
        mdiffs = np.diff(muon_array)
        muon_array = muon_array[np.concatenate((mdiffs > t_shutter, np.array([True])))]
        # clears closeby led hits
        led_array = np.array(sorted(self.led_times))
        ldiffs = np.diff(led_array)
        led_array = led_array[np.concatenate((ldiffs > t_shutter, np.array([True])))]
        # clears muons from closeby muon/led pairs
        m = 0
        l = 0
        while m < len(muon_array) and l < len(led_array):
            gap = led_array[l] - muon_array[m]
            if abs(gap) < t_shutter:
                muon_array[m] = 0
                l += 1
                m += 1
            elif gap < 0:
                l += 1
            else:
                m += 1
        
        muon_array = muon_array[muon_array > 0]
        return TelData(list(muon_array),list(led_array))
        
    def append(self, tel_data):
        self.muon_times += tel_data.muon_times
        self.led_times += tel_data.led_times
        return self

    # appends cleaned .root event trees to TelData
    def add_file(self, muon_name, led_name, t_shutter = 0):
        mFile = r.TFile(muon_name)
        lFile = r.TFile(led_name)
        mevents = mFile.Get("events")
        levents = lFile.Get("events")
        file_data = TelData([m.nano*1e-9 for m in mevents], [l.nano*1e-9 for l in levents])
        file_data.remove_off_time()
        if t_shutter > 0:
            file_data.remove_shutter_effect(t_shutter)
        self.append(file_data)

        return self

    def time_slice(self, t_i, t_f = float('inf')):
        return TelData([m for m in self.muon_times if m > t_i and m < t_f], \
                       [l for l in self.led_times if l > t_i and l < t_f])
    
    # subtracts time when experiment was off and returns the remaining time gap
    def t_tot(self):
        return sum(cleaned_diffs(self.led_times))

    # draws a histogram from a list of times and fits to exponential decay
    def muon_rate(self, max_time = None):
        return find_rate(self.muon_times, max_time)

    def led_rate(self, max_time = None):
        return find_rate(self.led_times, max_time)

    # returns number of times when a muon hit and LED hit were within t_coinc
    def find_matches(self, t_coinc, led_delay):
        muon_array = np.array(sorted(self.muon_times))
        led_array = np.array(sorted(self.led_times))
        matches = 0
        m = 0
        l = 0
        while m < len(muon_array) and l < len(led_array):
            gap = led_array[l] - muon_array[m] - led_delay
            if abs(gap) < abs(t_coinc) and np.sign(gap) == np.sign(t_coinc):
                l += 1
                m += 1
                matches += 1
            elif gap < 0:
                l += 1
            else:
                m += 1
        return matches

    # returns a histogram of matches exactly led_delay +/- t_coinc apart
    def shifted_diffs(self, t_max = 0.25, led_delay = 0.0, n_pts = 501):
        if n_pts % 2 == 0:
            n_pts += 1
        t_step = 2* t_max / float(n_pts - 1)
        t_coinc_array = np.linspace(-t_max + 0.5*t_step, t_max - 0.5*t_step, n_pts - 1)
        matches_array = np.zeros(n_pts)
        print "Calculating matches..."
        for i in xrange(n_pts):
            matches_array[i] = self.find_matches(-t_max + i*t_step, led_delay)
            if i % 50 == 0:
                print " %.1f%% (%d/%d)" %(100.*i/n_pts ,i, n_pts)
        print "Filling histogram..."
        diffs = abs(np.diff(matches_array))
        t_array = t_coinc_array + led_delay

        h = r.TH1F("h", "Matches by time interval", n_pts - 1, -t_max + led_delay, t_max + led_delay)
        for n in xrange(n_pts - 1):
            h.Fill(t_array[n], diffs[n])
        r.gROOT.FindObject("c1").Update()
        h.Draw()
        return h
        

    # returns a histogram of matches 
    def match_diffs(self, t_max = 1., n_bins = 250, n_frames = None, minus_bg = False):
        if n_frames != None:
            n_bins = int(BINS_PER_FRAME*(n_frames+0.5))
            t_max = (n_frames + 0.5)/FRAME_RATE
        n_pts = 2 * n_bins + 1
        t_step = 2* t_max / float(n_pts - 1)
        matches_array = np.zeros(n_pts)

        # create array of matches between t_i and t_(i+1)
        print "Calculating matches..."
        for i in xrange(n_pts):
            matches_array[i] = self.find_matches(-t_max + i*t_step, 0.0)
            if i % 50 == 0:
                print " %.1f%% (%d/%d)" %(100.*i/n_pts, i, n_pts)
        diffs = abs(np.diff(matches_array))
        
        # sort into relevant bins
        if n_frames != None:
            print "Sorting into frames..."
            x = np.arange(-n_frames, n_frames+1)
            y = frame_sort(diffs, n_frames)
        else:
            x = np.linspace(-t_max + 0.5*t_step, t_max - 0.5*t_step, n_pts - 1)
            y = diffs
            
        if minus_bg:
            x = -x[0:(len(x)+1)/2]
            for j in xrange((len(y)+1)/2):
                y[j] = abs(y[len(y)-j-1]) - abs(y[j])
            
        print "Filling histogram..."
        h = r.TH1F("h", "Matches", len(x), min(x), max(x)+1)
        for i in xrange(len(x)):
            h.Fill(x[i], y[i])
        return h
