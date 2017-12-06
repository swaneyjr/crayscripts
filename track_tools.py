#!/usr/bin/env python

import ROOT as r
import math

def compare_theta(t, name='hist', n_bins=100, cut=None):

  hist = r.TH1F(name,'Track angle excess, compared to 90 degree rotation',n_bins,-math.pi/2,math.pi/2)
  hist.GetXaxis().SetTitle('Track angle')
  hist.GetYaxis().SetTitle('Standard deviations above perpendicular')
  r_hist = r.TH1F('r_hist','r_hist',n_bins,-math.pi/2,math.pi/2)
  std_hist = r.TH1F('std_hist','std_hist',n_bins,-math.pi/2,math.pi/2)
  
  if cut:
    r.gROOT.cd()
    t = t.CopyTree(cut)
  for track in t:
    hist.Fill(track.theta)
    if track.theta<0:
      r_hist.Fill(track.theta+math.pi/2)
    else:
      r_hist.Fill(track.theta-math.pi/2)
  
  for b in xrange(n_bins):
    std_hist.SetBinContent(b, math.sqrt(hist.GetBinContent(b)+r_hist.GetBinContent(b)))
                           
  hist.Add(-1*r_hist)
  hist.Divide(std_hist)
  return hist
    

def track_density(frames, tracks, region=None, cut=''):
  if region:
    x0,y0,x1,y1 = region
  else:
    x0 = y0 = -float('inf')
    x1 = y1 = float('inf')
  
  total_tracks = tracks.GetEntries('Sum$(pix_x)/Length$(pix_x) >= x0 && Sum$(pix_x)/Length$(pix_x)' \
    + '&& Sum$(pix_y)/Length$(pix_y) >= y0 && Sum$(pix_y)/Length$(pix_y) < y1' + cut)
  
  total_frames = frames.GetEntries(cut)
  
  return total_tracks/float(total_frames)
  
  
  
  
