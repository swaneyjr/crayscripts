import ROOT as r
import math

def compare_theta(t, n_bins=100, cut=''):

  hist = r.TH1F('hist','hist',n_bins,-2,2)
  r_hist = r.TH1F('r_hist','r_hist',n_bins,-2,2)
  std_hist = r.TH1F('std_hist','std_hist',n_bins,-2,2)
  
  for track in t:
    hist.Fill(track.theta)
    if track.theta<=0:
      r_hist.Fill(track.theta+math.pi/2)
    else:
      r_hist.Fill(track.theta-math.pi/2)
  
  for b in xrange(n_bins):
    std_hist.SetBinContent(b, math.sqrt(hist.GetBinContent(b)+r_hist.GetBinContent(b)))
                           
  return (hist - r_hist).Divide(std_hist)
    

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
  
  
  
  
