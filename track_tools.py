import ROOT as r

def compare_theta(t, bins=100, cut=''):
  hist = r.TH1F('','',bins,-2,2)
  for evt in t:
    hist.Fill(t.hough_theta, )
    hist.Fill(t.hough_theta, )


def track_density(t, region=None, cut=''):
  if region:
    x0,y0,x1,y1 = region
  else:
    x0 = y0 = -float('inf')
    x1 = y1 = float('inf')
  
  total_tracks = t.GetEntries('Sum$(pix_x)/Length$(pix_x) >= x0 && Sum$(pix_x)/Length$(pix_x)' \\
    + '&& Sum$(pix_y)/Length$(pix_y) >= y0 && Sum$(pix_y)/Length$(pix_y) < y1' + cut)
    
  timestamps = [t.
  
  return total_tracks/float(total_frames)
  
  
  
  
