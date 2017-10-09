import csv 
from bisect import bisect_left

def toSeconds(hms, hms0):
	h, m, s = hms
	return 3600*(float(h)-float(h0)) + 60*(float(m)-float(m0)) + float(s) - float(s0)

with open("C:\Users\ishma\Dropbox (SteinLab)\spectra\(ID     186nm) 2016-07-26-tip15\ivpt\ivpt1.tsv",'r') as tsv:
    lines = [line.strip().split('\t') for line in tsv]
    
(i, vl1, vl2, p, y0, m0, d0, h0, m0, s0) = lines[0]
h0, m0, s0 = float(h0), float(m0), float(s0)

times = [toSeconds((h,m,s),(h0,m0,s0)) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
i = [float(i) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
vl1 = [float(vl1) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
vl2 = [float(vl2) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]
p =[float(p) for (i, vl1, vl2, p, y, m, d, h, m, s) in lines]

i = bisect_left(times, 5000)
print times[i-1]