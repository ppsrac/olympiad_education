import numpy as np
import matplotlib.pyplot as plt

f_hp=open("hip catalogue.txt")
lines1=f_hp.readlines()
HIP, mag, RA, Dec=[], [], [], []

for i, line in enumerate(lines1):
    if i>2:
        a=line.split('\t')
        HIP.append(int(a[0]))
        try:
            mag.append(float(a[1]))
        except:
            mag.append(float('nan'))
        try:
            RA.append(float(a[2]))
        except:
            RA.append(float('nan'))
        try:
            Dec.append(float(a[3][:-2]))
        except:
            Dec.append(float('nan'))
f_hp.close()