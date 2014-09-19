import matplotlib.pyplot as plt
import numpy as np
file  = "../../data/msjidlmh_140906_1057_50.csv"
x = np.loadtxt(file,delimiter=',',dtype=np.int32)
N = 500000
t = x[:,1]
dt =t[1:] - t[0:-1]
bw = (dt.max()-dt.min())/(N-1.)
bin1 = dt.min()+ bw*np.arange(N)
#define the array to hold the occurrence count
bincount= np.array([])
for bin in bin1:
    count = np.where((dt>=bin)&(dt<bin+bw))[0].size
    bincount = np.append(bincount,count)
fig = plt.figure()
ax = fig.add_subplot(111)
binc = bin1+0.5*bw
plt.plot(binc,bincount,drawstyle= 'steps-mid')
plt.title("Histogram", fontsize=12)
plt.xlabel("Interval[ticks]",fontsize=12)
plt.ylabel("Frequency",fontsize=12)
plt.xlim(0,4000)
plt.ylim(0,700)
ax.annotate('afterpulse', xy=(250, 250), xytext=(300, 300),arrowprops=dict(facecolor='black', shrink=0.05),)
ax.annotate('original pulse', xy=(100, 500), xytext=(300, 500),arrowprops=dict(facecolor='black', shrink=0.05),)
plt.title("Example of an afterpulse in on tick cycle",fontsize=12)
plt.show()