"""

	histW.py: 

		To make a histogram of the weights. 

"""

import numpy as np; 
import pylab as plt;

fIn = open("./dataSync.dat",'r');
dR = fIn.read();
fIn.close();
dL = dR.split('\n');
dL = [ll for ll in dL if ll!=''];

dl = [];
dPW = [];
for ll in dL: 
	dl += [int(ll.split("  ")[0])];
	dPW += [float(ll.split("  ")[1])];

print sum(dPW);


plt.figure();
plt.plot(dl,dPW);
plt.show();

