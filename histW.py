"""

	Copyright (C) 2012, 2013, 2015 Luis F Seoane. 

		Contact: luis.seoane@upf.edu, brigan@gmail.com


	This file is part of NeuralCryptography.

    NeuralCryptography is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NeuralCryptography is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with NeuralCryptography.  If not, see <http://www.gnu.org/licenses/>.

"""

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

