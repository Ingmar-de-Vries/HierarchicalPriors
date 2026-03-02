######################################################################################################
# Display head coordinate system (HCS) displacement with respect to device coordinate system (DCS) 
# and with respect to a reference coordinate system (RCS); the latter is specified in a given fiff 
# file, or automatically computed. Process all fiff files in a folder plotting  as darkest blue the 
# first one, as lightest blue the last one. Order is based on acquisition date/time from fiff file tag.
# Also plot Elekta's optimal suggested head position as a coordinate system (OCS).
# Note: --block option specify blocking or non blocking window (useful for bash wrapper).
#
# USAGE: 	python headPositionHistory.py --folder <folder> [--ref reference_fiff_file --block]
#
# E.G:		python headPositionHistory.py 201510130930 201510130930/19431104NTFR_201510130930_2012036_A15run01.fif
#
# Auhor:	Davide Tabarelli (davide.tabarelli@unitn.it)
#
# Revision:	28/09/2016
# Revision:	21/11/2016	:	Now prints raw filenames in order of acquisition time.
# Revision:	22/11/2016	:	If no ref is given computes the most suitable one (see blow for criterion).
# Revision:	28/11/2016	:	Convention (names) changed: now "subject" is "head"
# Revision:	26/01/2017	:	Python 2 code for fif listing in folder (BUG solve)
#
######################################################################################################


######################################################################################################
# SOME MATH: calculating head -> dev transform from dev -> head transform:
#
# y is point in HCS
# x is point in DCS
# ^(tr) means 'transposed'
# (R,T) is dev -> head transform from tag 222 of Elekta fif file
# dev -> head transform (R,T) => y = R x + T
# goal: calculate head -> dev transform (R', T') => y = R' x + T'
#
# y = R x + T
# y - T = R x
# R^(tr) (y - T) = x
# R^(tr) y - R^(tr) T = x
# therefore:
# R' = R^(tr)
# T' = - R^(tr) T
#
# Note: T' checked on a fiff file with the one computed by /neuro/bin/util/show_fiff 
######################################################################################################


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import mne
import datetime
import os
import sys
import warnings
import path
import glob
import argparse
import transforms3d		# For rotation -> euler angles conversion (see http://matthew-brett.github.io/transforms3d/). Requires installation with "pip install transforms3d"
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


##########################################################################################
# Let'define 3D arrow custom class (from http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector)
class Arrow3D(FancyArrowPatch):
	def __init__(self, a, b, *args, **kwargs):
		super().__init__((0,0), (0,0), *args, **kwargs)
		x = [a[0],b[0]];
		y = [a[1],b[1]];
		z = [a[2],b[2]];
		self._verts3d = x, y, z

	def do_3d_projection(self, renderer=None):
		x3d, y3d, z3d = self._verts3d
		x, y, z = proj3d.proj_transform(x3d, y3d, z3d, self.axes.M)
		self.set_positions((x[0],y[0]),(x[1],y[1]))       
        
		return np.min(z)

##########################################################################################

# Temporary disable all warnings
warnings.filterwarnings('ignore');

# Option handler
parser = argparse.ArgumentParser();
parser.add_argument('--folder', help='Fif files folder');
parser.add_argument('--ref', help='Reference run fif file');
parser.add_argument('--block', help='Whether to block execution or not after plot display', action='store_true', default=False);
args = parser.parse_args();
folder = 'SUB02';#args.folder;
blockFlag = args.block;
if args.ref == None:
	refFlag = False;
else:
	refFlag = True;
	refPath = args.ref;

# List of fiff files in given folder
if (sys.version_info > (3, 0)):
	# Python 3 code in this block
	fifs = path.Path(folder).files('*.fif');
else:
	# Python 2 code in this block
	fifs = glob.glob(folder + '*.fif');

Nfifs = len(fifs);

# Constants
L = 25;										# Axes length (mm)
B = [-30, 50, -30, 50, -100, 30];			# Plot boxes/axes limits (mm)
dcsC = [0,0,0];								# DCS axis color (rgb)
ocsC = [0,0.5,0];							# Optimal HCS axis color (rgb)
rcsC = [0.8,0,0];							# Reference HCS axis color (rgb)
hcsC = np.ones([Nfifs,3]);					# HCS axis colors (rgb)
hcsC[...,0] = np.linspace(0,0.75, Nfifs);
hcsC[...,1] = np.linspace(0,0.75, Nfifs);

# Prepare axes layout
fig = plt.figure(figsize=[12,8], facecolor='white');
gs = gridspec.GridSpec(5,4);
ax1 = plt.subplot(gs[0:5,0:3], projection='3d');
ax2 = plt.subplot(gs[0,3]);
ax3 = plt.subplot(gs[1,3], sharex=ax2);
ax4 = plt.subplot(gs[2,3], sharex=ax2);
ax5 = plt.subplot(gs[3,3], sharex=ax2);
ax6 = plt.subplot(gs[4,3], sharex=ax2);

# Prepare 3D coordinate systems view
ax1.set_xlim(B[0], B[1]);
ax1.set_ylim(B[2], B[3]);
ax1.set_zlim(B[4], B[5]);
# ax1.set_aspect('equal');
ax1.view_init(azim=45, elev=15);

# Points defining DCS (Device Coordinate System)
pd = np.array([0, 0, 0]);
xd = np.array([L, 0, 0]);
yd = np.array([0, L, 0]);
zd = np.array([0, 0, L]);

# Plot DCS
ax1.scatter(pd[0], pd[1], pd[2], c=dcsC, marker='o', s=15, edgecolors=dcsC);
x = Arrow3D(pd, xd, mutation_scale=15, lw=1, arrowstyle="->", color=dcsC);
y = Arrow3D(pd, yd, mutation_scale=15, lw=1, arrowstyle="->", color=dcsC);
z = Arrow3D(pd, zd, mutation_scale=15, lw=1, arrowstyle="->", color=dcsC);
ax1.add_artist(x);
ax1.add_artist(y);
ax1.add_artist(z);
ax1.text3D(xd[0]+3, xd[1], xd[2], 'x');
ax1.text3D(yd[0], yd[1], yd[2], 'y (NAS)');
ax1.text3D(zd[0], zd[1], zd[2], 'z');

# Points defining OCS (Elekta's optimal Subject Coordinate System)
po = np.array([0, 0, -40]);
xo = np.array([L, 0, -40]);
yo = np.array([0, L, -40]);
zo = np.array([0, 0, L-40]);

# Plot OCS
ax1.scatter(po[0], po[1], po[2], c=ocsC, marker='o', s=15, edgecolors=ocsC);
x = Arrow3D(po, xo, mutation_scale=15, lw=1, arrowstyle="->", color=ocsC);
y = Arrow3D(po, yo, mutation_scale=15, lw=1, arrowstyle="->", color=ocsC);
z = Arrow3D(po, zo, mutation_scale=15, lw=1, arrowstyle="->", color=ocsC);
ax1.add_artist(x);
ax1.add_artist(y);
ax1.add_artist(z);


# Open all fifs and preload some informations (like recording date and time etc...)
rawH = [];
acqTimeH = np.array([]);
for i in range(0, Nfifs):
	
	rawH.append( mne.io.read_raw_fif(fifs[i], preload=False) );	
	d = rawH[i].info['meas_date'];
	print(d.hour*60+d.minute);									# Open all files
	acqTimeH = np.append( acqTimeH, datetime.datetime.fromtimestamp(d.hour*60+d.minute) );	# Get all recording datetimes (for sorting)


# Compute fif files HCS plot order according to acquisition time
hcsOrder = acqTimeH.argsort();

print(hcsOrder)

# Plot all HCS in folder
R = np.zeros([Nfifs, 3, 3]);
T = np.zeros([Nfifs, 3]);
ph = np.zeros([Nfifs, 3]);
xh = np.zeros([Nfifs, 3]);
yh = np.zeros([Nfifs, 3]);
zh = np.zeros([Nfifs, 3]);
for i in range(0, Nfifs):
	
	# Defined order
	j = hcsOrder[i];
	
	# Get transform HCS -> DCS
	transform = rawH[j].info['dev_head_t'].get('trans');
	Rdh = transform[0:3,0:3];
	Tdh = transform[0:3,3] * 1000.0;	# Moltiplico per 1000 perche' e' in metri
	
	# From dev->head to head->dev
	R[i,:,:] = Rdh.transpose();
	T[i,:] = - Rdh.transpose().dot(Tdh);
	
	# Add Subject Coordinate System (HCS)
	ph[i,:] = R[i].dot(pd) + T[i];
	xh[i,:] = R[i].dot(xd) + T[i];
	yh[i,:] = R[i].dot(yd) + T[i];
	zh[i,:] = R[i].dot(zd) + T[i];
	ax1.scatter(ph[i,0], ph[i,1], ph[i,2], c=hcsC[j], marker='o', s=15, edgecolors=hcsC[j]);
	x = Arrow3D(ph[i], xh[i], mutation_scale=15, lw=1, arrowstyle="->", color=hcsC[j]);
	y = Arrow3D(ph[i], yh[i], mutation_scale=15, lw=1, arrowstyle="->", color=hcsC[j]);
	z = Arrow3D(ph[i], zh[i], mutation_scale=15, lw=1, arrowstyle="->", color=hcsC[j]);
	ax1.add_artist(x);
	ax1.add_artist(y);
	ax1.add_artist(z);



# Plot refecence HCS
if refFlag:
	# If defined plot reference HCS use it ...
	
	# Open the given file
	rawREF = mne.io.read_raw_fif(refPath, preload=False);
	
else:
	# ... otherwise calculate it as follows.
	
	# Among all HCS origin points find the one whose minimize the sum of distances from all the other HCS origins.
	# This HCS is suggested as the most suitable as a reference for MaxMove.
	# Consider that for positive defined quantities minimizing all quantities simultaneously on a N-dim space or minimizing the sum is the same.
	# The problem is the same as finding the geometric median of set of given points (en.wikipedia.org/wiki/Geometric_median), but constraining the 
	# median point to be itself one of the given points. For this reason is computationally more easy. It principle it can be found as the neareast 
	# neighbour of the true geometric median (wasting a lot of computational power ....)
	sum_dist = np.zeros(Nfifs);
	for i in hcsOrder:
		for j in hcsOrder:
			sum_dist[i] += np.linalg.norm(ph[i]-ph[j]);	

	rawREF = rawH[hcsOrder[np.argmin(sum_dist)]];
	
	
# Get transform refSCS -> DCS
transform = rawREF.info['dev_head_t'].get('trans');
refRdh = transform[0:3,0:3];
refTdh = transform[0:3,3] * 1000.0;	# Moltiplico per 1000 perche' e' in metri
# From dev->head to head->dev
refR = refRdh.transpose();
refT = - refRdh.transpose().dot(refTdh);

# Add Reference Coordinate System (refSCS)
pr = refR.dot(pd) + refT;
xr = refR.dot(xd) + refT;
yr = refR.dot(yd) + refT;
zr = refR.dot(zd) + refT;
ax1.scatter(pr[0], pr[1], pr[2], c=rcsC, marker='o', s=15, edgecolors=rcsC);
xR = Arrow3D(pr, xr, mutation_scale=15, lw=1, arrowstyle="->", color=rcsC);
yR = Arrow3D(pr, yr, mutation_scale=15, lw=1, arrowstyle="->", color=rcsC);
zR = Arrow3D(pr, zr, mutation_scale=15, lw=1, arrowstyle="->", color=rcsC);
ax1.add_artist(xR);
ax1.add_artist(yR);
ax1.add_artist(zR);

	

# Calculate movement parameters
d1 = np.zeros(Nfifs);
d2 = np.zeros(Nfifs);
a = np.zeros(Nfifs);
b = np.zeros(Nfifs);
c = np.zeros(Nfifs);
n = np.linspace(1,Nfifs,Nfifs);
for j in hcsOrder:
	
	# Distance (L2-norm) between DCS and HCS origin
	d1[j] = np.linalg.norm(pd-ph[j]);

	# Distance (L2-norm) between REF and HCS origin
	d2[j] = np.linalg.norm(pr-ph[j]);
	
	# Rotation angles follow Tait-Bryan convention (see http://www.ce.unipr.it/people/medici/geometry/node155.html and http://matthew-brett.github.io/transforms3d/reference/transforms3d.taitbryan.html)
	#
	# 	a  =>	rotation around z axis (yaw or 'imbardata' as Girolamo Cardano's nautical angles). Corresponds to head rotation.
	# 	b  =>	rotation around y axis (pitch or 'beccheggio' as Girolamo Cardano's nautical angles). Corresponds to head left/right tilting.
	# 	c  =>	rotation around x axis (roll or 'rollio' as Girolamo Cardano's nautical angles). Corresponds to head front/back inclination.
	#
	a[j], b[j], c[j] = np.degrees(transforms3d.taitbryan.mat2euler(R[j]));

# Plot movement parameters
ax2.plot(n, d1, ':');
ax2.scatter(n, d1, c=hcsC, edgecolors=hcsC);
ax2.set_xticklabels([]);
ax2.grid();
ax2.set_ylabel('mm', fontsize=9);
#ax2.set_title('Distance from device coord. system', fontsize=9);
ax2.set_title('Distance from DCS', fontsize=9);
ax2.tick_params(labelsize=9);

ax3.plot(n, d2, ':');
ax3.scatter(n, d2, c=hcsC, edgecolors=hcsC);
ax3.grid();
ax3.set_ylabel('mm', fontsize=9);
#ax3.set_title('Distance from reference coord. system', fontsize=9);
ax3.set_title('Distance from RCS', fontsize=9);
ax3.tick_params(labelsize=9);

ax4.plot(n, a, ':');
ax4.scatter(n, a, c=hcsC, edgecolors=hcsC);
ax4.grid();
ax4.set_ylabel('degrees', fontsize=9);
ax4.set_title('Head rotation (around z axis)', fontsize=9);
ax4.tick_params(labelsize=9);

ax5.plot(n, b, ':');
ax5.scatter(n, b, c=hcsC, edgecolors=hcsC);
ax5.grid();
ax5.set_ylabel('degrees', fontsize=9);
ax5.set_title('Head left/right tilting (around y axis)', fontsize=9);
ax5.tick_params(labelsize=9);

ax6.plot(n, c, ':');
ax6.scatter(n, c, c=hcsC, edgecolors=hcsC);
ax6.grid();
ax6.set_ylabel('degrees', fontsize=9);
ax6.set_title('Head front/back inclination (around x axis)', fontsize=9);
ax6.tick_params(labelsize=9);


# Plot general legend
fig.suptitle('Head movement estimation for *.fif files in folder ' + folder.strip('/').split('/')[-1], x=0.25, y=0.98, fontsize=12);
fig.text(s='Device coordinate system (DCS)',x=0.03, y=0.9, ha='left', fontsize=10, color=dcsC);
fig.text(s='Elekta optimal suggested head position (OCD)',x=0.03, y=0.88, ha='left', fontsize=10, color=ocsC);
fig.text(s='Head coordinate system (HCS)',x=0.03, y=0.86, ha='left', fontsize=10, color=hcsC[0]);
fig.text(s='Reference run coordinate system (RCS)',x=0.03, y=0.84, ha='left', fontsize=10, color=rcsC);


# Print to console raw data filenames according to acquisition time
print('\n\n\n\n')
print('\nFilenames ordered according to acquisition time:\n')
for i in range(0, Nfifs):
	print('\t # ' + repr(i+1) + '\t' + fifs[hcsOrder[i]].name)

if not refFlag:
	print('\n')
	print('\nSuggested (computed) reference for MaxMove realignment: \n')
	print('\t #' + repr(np.argmin(sum_dist)+1) + '\t' + fifs[hcsOrder[np.argmin(sum_dist)]].name)


# Display graph
fig.tight_layout(pad=0.5);
plt.show(block=blockFlag);				# Set blockFlag to False if you run in a ipython window
#fig.subplots_adjust( left=0.02, bottom=0.02, right=0.98, top=0.98, hspace=0.0 );


