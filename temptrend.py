#!/usr/bin/python
#	Compute temperature trends
#@ Created by P.Inacio <pedromiragaia@gmail.com>

# load packages
import numpy as np
import datetime
import grids
import dates
from utils import call
import os
import argparse

# profiling
import time

# functions
def main():

	# Build parser and parse the input arguments
	args = parse_args()

	# defaults
	ldates = dates.date_range(dates.parse_date_range(args.date_range))
	print ldates

	# load data
	print "Loading grids ..."
	tick = time.time()
	ldates = [ghcn_defaults('initial date') + x*days for x in range(30)]
	lgrds = map(ghcn_load, ldates)
	tock = time.time() - tick
	print '	Elapsed time: %f' % (tock,)

	# accumulate all grids 
	print "Accumulating LS system ..."
	tick = time.time()
	sumA = reduce(np.add, map(grd_to_trend_ls, lgrds))
	tock = time.time() - tick
	print '	Elapsed time: %f' % (tock,)

	# solve the ls of equations
	print "Computing trend ..."
	trend, bias = solve_trend_ls(sumA)

	# now create corresponding grids
	aux = lgrds[0]
	gtrend = grids.Grid(xlim=aux.xlim,ylim=aux.ylim,dx=aux.dx,dy=aux.dy,z=np.reshape(trend,(aux.ny, aux.nx)))
	gbias = grids.Grid(xlim=aux.xlim,ylim=aux.ylim,dx=aux.dx,dy=aux.dy,z=np.reshape(bias,(aux.ny, aux.nx)))

	# save to file
	tfile = 'trend_from_%s_to_%s' % (str(lgrds[0].date), str(lgrds[-1].date))
	gtrend.save(tfile+'.pkl')
	gtrend.plot(tfile+'.pdf', sym_colorbar=True)
	print "Saved trend: %s" % (tfile,)

	tbias = 'bias_from_%s_to_%s' % (str(lgrds[0].date), str(lgrds[-1].date))
	gbias.save(tbias+'.pkl')
	gbias.plot(tbias+'.pdf', sym_colorbar=True)

	return
	
# PLOT FIT ON A GRID
# # plot debug

# # Data for plotting
# # import matplotlib
# import matplotlib.pyplot as plt

# t = list(x[0] for x in data[6911])
# s = list(x[1] for x in data[6911])

# mdl = fit[:,6911]
# f = [x*mdl[0] + mdl[1] for x in t]

# fig, ax = plt.subplots()
# ax.plot(t, s)
# ax.plot(t, f)
# # ax.set(xlabel='time (s)', ylabel='voltage (mV)',
# #        title='About as simple as it gets, folks')
# # ax.grid()
# fig.savefig("test.pdf")

def parse_args():
	"""
	returns a parser to take care of the input arguments
	"""

	# Define the input parser
	desc = "computes long term temperatur anomaly trend for the GHNC dataset"
	epilog = """
datarange input argument is of the format:
		   YYYY[MM[DD]][:YYYY[MM[DD]]]
Where the date before the optional ':'' represents the lower bound of
the range and the optional date after the : represents the upper
bound. The optional elements of the date default to the lowest possible
value for the lower bound and to the maximum possible for the upper
one. For example,
	2006    is equivalent to 2006/01/01:2006/12/31
	2006/02 is equivalent to 2006/02/01:2006/02/28

prod input arguments can be any of:
"""

	parser = argparse.ArgumentParser(description=desc, epilog=epilog,
						formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("daterange",
						help="range of dates to make available locally")

	return parser.parse_args()


def grd_to_trend_ls(dgrd):
	'''
	Given a list set of temperature observations y at time indexes t,
	I want to fit the linear model y = m.t + b to the data where,
		- m is the trend I would like to estimate
		- b is a linear offset constant
	
	The simplest approach solution is to run
		p = np.polyfit(t,y,1)
	on each grid note. This works well, but it is not scalable to handle
	very large datasets distributed in various computers.
	
	Therefore, I will implement a method which allows m to be estimated 
	by combining computations of multiple partitions of the dataset. 
	Given the model,
		y = m.t + b
	We can write in least squares form as,
		y = A x
	where,
		A = [[t_0 1],[t_1 1],...,[t_N,1]]
		y = [y_0, y_1, ..., y_N]'
		x = [m, b]'
	
	The LS solution to this problem is obtained as,
		A' y = A' A x => Y = N x
	
	Notice that if we now get another timeseries of data, we can simply accumulate
	the Y and N matrices, and in fact, the solution for the whole dataset
	can be obtained from solving
		sum(Yi) = sum(Ni) x
	for each i partition of the total dataset.
	
	And now one can estimate x by solving this linear system of equations.
	This is relevant because if we partition the dataset we can compute Y and N 
	for each partition and then accumulate them to compute a global estimate. 
	This approach is scalable because one can have data distributed over various
	workers. 
	
	A little algebra shows that,
		N = [ [t**2, t], [t, 1]]
		Y = [ t*y, y ]
	
	Generalizing an input grid at time t with (N,M) nodes, I will create an
	accumulation matrix A shaped as (6,NxM) where each column contains,
		A = [ N_11 N_12 N_21 N_22 Y_11 Y_12 ]
	
	This function outputs matrix A for the input grid
	'''

	# defaults
	seconds_in_year = 365*24*3600

	# elapsed time since beggining dataset
	ghcn_initial_date = ghcn_defaults('initial date')
	t = (dgrd.date - ghcn_initial_date).total_seconds()/seconds_in_year
	
	# grid size
	n = dgrd.nx*dgrd.ny

	# flatten the temperature 
	Y = dgrd.z.flatten()

	# vetorize the time index to build matrix
	T = np.tile(t, (1, n))

	# build accumulation matrix
	B = np.vstack((T**2, T, T, np.ones((1,n)), T*Y, Y))

	# build bolean index
	I = np.tile(np.isfinite(Y),(6, 1))

	# return zero where no data exists
	A = np.where(I,B,0)

	return A

def solve_trend_ls(sumA):
	'''
	Now take the accumulated A matrix, reshape it and solve for each grid 
	node
		Y = N x

	Notice that the third row of A, by coincidence, tells us the number of
	points accumlated for that grid node.
	I will skip nodes with one or less data points.
	'''

	# loop the grid nodes in the columns
	fit = np.zeros((2, sumA.shape[1]))
	for c in range(sumA.shape[1]):

		# compute trend only if there is more than one data point 
		# accumulated
		if sumA[3,c] > 1:
			# print 'c=', c
			col = sumA[:,c].copy()
			# print 'col=', col
			N = np.reshape(col[0:4],(2,2))
			# print 'N=', N
			Y = col[4:]
			# print 'Y=', Y
			X = np.linalg.solve(N,Y)
			# print 'X=', X
			fit[:,c] = X

	trend = fit[0,:]
	bias = fit[1,:]

	return trend, bias

def ghcn_defaults(key):
	'''
	return paths for different ghcn dataset files
	'''

	if key == 'initial date':
		return datetime.date(1950,01,01)
	elif key == 'final date':
		# TODO: check this
		return datetime.date(2020,06,01)
	elif key == 'xlim':
		return [-180.0, 176.25]
	elif key == 'ylim':
		return [-90, 90.0]
	else:
		raise ValueError('wrong key')

def ghcn_path(key, date=None):
	'''
	return paths for different ghcn dataset files
	'''

	if key == 'url':
		return 'ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/grid/years/%d.tmax' % date.year
	elif key == 'base':
		return os.path.abspath(os.path.dirname(__file__))
	elif key == 'yearly dir':
		return os.path.join(ghcn_path('base'),'data')
	elif key == 'yearly file':
		aux = 'tmp_%d.tmax' % (date.year,)
		return os.path.join(ghcn_path('yearly dir'), aux)
	elif key == 'daily dir':
		aux = '%d' % (date.year,)
		return os.path.join(ghcn_path('base'),'data','daily',aux)
	elif key == 'daily file':
		aux = '%d%02d%02d.pbz2' % (date.year,date.month,date.day)
		return os.path.join(ghcn_path('daily dir', date), aux)
	else:
		raise ValueError('wrong key')

def ghcn_load(date):
	'''
	return a tmax grid for the specified date
	'''

	# check if the daily file is there
	dfile = ghcn_path('daily file',date)

	# if file does not exist need to preprocess this year
	if not os.path.isfile(dfile):
		ghcn_preprocess(date)

	# now load the daily file
	return grids.load(dfile)


def ghcn_preprocess(date):
	'''
	read the file and pickle all daily grids, return the requested grid.
	'''

	# check if the daily file is there
	yfile = ghcn_path('yearly file',date)

	# make sure the file exists
	if not os.path.isfile(yfile):
		ghcn_preload(date)

	# read the file and create daily grids
	# 1st column: Month
	# 2nd column: Day
	# 3rd column: Grid box ID (value range: 1 to 7004, grid spacing = 3.75 deg x 2.5 deg)
	# 4th column: Longitude of lower left corner of grid box (degrees)
	# 5th column: Latitude of lower left corner of grid box (degrees)
	# 6th column: Temperature anomaly (whole degrees Celsius) 
	raw = np.loadtxt(yfile)
	# print 'DEBUG: only reading partial file'
	# raw = np.loadtxt(yfile, max_rows=20000)

	# loop the raw data and create a grid for each day
	rmonth,rday = (int(z) for z in raw[0,0:2])
	istart = 0
	for i, x in enumerate(raw):

		# find index at which day changes
		if rday != x[1]:

			# make new grid
			rdate = datetime.date(date.year,rmonth,rday)
			ghnc_pickle_grid(raw[istart:i,3:],rdate)

			# new day
			istart = i
			rmonth,rday = (int(z) for z in x[0:2])

	# save the last day of the year
	rdate = datetime.date(date.year,rmonth,rday)
	ghnc_pickle_grid(raw[istart:,3:],rdate)

	# post
	ghnc_postload(date)

def ghnc_pickle_grid(raw, rdate):
	'''
	Create a Grid from data raw with date rdate and store it in 
	the data folder 
	'''

	dgrd = grids.Grid(xyz=raw, date=rdate,
			xlim=ghcn_defaults('xlim'),ylim=ghcn_defaults('ylim'))

	# create directory for daily files it not existent
	ddir = ghcn_path('daily dir', rdate)
	if not os.path.isdir(ddir):
		os.makedirs(ddir)

	dgrd.save(ghcn_path('daily file', rdate))
	print 'INFO: stored grid: '+str(rdate)

def ghcn_preload(date):
	'''
	download the yearly ghcn grid file and return a path to that file
	'''

	# check data range
	ghcn_in = ghcn_defaults('initial date')
	ghcn_fn = ghcn_defaults('final date')
	if date < ghcn_in or date > ghcn_fn:
		raise ValueError('year must be between', ghcn_in ,'and ', ghcn_fn)

	# create subdirectories if not available
	ydir = ghcn_path('yearly dir', date)
	if not os.path.isdir(ydir):
		os.makedirs(ydir)

	# download yearly file if not available
	rfile = ghcn_path('url', date)
	lfile = ghcn_path('yearly file', date)

	if not os.path.isfile(lfile):
		call('wget -O'+' '.join([lfile, rfile]), live=True)

def ghnc_postload(date):
	'''
	delete the downloaded text file as we do not need it anymore
	'''

	os.remove(ghcn_path('yearly file', date))

# this idiom means the below code only runs when executed from command line
if __name__ == '__main__':
	main()
