#!/usr/bin/env python
#	utility to compute temperature anomaly trends from the GHCN dataset
# 	run "temptrend -h" for instruction on how to use this utility
#
# Created by P.Inacio <pedromiragaia@gmail.com>

# load packages
# plot
import matplotlib
# 	force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
import ghcn
import datetime
import grids
import dates
from utils import call
import os
import argparse
import itertools

# profiling
import time


# functions
def main():

	# Build parser and parse the input arguments
	args = parse_args()

	# parse date range
	ldates = dates.range_date(*dates.parse_date_range(args.daterange))

	# 
	tfile = ghcn.path('trend', ldates)
	bfile = ghcn.path('bias', ldates)

	if ((os.path.isfile(tfile) and os.path.isfile(bfile)) and 
			not args.recompute):

		print "Loading existing results ..."
		trend = grids.load(tfile)
		bias = grids.load(bfile)

	else:

		trend, bias = compute_trend(ldates)

	# plot grid
	tplot = ghcn.path('trend plot',date=ldates)
	trend.plot(tplot, sym_colorbar=True, transpose=True)
	print "Saved plot: %s" % (tplot,)

	# plot timeseries
	if args.timeseries:

		# get lon, lat pair
		lon, lat = args.timeseries

		# load timeseries for a coordinate
		print "Loading timeseries for (%.2f,%.2f) ..." % (lon, lat)
		s = np.array(map(lambda x: x[lon,lat], map(ghcn.load, ldates)))

		# get time idx
		mdate = ldates[len(ldates)/2]
		t = np.array([(x - mdate).days for x in ldates])

		# get fit
		f = trend[lon,lat]/365.0*t + bias[lon,lat]

		# plot
		plot_timeseries(t, s, f, lon, lat, ldates)
	
def parse_args():
	"""
	returns a parser to take care of the input arguments
	"""

	# Define the input parser
	desc = "computes long term temperature anomaly trend for the GHNC dataset"
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
"""

	parser = argparse.ArgumentParser(description=desc, epilog=epilog,
						formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("daterange",
						help="range of dates to make available locally")
	parser.add_argument('-t',"--timeseries",nargs=2,metavar=('lon','lat'),type=float,
						help="plot timeseries for the lon lat pair of coordinates")
	parser.add_argument('-r',"--recompute",default=False,action='store_true',
						help="force recompute trend")

	return parser.parse_args()


def compute_trend(ldates):
	'''
	request trend computation for the input list of dates
	'''

	# accumulate all grids 
	# notice that the reduce(map ...) allows for minimal memory usage
	# and for code parelellization
	tick = time.time()
	sumA = reduce(np.add, map(grd_to_trend_ls, 
		map(ghcn.load, ldates), itertools.repeat(ldates[len(ldates)/2],len(ldates))))
	tock = time.time() - tick
	print '	Elapsed time: %f' % (tock,)

	# solve the ls of equations
	print "Computing trend ..."
	trend, bias = solve_trend_ls(sumA, thrs=len(ldates)/10.0)

	# now create corresponding grids
	aux = ghcn.load(ldates[0])
	trend = grids.Grid(xlim=aux.xlim,ylim=aux.ylim,dx=aux.dx,dy=aux.dy,
		z=np.reshape(trend,(aux.nx, aux.ny)),
		zlabel='temperature anomaly trend',zunit='deg/yr',
		xlabel='longitude',xunit='deg',
		ylabel='latitude',yunit='deg')
	bias = grids.Grid(xlim=aux.xlim,ylim=aux.ylim,dx=aux.dx,dy=aux.dy,
		z=np.reshape(bias,(aux.nx, aux.ny)))

	if not os.path.isdir(ghcn.path('results')):
		os.makedirs(ghcn.path('results'))

	# save computed results
	trend.save(ghcn.path('trend',ldates))
	bias.save(ghcn.path('bias',ldates))

	return trend, bias	

def grd_to_trend_ls(dgrd, mdate):
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

	# elapsed time since middle of the dataset
	t = (dgrd.date - mdate).total_seconds()/seconds_in_year
	
	# grid size
	n = dgrd.nx*dgrd.ny

	# 
	Z = dgrd.z.reshape(n,1)
	T = np.tile(t,(n,1))

	A = np.stack((
			np.hstack((T**2, T)),
			np.hstack((T, np.ones((n,1)))),
			np.hstack((Z*T, T))),axis=2)

	# set invalid grid nodes to 0
	I = np.isnan(Z).flatten()
	A[I,:,:] = 0

	return A

def solve_trend_ls(sumA, thrs=1):
	'''
	Now take the accumulated A matrix, reshape it and solve for each grid 
	node
		Y = N x

	Notice that the third row of A, by coincidence, tells us the number of
	points accumlated for that grid node.
	I will skip nodes with one or less data points.
	'''

	# exclude nodes with less than thrs number of points
	I = (sumA[:,1,1] > thrs).flatten()

	N = sumA[:,:,0:2]
	Y = sumA[:,:,2]

	X = np.full((len(sumA),2),np.nan)
	X[I,:] = np.linalg.solve(N[I,:,:],Y[I,:])

	trend, bias = np.hsplit(X,2)

	return trend, bias

def plot_timeseries(t, s, f, lon, lat, ldates):
	'''
	create a plot for the timeseries s and the fitter linear model t
	'''
	
	fig, ax = plt.subplots()
	ax.plot(t, s)
	ax.plot(t, f)
	ax.set(xlabel='# days', ylabel='temperature anomaly [deg]',
	        title='Temperature anomaly for (%.2f,%.2f) \nbetween %s and %s' % (lon, lat, 
	        	str(ldates[0]), str(ldates[-1])))
	ax.grid()
	tsplot = ghcn.path('timeseries',date=ldates,ix=0,iy=0)
	fig.savefig(tsplot)
	print "Saved plot: %s" % (tsplot,)

# this idiom means the below code only runs when executed from command line
if __name__ == '__main__':
	main()
