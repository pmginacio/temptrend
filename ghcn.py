import os
import datetime
import grids
import numpy as np
import utils

def defaults(key):
	'''
	return paths for different ghcn dataset files

	note: last time I checked data was available until two days before
	'''
	before_yesterday = datetime.date.today() - datetime.timedelta(days=2)

	if key == 'initial date':
		return datetime.date(1950,01,01)
	elif key == 'final date':
		return before_yesterday
	elif key == 'xlim':
		return [-180.0, 176.25]
	elif key == 'ylim':
		return [-90, 90.0]
	else:
		raise ValueError('wrong key')

def path(key, date=None, ix=None, iy=None):
	'''
	return paths for different ghcn dataset files
	create directories when they are requested
	'''

	if key == 'url':
		return 'ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/grid/years/%d.tmax' % date.year
	elif key == 'base':
		return os.path.abspath(os.path.dirname(__file__))
	elif key == 'data':
		return utils.makedir_if_not_exist(os.path.join(path('base'),'data'))
	elif key == 'yearly file':
		aux = 'tmp_%d.tmax' % (date.year,)
		return os.path.join(path('data'), aux)
	elif key == 'daily dir':
		aux = '%d' % (date.year,)
		ddir = os.path.join(path('base'),'data','daily',aux)
		return utils.makedir_if_not_exist(ddir)
	elif key == 'daily file':
		aux = '%d%02d%02d.pbz2' % (date.year,date.month,date.day)
		return os.path.join(path('daily dir', date), aux)
	elif key == 'computed':
		return utils.makedir_if_not_exist(os.path.join(path('data'),'computed'))
	elif key == 'trend':
		return os.path.join(path('computed'), 
			'trend_from_%s_to_%s.pbz2' % (str(date[0]), str(date[-1])))
	elif key == 'bias':
		return os.path.join(path('computed'),
			'bias_from_%s_to_%s.pbz2' % (str(date[0]), str(date[-1])))
	elif key == 'results':
		return os.path.join(path('base'),'results')
	elif key == 'trend plot':
		return os.path.join(path('results'), 
			'trend_from_%s_to_%s.pdf' % (str(date[0]), str(date[-1])))
	elif key == 'timeseries':
		return os.path.join(path('results'), 
			'trend_from_%s_to_%s_%d_%d.pdf' % (str(date[0]), str(date[-1]),ix,iy))
	else:
		raise ValueError('wrong key: %s' % key)

def load(date):
	'''
	return a tmax grid for the specified date
	'''

	# get the daily grid filename
	dfile = path('daily file',date)

	# if file does not exist need to preprocess this year
	if not os.path.isfile(dfile):
		preprocess(date)

	# now load the daily file
	return grids.load(dfile)


def preprocess(date):
	'''
	read the yearly file and pickle all daily grids
	'''

	# check if the daily file is there
	yfile = path('yearly file',date)

	# make sure the file exists
	if not os.path.isfile(yfile):
		preload(date)

	print "Loading daily grids ..."
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
			pickle_grid(raw[istart:i,3:],rdate)

			# new day
			istart = i
			rmonth,rday = (int(z) for z in x[0:2])

	# save the last day of the year
	rdate = datetime.date(date.year,rmonth,rday)
	pickle_grid(raw[istart:,3:],rdate)

	# post
	postload(date)

def pickle_grid(raw, rdate):
	'''
	Create a Grid from data raw with date rdate and store it in 
	the data folder 
	'''

	dgrd = grids.Grid(xyz=raw, date=rdate,
			xlim=defaults('xlim'),ylim=defaults('ylim'))

	# save
	dgrd.save(path('daily file', rdate))

def preload(date):
	'''
	download the yearly ghcn grid file and return a path to that file
	'''

	# check data range
	ghcn_in = defaults('initial date')
	ghcn_fn = defaults('final date')
	if date < ghcn_in or date > ghcn_fn:
		raise ValueError('year must be between', ghcn_in ,'and ', ghcn_fn)

	# download yearly file if not available
	rfile = path('url', date)
	lfile = path('yearly file', date)

	if not os.path.isfile(lfile):
		print "Downloading %s ..." % rfile
		utils.call('wget -O'+' '.join([lfile, rfile]), live=True)

def postload(date):
	'''
	delete the downloaded text file as we do not need it anymore
	'''

	os.remove(path('yearly file', date))