#!/usr/bin/python
# -*- encoding: utf8
 
# automatically retrieve weather data from the Norwegian weather service
# see: http://eklima.met.no/wsKlima/start/start_en.html


import xml.etree.ElementTree as et
from urllib2 import urlopen

from datetime import datetime as dt, timedelta as td
from copy import copy


AVAIL_TSTYPES = {'hourly': 2, 'daily': 0, 'monthly': 1}
NOBS_MAX = 400000

def _get(kworder, **kwargs):
	''' Basic method for accessing the eklima data portal. 
 All keyword arguments are transformed into http arguments.
 Returns the parsed XML response as a minidoc instance. 
 
 Unfortunately, the eklima API is sensitive to the order of 
 the parameters and expects also empty parameters to be present. 
 Hence the kworder argument. It is the list of (possibly unset) 
 arguments in the order expected by the eklima API. '''
	
	url = 'http://eklima.met.no/metdata/MetDataService'
	urlarg = ''
	for kw in kworder:  
		argv = kwargs.get(kw, '')
		if type(argv) == list or type(argv) == dict:
			argv = '%2C'.join(argv)
		urlarg += '&%s=%s' % (kw, str(argv))
	
	u = urlopen(url+'?'+urlarg[1:])
	s = u.read()
	u.close()

	return et.fromstring(s)


def __strftime_filename(date):
	''' Take a datetime object and return a formatted string 
 in the date format used as part of filenames. '''
	return '%04d%02d%02d' % (date.year, date.month, date.day)


def __strftime_prettyprint(date):
	''' Take a datetime object and return a formatted string 
 in the date format for nicely formatted output. '''
	return '%04d-%02d-%02dT%02d:%02d' % (date.year, date.month, date.day, date.hour, date.minute)


def __strftime_eklima(date):
	''' Take a datetime object and return a formatted string 
 in the date format the eklima data portal expects its it to be. '''
	return '%02d-%02d-%02d' % (date.year, date.month, date.day)


def _get_stations_from_tstype_elems(tstype, elems):
	''' Get stations by given time series type and given 
 element codes '''
	
	order = ['invoke', 'timeserietypeID', 'elem_codes', 'username']
	request = {'invoke': 'getStationsFromTimeserieTypeElemCodes',
		   'timeserietypeID': AVAIL_TSTYPES[tstype], 
		   'elem_codes': elems,
	}
	
	tree = _get(order, **request)
	stnrs  = tree.iter('stnr')
	names  = tree.iter('name')
	lats   = tree.iter('latDec')
	lons   = tree.iter('lonDec')
	elevs  = tree.iter('amsl')
	estds  = tree.iter('fromYear')
	demds  = tree.iter('toYear')
	
	stations = {}
	for stnr,name,lat,lon,elev,estd,demd in zip(stnrs,names,lats,lons,elevs,estds,demds):
		stnr = stnr.text
		stations[stnr] = {'name': name.text, 
				  'lat': float(lat.text), 
				  'lon': float(lon.text), 
				  'elev': int(elev.text),
				  'estd': int(estd.text),
				  'demd': int(demd.text), }

	return stations


def _get_data(tstype, dfrom, dto, stations, elems, hours=[0,6,12,18], months=''):
	''' Get data for given stations and elements. '''
	
	if type(elems) == str:
		elems = [elems,]
	elems = map(lambda elem: elem.lower(), elems)

	dates = []
	data = {}
	reported = set([])
	for elem in elems:
		data[elem] = {}
		data['quality_'+elem] = {}
	

	tdiff = dto - dfrom

	nobs = (tdiff.days+1)*len(hours)*len(stations)*len(elems)
	days_max = NOBS_MAX / (len(hours)*len(stations)*len(elems))
	print 'Requested approximately %d observations' % nobs
	if nobs > NOBS_MAX:
		print 'Max days per request: %d' % days_max

	dstart = dfrom 
	while dstart < dto:
		dend = min(dstart+td(days_max-1,86399), dto)
		print 'Requesting period: %s -- %s' % (str(dstart), str(dend))

		order = ['invoke', 'timeserietypeID', 'format', 'from', 'to', 'stations', 
				'elements', 'hours', 'months', 'username']
		request = {'invoke': 'getMetData',
			   'timeserietypeID': AVAIL_TSTYPES[tstype], 
			   'from': __strftime_eklima(dstart),
			   'to': __strftime_eklima(dend),
			   'stations': stations,
			   'elements': elems, 
			   'hours': hours, 
			   'months': months,
		}
		
		tree = _get(order, **request)


		# Iteration over all time stamps
		obsroot = list(tree.iter('timeStamp'))[0]
		for obs in obsroot.findall('item'):
			date = dt.strptime(obs.find('from').text, '%Y-%m-%dT%H:%M:%S.000Z')
			dates.append(date)
			for elem in elems:
				data[elem][date] = {}
				data['quality_'+elem][date] = {}
			
			remain_st = stations.keys()
			# Iteration over all weather stations
			for st in obs.find('location').findall('item'):
				stnr = st.find('id').text
				reported.add(stnr)
				remain_st.pop(remain_st.index(stnr))

				remain_elems = copy(elems)
				# Iteration over all elements
				for elem in st.find('weatherElement').findall('item'):
					elemcode = elem.find('id').text.lower()
					value    = float(elem.find('value').text)
					quality  = int(elem.find('quality').text)

					data[elemcode][date][stnr] = value
					data['quality_'+elemcode][date][stnr] = quality
					remain_elems.pop(remain_elems.index(elemcode))
				
				# Iteration over non-observed elements
				for elemcode in remain_elems:
					data[elemcode][date][stnr] = -99999.0
					data['quality_'+elemcode][date][stnr] = -99999.0

			# Iteration over all weather stations that did not report at this time
			for stnr in remain_st:
				for elemcode in elems:
					data[elemcode][date][stnr] = -99999.0
					data['quality_'+elemcode][date][stnr] = -99999.0
		
		# Remove stations that never reported
		for stnr in stations.keys():
			if stnr not in reported:
				stations.pop(stnr)
				for elemcode in elems:
					for date in dates:
						del data[elemcode][date][stnr]
						del data['quality_'+elemcode][date][stnr]

		dstart = dstart + td(days_max)

	return stations, dates, data


def _pretty_print(filename, stations, dates=[], data=[]):
	''' Pretty print the data for ONE element of the data 
 array returned by _get_data() '''

	f = file(filename, 'w')
	f.write('# Eklima data downloaded by eklima.py\n')
	f.write('#\n')
	f.write('# Download time: %s\n' % __strftime_prettyprint(dt.now()))
 	if len(dates) > 0:
		first = dates[0]
		last  = dates[-1]
		f.write('# First date:    %s\n' % __strftime_prettyprint(first))
		f.write('# Last date:     %s\n' % __strftime_prettyprint(last))
		stnrs = data[first].keys()
	else:
		stnrs = stations.keys()
	f.write('#\n')
	f.write('# Station data:\n')
	f.write('# %s\t%32s\t%s\t%s\t%s\t%s\t%s\n' % ('no', 'name', 'lat', 'lon', 'elev', 'since', 'until') )
	for stnr in stnrs:
		demd = str(stations[stnr]['demd'])
		if demd == '0': demd = ''
		f.write((u'# %s\t%32s\t%.4f\t%.4f\t%d\t%d\t%s\n' % (stnr, stations[stnr]['name'], 
			  stations[stnr]['lat'], stations[stnr]['lon'], stations[stnr]['elev'], 
			  stations[stnr]['estd'], demd) 
			).encode('utf8')  )
	if len(dates) > 0:
		f.write('#\n')
		f.write('# Weather data:\n')
		f.write('# Datetime      \t%9s\n' % '\t'.join(map(lambda x: '%9s' % x, data[first].keys())) )
		for date in dates: 
#			f.write('%s\t%s\n' % (__strftime_prettyprint(date), 
#				'\t'.join(map(lambda x: '%9.1f' % x, data
			f.write('%s ; \t%s\n' % (__strftime_prettyprint(date), 
				'\t'.join(map(lambda x: '%9.1f' % x, data
[date].values())) )  )
	f.close()

	return


def _write_mat(filename, stations, dates, data, elem):
	import numpy as np
	import scipy.io.matlab as mio

	matdates = map(lambda date: np.array([date.year, date.month, date.day, 
			date.hour, date.minute, date.second], dtype='f8'), dates)
	
	a = {'dates': matdates}			

	print elem 

	for stnr in stations:

		stnum = 'st%s' % stnr

		a[stnum]={}
		a[stnum][elem]=[]
		a[stnum]['quality_'+elem]=[]
		

		a[stnum]['name']=[]
		a[stnum]['name'].append(stations[stnr]['name'])
		a[stnum]['lat']=[]
		a[stnum]['lat'].append(stations[stnr]['lat'])
		a[stnum]['lon']=[]
		a[stnum]['lon'].append(stations[stnr]['lon'])
		a[stnum]['elev']=[]
		a[stnum]['elev'].append(stations[stnr]['elev'])
		a[stnum]['estd']=[]
		a[stnum]['estd'].append(stations[stnr]['estd'])
		a[stnum]['demd']=[]
		a[stnum]['demd'].append(stations[stnr]['demd'])
			
			#data['quality_'+elem]
			
		for date in dates:
			a[stnum][elem].append(data[elem][date][stnr])
			a[stnum]['quality_'+elem].append(data['quality_'+elem][date][stnr])	
				
	mio.savemat(filename,a)

	return


def download_hourly_by_elems(elems, dfrom, dto, additional_elems=[], hours='', months='', stations=[], 
		verbose=True, savetext=False, savemat=True):
	''' Download time series of elems and additional elems 
 for all stations that provide hourly measurements of any elem in elems. '''

	if not savetext and not savemat:
		print 'WARNING: both savetext and savemat is turned off!'

	# Find matching stations by elems
	if len(stations) == 0:
		stations = _get_stations_from_tstype_elems('hourly', elems)
	if verbose:
		print 'Requesting %d stations' % len(stations)
	
	if not type(elems) == list:
		elems = [elems, ]
	if not type(additional_elems) == list:
		additional_elems = [additional_elems, ]
	elems.extend(additional_elems)
	elems = map(lambda x: x.lower(), elems)

	if hours == '':
		hours = map(str, range(24))
	elif type(hours) == list:
		hours = map(str, hours)

	if type(months) == list:
		months = map(str, months)
	
	# Retrieve data for matching stations and for both elems and additional_elems
	stations, dates, data = _get_data('hourly', dfrom, dto, stations, elems, hours, months)
	if verbose:
		print 'Retrieved %d time steps from %d stations' % (len(dates), len(stations))
	
	# Save the data
	if savetext:
		for elem in data:
			_pretty_print('%s_%s-%s.txt' % (elem, __strftime_filename(dfrom), 
					__strftime_filename(dto)), stations, dates, data[elem])
	if savemat:
			for elem in elems:
				#print elems[0]
				_write_mat('%s_%s-%s.mat' % (elem, __strftime_filename(dfrom), __strftime_filename(dto)),
					stations, dates, data, elem)
		#else:
		#	_write_mat('eklima_%s-%s.mat' % (__strftime_filename(dfrom), __strftime_filename(dto)), 
		#			stations, dates, data)
	#print len(elems)
	
	return stations


# the end
