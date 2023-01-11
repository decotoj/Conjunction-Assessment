import spacetrack.operators as op
import datetime as dt
from collections import defaultdict
from datetime import datetime, timedelta
from os import listdir
from os.path import isfile, join
import numpy as np
from astropy.time import Time
from sgp4.api import Satrec
from astropy.coordinates import CartesianRepresentation
from astropy.coordinates import TEME, EarthLocation, AltAz
from astropy import units as u

# Pull Catalog From Spacetrack
def spaceTrackCatalogPull(dt0, out_file, st, nDays=1):

    epochs = op.inclusive_range(dt0.replace(tzinfo=None), dt0.replace(tzinfo=None) + dt.timedelta(days=nDays))
    lines = st.tle_publish(publish_epoch=epochs)

    out = []
    for line in lines:
        out.append(str(line['TLE_LINE1']) + '\n')
        out.append(str(line['TLE_LINE2']) + '\n')

    with open(out_file, 'w') as f:
        f.writelines(out)

    return out

# Load TLE Cat
def parseTLEFile(tleFile):
    tles = defaultdict(list)
    with open(tleFile, 'r') as f:
        lntle = f.readlines()
    for i in range(0,len(lntle),2):
        k = int(float(lntle[i][2:7]))

        tles[k] = {'name': k, 1: lntle[i], 2: lntle[i+1]}
    return tles

# Pull A Data Range of SpaceTrack TLE Catalog & Organize By Date
def pullRangeCatalogs(dt0, dt1, st, catalog_path):

    # Determine TLE Files Existing
    files = [q for q in [join(catalog_path,f) for f in listdir(catalog_path) if isfile(join(catalog_path, f))] if '.txt' in q]

    # Determine Files Missing
    labels = []
    epochs = []
    while dt0 < dt1 and dt0<datetime.utcnow():
        
        b = dt0.strftime('%Y%m%d')
        b = f'{catalog_path}{b}_TLE.txt'
        if b not in files:
            labels.append(b)
            epochs.append(dt0)
        dt0 += timedelta(days=1)

    # Pull Data and Build Missing Files from SpaceTrack
    for i in range(len(labels)):
        print('pulling', labels[i])
        _ = spaceTrackCatalogPull(epochs[i], labels[i], st, nDays=1)
    
    return

# Propagate Entire TLE Catalog Over Range of Epochs and Return TEME Ephemeris as Numpy Array w/ Dimensions RSO->Epoch->Position
def buildTLEEphem(file, epochs):

    # Parse TLE File
    tles= parseTLEFile(file)

    # Format Epochs for SGP4
    epochsAP = Time(epochs)
    jd1 = [t.jd1 for t in epochsAP]
    jd2 = [t.jd2 for t in epochsAP]

    # SGP4 Prep RSOs    
    satellite = [Satrec.twoline2rv(tles[scc][1], tles[scc][2]) for scc in tles.keys()]

    teme_p = np.asarray([rso.sgp4_array(np.array(jd1), np.array(jd2))[1] for rso in satellite])

    return teme_p, tles.keys()

# Given TEME Catalog Ephemeris Return Az/El/Epoch/SCC for RSOs Passes Above Elevation Threshold
def rsoInViewGround(teme, scc, epochs, site, elLimit=10):

    # Get AzEl From Ground Site
    x = teme[:,:,0].flatten()
    y = teme[:,:,1].flatten()
    z = teme[:,:,2].flatten()
    epochsF = np.tile(epochs, teme.shape[0])
    epochsF = epochsF.flatten()
    epochsAP = Time(epochsF)
    temeAP = CartesianRepresentation(x*u.km, y*u.km, z*u.km)
    temeAP = TEME(temeAP, obstime=epochsAP)
    groundsite = EarthLocation(site[0]*u.km, site[1]*u.km, site[2]*u.km)
    elazAP = temeAP.transform_to(AltAz(obstime=epochsAP, location=groundsite))
    
    # Filter Data to Those in View Per Elevation Constraint
    el = elazAP.alt[:].value
    az = elazAP.az[:].value
    scc = np.repeat(np.asarray(list(scc)), teme.shape[1])

    r = np.argwhere(el>elLimit) # rows meeting elevation filter

    return np.squeeze(el[r]), np.squeeze(az[r]), np.squeeze(epochsF[r]), np.squeeze(scc[r])

