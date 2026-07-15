import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from stoemer import e_dipole, stormer_approximation
from numba import prange
from random import uniform as rand
import typing as T
import argparse


R = 6371  # Earth radius
ISS_r = 400 + R  # ISS radius
incl = np.radians(51.64)

# This program generates a number of plots for the TIGERISS Pioneers proposal.
#       3D Plot of ISS orbits
#       ISS overlay on Stoermer Rigidity
#       ISS plot of # orbits

def hav(ang):
    return np.sin(ang/2)**2

def orientation(t0, t1):  # Function to generate pointing of satelitte based on current motion.
    # HEading equation taken from https://www.movable-type.co.uk/scripts/latlong.html - under topic Bearing
    thetainit, phiinit = t0
    thetaf, phif = t1
    dtheta = np.radians(thetaf-thetainit)
    if phiinit > phif:
        direction, direct_i = -1, "South"
    else:
        direction, direct_i = 1, "North"
    # haver = hav(np.radians(phif-phiinit))+np.cos(np.radians(phiinit))*np.cos(np.radians(phif))*hav(np.radians(thetaf-thetainit))
    # c = 2 * np.arctan2(np.sqrt(haver),np.sqrt(1-haver))
    # d = 2* ISS_r * c
    heading = np.arctan2(np.cos(np.radians(phiinit))*np.sin(np.radians(phif))
                         - np.sin(np.radians(phiinit))*np.cos(np.radians(phif))* np.cos(dtheta),
                          np.sin(dtheta)*np.cos(np.radians(phif)))
    pointing = np.degrees(heading)
    return pointing





def ISS_pointrange(points_wanted=5001, orbits_wanted=1):  # Quick function to create time range of points
    prat = points_wanted  # points per orbit
    orbits = orbits_wanted  # number of orbits
    tmaxpoint = orbits * 2 * np.pi  # set max value for time in units of radians
    timespan = np.linspace(0,tmaxpoint,int(orbits*prat))
    #timespan = np.linspace(-3*np.pi/2, np.pi/2, int(orbits*prat))  # generating points in time.
    print("Points:", int(orbits*prat))  # print number of points
    return tmaxpoint, timespan








def orbit_posxyz(timespan, orbittimemin, incl_sat, orbitrad=ISS_r):  # Equation of motion for a satellite: Taken from
    # Goldstein Classical Mechanics 3rd ed. Obviously only works for LEO circular.
    oneorbit = orbittimemin
    times = timespan
    phi = times*oneorbit  # multiplying by orbits per day (1 orbit is 92.68 min)
    xi = times
    x_sat = orbitrad * (np.cos(phi)*np.cos(xi)+np.sin(phi)*np.sin(xi)*np.cos(incl_sat))
    y_sat = orbitrad * (np.sin(phi)*np.cos(xi)-np.cos(phi)*np.sin(xi)*np.cos(incl_sat))
    z_sat = orbitrad * (np.sin(incl_sat)*np.sin(xi))
    return x_sat, y_sat, z_sat




def genISSridityplot(x_pos, y_pos, z_pos):  # Function to generate monte carlo contour plot of stoermer
    # approx with ISS orbit.
    # Generating dataset of stoermer apporixmation values
    c_lat = np.floor(np.degrees(-np.pi/3))
    max_lat = np.ceil(np.degrees(np.pi/3))
    init_lon = 0
    c_lon = init_lon
    avgrig = []
    events = []
    total = 0
    print(c_lat, max_lat)
    edip = e_dipole(2015)
    while c_lat < max_lat:
        while c_lon < 360:
            for _ in prange(2000):
                total += stormer_approximation(rand(c_lon, c_lon+1), rand(c_lat, c_lat+1), 0, 90, 400, edip)  # calculate
                # random lat and lon (in bounds) for stoermer approximation
            if total > 0:
                rigid = total/2000
            else:
                rigid = 0
            print(c_lat, c_lon, rigid)
            avgrig.append(rigid)
            events.append(2000)
            c_lon += 1
            total = 0
        c_lat += 1
        c_lon = init_lon
    y_range = np.floor(np.degrees(-np.pi/3))
    resolution = 1
    x_s = np.arange(0, 360, 1)
    y_s = np.arange(y_range, max_lat, 1)
    #  X_s, Y_s = np.meshgrid(x_s, y_s)
    fig1, ax1 = plt.subplots()
    avgrig = np.reshape(avgrig, [len(y_s), len(x_s)])
    cs = ax1.contour(avgrig, 20, extent=(0, 360-resolution, y_range, -y_range-resolution))
    fig1.set_size_inches(11.6, 7)
    ax1.clabel(cs, inline=1, fontsize=12)
    plt.xlim((0, 360))
    plt.ylim((-60, 60))
    ax1.set_title('Rigidity Plot with ISS overlay')
    ax1.set_ylabel("Latitude", size=20)
    ax1.set_xlabel("(East) Longitude", size=20)
    # Converting point data from 3-d to lat-lon plot
    THETA_pos = np.degrees(np.arctan2(x_pos, y_pos))  # longitude
    PHI_pos = np.degrees(np.arctan(z_pos/np.sqrt(x_pos**2+y_pos**2)))  # latitude
    temp_x = []
    temp_y = []
    a1 = 0  # counting plot loops val
    val1 = 0  # counting position in list
    place1 = 0  # placeholder
    plotter_ck = True
    # Following loop makes sure each 360 degrees of longitude are a different color in the plot.
    for j in range(len(THETA_pos)):
        if np.floor(THETA_pos[j] % 360) == 0 and plotter_ck is True:
            a1 += 1
            place1 = j
            if a1 % 20 == 0:
                print("cycle", a, j)
            ax1.plot(temp_x, temp_y)
            temp_x, temp_y, temp_y_ = [], [], []
            plotter_ck = False
        else:
            val1 += 1
            if val1 > place1 + 10:
                plotter_ck = True
            temp_x.append(THETA_pos[j] % 360)
            temp_y.append(PHI_pos[j])
            if j > len(THETA_pos)-2:
                ax1.plot(temp_x, temp_y, label="Data")
                print("Done:", a1)
    plt.hlines([-np.degrees(incl)-.1, np.degrees(incl)+.1], xmin=0, xmax=360, colors="black", linestyles='dashed')
    # add dashed line for min max latitude
    plt.show()  # Creating merged ISS orbit plot with stoermer rigidities




def geo2mag(lat,lon):
    dlong = np.radians(360-(180-107.4)) # South magnetic pole 1995
    dlat = np.radians(80.4) # south mag
    #dlong = np.radians(288.59)
    #lat = np.radians(79.3)
    Rval = 1# ISS_r/6371
    glat = np.radians(lat)
    glon = np.radians(lon)
    glat_shape = glat.shape
    glon_shape = glon.shape
    coord = [glat,glon,Rval]
    x_pos = Rval * np.cos(coord[0])*np.cos(coord[1])
    y_pos = Rval * np.cos(coord[0])*np.sin(coord[1])
    z_pos = Rval * np.sin(coord[0])
    p_arr = np.array((x_pos,y_pos,z_pos))
    #print(p_arr)
    geolong2maglong = np.zeros((3,3))
    geolong2maglong[0, 0] = np.cos(dlong)
    geolong2maglong[1, 0] = np.sin(dlong)
    geolong2maglong[0, 1] = -np.sin(dlong)
    geolong2maglong[1, 1] = np.cos(dlong)
    geolong2maglong[2, 2] = 1.0
    #print(geolong2maglong)
    #out = geolong2maglong * p_arr
    out = geolong2maglong.T @ np.array([x_pos, y_pos, z_pos])
    #print(out)
    tomaglat = np.zeros((3, 3))
    tomaglat[0, 0] = np.cos(np.pi/ 2 - dlat)
    tomaglat[2, 0] = -np.sin(np.pi/ 2 - dlat)
    tomaglat[0, 2] = np.sin(np.pi / 2 - dlat)
    tomaglat[2, 2] = np.cos(np.pi / 2 - dlat)
    tomaglat[1, 1] = 1.0
    #print(final)
    out = tomaglat.T @ out
    mlat = np.arctan2(out[2], np.sqrt(out[0] ** 2 + out[1] ** 2))
    mlat = np.degrees(mlat)
    mlon = np.arctan2(out[1], out[0])
    mlon = np.degrees(mlon)
    '''mlat = np.arctan2(final[2,2], np.sqrt(final[0,0] ** 2 + final[1,1] ** 2))
    mlat = np.degrees(mlat)
    mlon = np.arctan2(final[1,1], final[0,0])
    mlon = np.degrees(mlon)'''
    return mlat.reshape(glat_shape), mlon.reshape(glon_shape)



def geomagplot(points,res,x_pos, y_pos, z_pos):
    c_lat = np.floor(-90)#np.degrees(-np.pi / 2))
    max_lat = np.ceil(np.degrees(np.pi / 2))
    init_lon = 0
    c_lon = init_lon
    avgrig = []
    events = []
    total_n = 0
    print(c_lat, max_lat)
    while c_lat < max_lat-.0001:
        while c_lon < 360:
            for _ in prange(points):
                tup =  geo2mag(rand(c_lat, c_lat + res),rand(c_lon, c_lon + res)) # calculate
                total_n += tup[0]
                # random lat and lon (in bounds) for stoermer approximation
                rigid = total_n / points
            print(c_lat, c_lon, rigid)
            avgrig.append(rigid)
            events.append(2000)
            c_lon += res
            total_n = 0
        c_lat += res
        c_lon = init_lon
    y_range = np.floor(-90)#np.degrees(-np.pi / 2))
    resolution = res
    x_s = np.arange(0, 360, res)
    y_s = np.arange(y_range, max_lat, res)
    #  X_s, Y_s = np.meshgrid(x_s, y_s)
    fig1, ax1 = plt.subplots()
    avgrig = np.reshape(avgrig, [len(y_s), len(x_s)])
    cs = ax1.contour(avgrig, 25, extent=(0, 360 - resolution, 90, -90))
    fig1.set_size_inches(11.6, 7)
    ax1.clabel(cs, inline=1, fontsize=12)
    plt.xlim((0, 360))
    plt.ylim((-90, 90))
    ax1.set_title('Magnetic Pointing Plot with ISS overlay')
    ax1.set_ylabel("Latitude", size=20)
    ax1.set_xlabel("(East) Longitude", size=20)

    # Converting point data from 3-d to lat-lon plot
    THETA_pos = np.degrees(np.arctan2(x_pos, y_pos))  # longitude
    PHI_pos = np.degrees(np.arctan(z_pos / np.sqrt(x_pos ** 2 + y_pos ** 2)))  # latitude
    temp_x = []
    temp_y = []
    a1 = 0  # counting plot loops val
    val1 = 0  # counting position in list
    place1 = 0  # placeholder
    plotter_ck = True
    # Following loop makes sure each 360 degrees of longitude are a different color in the plot.
    for j in range(len(THETA_pos)):
        if np.floor(THETA_pos[j] % 360) == 0 and plotter_ck is True:
            a1 += 1
            place1 = j
            if a1 % 20 == 0:
                print("cycle", a, j)
            ax1.plot(temp_x, temp_y)
            temp_x, temp_y, temp_y_ = [], [], []
            plotter_ck = False
        else:
            val1 += 1
            if val1 > place1 + 10:
                plotter_ck = True
            temp_x.append(THETA_pos[j] % 360)
            temp_y.append(PHI_pos[j])
            if j > len(THETA_pos) - 2:
                ax1.plot(temp_x, temp_y, label="Data")
                print("Done:", a1)
    plt.hlines([-np.degrees(incl) - .1, np.degrees(incl) + .1], xmin=0, xmax=360, colors="black", linestyles='dashed')
    # add dashed line for min max latitude
    plt.show()  # Creating merged ISS orbit plot with stoermer rigidities
