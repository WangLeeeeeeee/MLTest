#!/usr/python3

import numpy as np
import math
import progressbar
from astropy.time import Time

GPS_week = 1808
GPS_I = np.zeros((8640000,32),dtype=int);
GPS_Q = np.zeros((8640000,32),dtype=int);

def gps2utc(GPS_week, GPS_tow):
    GPS_sec = GPS_week*7*86400 + GPS_tow 
    TimeGPS_cur = Time(GPS_sec,format='gps')
    TimeUTC_cur = Time(TimeGPS_cur,format='iso',scale='utc')
    UTC_hou = int(TimeUTC_cur.value[11:12])
    UTC_min = int(TimeUTC_cur.value[14:15])
    UTC_sec = float(TimeUTC_cur.value[17:22])
    return UTC_hou, UTC_min, UTC_sec

def determineSV(SVID):
    if SVID_cur >= 0 and SVID_cur <= 37: 
        Sat_num = SVID; # GPS
        Sat_con = 'GPS'
    elif SVID_cur >= 38 and SVID_cur <= 61:
        Sat_num = SVID - 37 # GLONASS
        Sat_con = 'GLONASS'
    elif SVID_cur >= 63 and SVID_cur <= 68:
        Sat_num = SVID - 38 # GLONASS
        Sat_con = 'GLONASS'
    elif SVID_cur >= 71 and SVID_cur <= 106:
        Sat_num = SVID - 70 # GALILEO
        Sat_con = 'GALILEO'
    elif SVID_cur >= 120 and SVID_cur <= 140: 
        Sat_num = SVID - 100
        Sat_con = 'SBAS'
    elif SVID_cur >= 141 and SVID_cur <= 180: 
        Sat_num = SVID - 140
        Sat_con = 'BeiDou'
    elif SVID_cur >= 181 and SVID_cur <= 187: 
        Sat_num = SVID - 180
        Sat_con = 'QZSS'
    else:
        print('Unknow SVID: ',SVID)
    return Sat_num, Sat_con

def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

def fileLineCount(filepath):
    with open(filepath,'rb') as fp:
        c_generator = _count_generator(fp.raw.read)
        lineCount = sum(buffer.count(b'\n') for buffer in c_generator)
    return lineCount
        
        
filepath = "/home/wangleu/Documents/share/ismroutput.txt"  
lineCount = fileLineCount(filepath) 
with open(filepath,'r') as filestream:
    # Read IQ data and save them to matrix
    # the rwo is time in 10ms and column is PRN
    bar = progressbar.ProgressBar(maxval=lineCount,\
            widgets=[progressbar.Bar('=','[', ']'),' ',progressbar.Percentage()])
    bar.start()
    lineNum = 0;
    timeStartFlag = 0;
    for line in filestream:
        bar.update(lineNum+1)
        currentline = line.split(",")
        GPS_tow = float(currentline[0])
        if timeStartFlag == 0:
            # Determine the first epoch
            utc_hou,utc_min,utc_sec = gps2utc(GPS_week,GPS_tow) # cost too much time
            UTC_curSec = int((utc_hou*3600 + utc_min*60 + utc_sec)*100)
            timeIndex = UTC_curSec
            timePre = GPS_tow
            timeStartFlag = 1
        else:
            if GPS_tow != timePre:
                timeIndex = timeIndex + int((GPS_tow-timePre)*100)
                timePre = GPS_tow
        #utc_hou,utc_min,utc_sec = gps2utc(GPS_week,GPS_tow) # cost too much time
        #UTC_curSec = int((utc_hou*3600 + utc_min*60 + utc_sec)*100)
        SVID_cur = int(currentline[1])
        Carrier_cur = float(currentline[3])
        I_cur = int(currentline[4])
        Q_cur = int(currentline[5])
        Sat_num,Sat_con = determineSV(SVID_cur)
        if Sat_con == 'GPS':
            GPS_I[timeIndex,Sat_num-1] = I_cur
            GPS_Q[timeIndex,Sat_num-1] = Q_cur
    bar.finish()
            
    # Calculate S4
    # 1. Calculate SI_raw by NBP and WBP
    GPS_SIraw = zeros((8640000,32),dtype=int)
    for i in range(32):
        for j in range(8640000-1):
            if GPS_I[j,i] and GPS_Q[j,i] and GPS_I[j+1,i] and GPS_Q[j+1,i]:
                GPS_SIraw[j,i] = math.pow(GPS_I[j,i]+GPS_I[j+1,i],2) + math.pow(GPS_Q[j,i]+GPS_Q[j+1,i],2) - math.pow(GPS_I[j,i],2) - math.pow(GPS_I[j+1,i],2) - math.pow(GPS_Q[j,i],2) - math.pow(GPS_Q[j+1,i],2)
    # delete the GPS_I and GPS_Q to save memory
    del(GPS_I)
    del(GPS_Q)
    # 2. Compute the detrended SI_raw, using the mean over 60s
    GPS_SIdet = zeros((8640000,32),dtype=int)
    for i in range(32):
        for j in range(8640000-6000-1):
            if np.all(GPS_SIraw[j:j+5999,i]):
                SI_trend = np.average(GPS_SIraw[j:j+5999,i])
                GPS_SIdet[j:j+5999,i] = np.divide(GPS_SIraw[j:j+5999,i],SI_trend)
    # 3. Compute S4 index
    del(GPS_SIraw)
    GPS_S4 = zeros((1440,32),dtype=int)
    for i in range(10):
        while j < 86400:
            if np.all(GPS_SIdet[j:j+5999,i]):
                SI_sq_ave = np.average(np.square(GPS_SIdet[j:j+5999,i]))
                SI_ave_sq = math.pow(np.average(GPS_SIdet[j:j+5999,i]),2)
                GPS_S4[int(j/6000),i] = math.sqrt((SI_sq_ave-SI_ave_sq)/SI_ave_sq)
                j = j + 5999                
            else:
                j += 1



