#!/usr/python3

import numpy as np
import math
import progressbar
from astropy.time import Time
import time
import pickle
from matplotlib import pyplot as plt

GPS_week = 1808
GPS_PRN_Num = 32 
GPS_I_L1 = np.zeros((8640000,GPS_PRN_Num),dtype=int);
GPS_Q_L1 = np.zeros((8640000,GPS_PRN_Num),dtype=int);
#GPS_I_L2 = np.zeros((8640000,GPS_PRN_Num),dtype=int);
#GPS_Q_L2 = np.zeros((8640000,GPS_PRN_Num),dtype=int);
#GPS_I_L5 = np.zeros((8640000,GPS_PRN_Num),dtype=int);
#GPS_Q_L5 = np.zeros((8640000,GPS_PRN_Num),dtype=int);

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
        Sat_num = SVID
        Sat_num = 'Unknown'
    return Sat_num, Sat_con

def determineST(SNum):
    if SNum == 0 or SNum == 6 or SNum == 8 or SNum == 24:
        Sig_type = 'L1CA'
    elif SNum == 1 or SNum == 9:
        Sig_type = 'L1P'
    elif SNum == 2 or SNum == 10:
        Sig_type = 'L2P'
    elif SNum == 3 or SNum == 7:
        Sig_type = 'L2C'
    elif SNum == 4 or SNum == 15 or SNum == 25 or SNum == 26:
        Sig_type = 'L5'
    elif SNum == 5 or SNum == 32:
        Sig_type = 'L1C'
    elif SNum == 11:
        Sig_type = 'L2CA'
    elif SNum == 12:
        Sig_type = 'L3'
    elif SNum == 13:
        Sig_type = 'B1C'
    elif SNum == 14:
        Sig_type = 'B2a'
    elif SNum == 17:
        Sig_type = 'E1'
    elif SNum == 19:
        Sig_type = 'E6'
    elif SNum == 20:
        Sig_type = 'E5a'
    elif SNum == 21:
        Sig_type = 'E5b'
    elif SNum == 22:
        Sig_type = 'E5 AltBOC'
    elif SNum == 27:
        Sig_type = 'L6'
    elif SNum == 28:
        Sig_type = 'B1I'
    elif SNum == 29:
        Sig_type = 'B2I'
    elif SNum == 30:
        Sig_type = 'B3I'
    elif SNum == 33:
        Sig_type = 'L1S'
    elif SNum == 34:
        Sig_type = 'B2b'
    elif SNum == 38:
        Sig_type = 'L1CB'
    elif SNum == 39:
        Sig_type = 'L5S'
    else:
        Sig_type = 'Unknown'
    return Sig_type

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

def Nonzero_runs(a):
    # create an arrary that is 1 where a is 0, and pad each end with an extra 0
    isnonzero = np.concatenate(([0], np.not_equal(a,0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(isnonzero))
    # Runs start and end where absdiff is 1
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

#### Calculate SI_raw by NBP and WBP        
def ComputeSIraw(GPS_I,GPS_Q,GPS_SIraw):
    for i in range(GPS_PRN_Num):
        arc = Nonzero_runs(GPS_I[:,i])
        for j in range(len(arc)):
            if arc[j,1] - arc[j,0] >= 2:
                GPS_I_pre = GPS_I[arc[j,0]:arc[j,1]-2,i]
                GPS_I_las = GPS_I[arc[j,0]+1:arc[j,1]-1,i]
                GPS_Q_pre = GPS_Q[arc[j,0]:arc[j,1]-2,i]
                GPS_Q_las = GPS_Q[arc[j,0]+1:arc[j,1]-1,i]
                GPS_SIraw[arc[j,0]:arc[j,1]-2,i] = np.square(GPS_I_pre+GPS_I_las) + np.square(GPS_Q_pre+GPS_Q_las) - np.square(GPS_I_pre) - np.square(GPS_I_las) - np.square(GPS_Q_pre) - np.square(GPS_Q_las)


#### Compute the detrended SI_raw, using the mean over 60s 
def ComputeS4(GPS_SIraw,GPS_S4):
    for i in range(GPS_PRN_Num):
        arc = Nonzero_runs(GPS_SIraw[:,i])
        for j in range(len(arc)):
            if arc[j,1] - arc[j,0] >= 6000:
                GPS_SIraw[arc[j,0]:arc[j,1]-1,i] = np.convolve(GPS_SIraw[arc[j,0]:arc[j,1]-1,i],np.ones(6000)/6000,mode='same')
    
    for i in range(GPS_PRN_Num):
        MinCount = -1 
        for j in range(0,8640000,6000):
            MinCount += 1
            if GPS_SIraw[j,i] != 0:
                if np.all(GPS_SIraw[j:j+5999,i]):
                    SI_sq_ave = np.average(np.square(GPS_SIraw[j:j+5999,i]))
                    SI_ave_sq = math.pow(np.average(GPS_SIraw[j:j+5999,i]),2)
                    GPS_S4[MinCount,i] = math.sqrt((SI_sq_ave-SI_ave_sq)/SI_ave_sq)


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
        lineNum += 1
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
        SNum = int(currentline[2])
        Carrier_cur = float(currentline[3])
        I_cur = int(currentline[4])
        Q_cur = int(currentline[5])
        Sat_num,Sat_con = determineSV(SVID_cur)
        if Sat_con == 'GPS':
            bar.update(lineNum)
#            if lineNum >= 0.01 * lineCount: # Output some data for test
#                break
            SigType = determineST(SNum); 
            if SigType == 'L1CA':
                GPS_I_L1[timeIndex,Sat_num-1] = I_cur
                GPS_Q_L1[timeIndex,Sat_num-1] = Q_cur
#            elif SigType == 'L2C':
#                GPS_I_L2[timeIndex,Sat_num-1] = I_cur
#                GPS_Q_L2[timeIndex,Sat_num-1] = Q_cur
#            elif SigType == 'L5':
#                GPS_I_L5[timeIndex,Sat_num-1] = I_cur
#                GPS_Q_L5[timeIndex,Sat_num-1] = Q_cur
#            else:
#                print('Not the one we want')            
    bar.finish()
    np.savez_compressed('GPS_IQ_L1.npz',GPS_I_L1 = GPS_I_L1, GPS_Q_L1 = GPS_Q_L1)
            
# Load saved GPS L1 IQ data
obj = np.load('GPS_IQ_L1.npz')
namelist = obj.zip.namelist()
obj.zip.extract(namelist[0])
GPS_I_L1 = np.load(namelist[0],mmap_mode='r+')
obj.zip.extract(namelist[1])
GPS_Q_L1 = np.load(namelist[1],mmap_mode='r+')

# Compute S4
GPS_SIraw_L1 = np.zeros((8640000,GPS_PRN_Num),dtype=int)
ComputeSIraw(GPS_I_L1,GPS_Q_L1,GPS_SIraw_L1); 
del(GPS_I_L1)
del(GPS_Q_L1)
GPS_S4_L1 = np.zeros((1440,GPS_PRN_Num),dtype=float)
ComputeS4(GPS_SIraw_L1,GPS_S4_L1)
del(GPS_SIraw_L1)

# Svae S4
np.savez_compressed('GPS_S4_L1.npz',GPS_S4_L1 = GPS_S4_L1)

# Plot the S4 for a specific satellite
plt.title('S4 index test')
plt.plot(range(0,1440), GPS_S4_L1[:,6], label='PRN7')
plt.legend()
plt.show()
