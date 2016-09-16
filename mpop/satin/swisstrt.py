import Image
import logging
import glob
import os
from ConfigParser import ConfigParser
import glob
import numpy as np
import numpy.ma as ma
from mpop import CONFIG_PATH

import pyresample
import datetime
from time import gmtime

import logging
#LOG = logging(__name__)



def load(satscene, *args, **kwargs):
    """Loads the *channels* into the satellite *scene*.
    """

    print "... load TRT data (trt.py)"


    # Dataset information
    #
    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

    values = {"satname": satscene.satname,
              "number": satscene.number,
              "instrument": satscene.instrument_name,
              "satellite": satscene.fullname
              }
    # print values
    projectionName = conf.get("radar-level2", "projection")

    gmt = gmtime()
    now = datetime.datetime(gmt.tm_year, gmt.tm_mon, gmt.tm_mday, gmt.tm_hour, gmt.tm_min, 0) 
    dtime = now - satscene.time_slot
    print "... desired time slot: "+ str(satscene.time_slot)+ ", current time: "+ str(now)
    print "    time difference: "+ str(dtime)
    NEAR_REAL_TIME=False

    archive_time = datetime.timedelta(days=2) # last 2( or 3) days are in the archive
    if dtime < archive_time:
        print "    use the near real time data "
        NEAR_REAL_TIME=True
        filename = os.path.join(
            satscene.time_slot.strftime(conf.get("radar-level2", "nearrealtime_dir")),
            satscene.time_slot.strftime(conf.get("radar-level2", "filename", raw=True)) % values)
    else: 
        print "    use the all time archive "
        data_dir = satscene.time_slot.strftime(conf.get("radar-level2", "archive_dir"))
        filename = os.path.join(
            satscene.time_slot.strftime(conf.get("radar-level2", "archive_dir")),
            satscene.time_slot.strftime(conf.get("radar-level2", "filename", raw=True)) % values)


    # Load data from txt file 
    print "... read data from", str(filename)
    if len(glob.glob(str(filename))) == 0:
        "*** WARNING, no file "+str(filename)+" found!"
        filename=""
    elif len(glob.glob(str(filename))) > 1:
        "*** WARNING, more than one file "+str(filename)+" found!"
        filename = glob.glob(str(filename))[0]
    else: 
        filename = glob.glob(str(filename))[0]

    # get projection 
    satscene.area = pyresample.utils.load_area(os.path.join(CONFIG_PATH, "areas.def"), projectionName)
   
    # Read TRTcells
    print "... read TRT cells from file: " + filename
    satscene.traj_IDs, satscene.TRTcells, satscene['TRTcells'] = readRdt(filename, **kwargs)


# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# class definition of a cell
class TRTcell:
    """TRT cell as class"""
    date=[]     # yyyymmddHHMM
    lon=[]      # 6.3862   [longitude in degree E]
    lat=[]      # 46.8823  [latitude in degree N]
    jCH=[]      # 20       [longitude position in Swiss projection]
    iCH=[]      # 80       [latitude position in Swiss projection]
    ell_L=[]    # 6.74     [ellipse long in km]
    ell_S=[]    # 3.37     [ellipse short in km]
    angle=[]    # 4.9      [ellipse orientation in degree from north clockwise]
    area=[]     # 46       [area of the cell (see det) in km2]
    vel_x=[]    # 1.01     [zonal velocity in km/h]
    vel_y=[]    # -5.49    [meridional velocity in km/h]
    det=[]      # 36       [detection threshold in dBZ (36...48)]
    RANKr=[]    # 3        [RANK = round(2*VIL(n) + 2*ET45m(n) + 1*dBZmax(n) + 2*area55dBZ(n)* / 7) ] (n)== normalized[0:4] in 1/10
    CG_minus=[] # 0        [cloud to ground flashs in the last 5 min with negative polarization]
    CG_plus=[]  # 0        [cloud to ground flashs in the last 5 min with positive polarization]
    CG=[]       # 0        [total cloud to ground flashs in the last 5 min]
    perc_CG_plus=[] # 0    [precentage of positive cloud to ground flashs in the last 5 min]
    ET45 =[]    # 1.8      [maximum height of ECHO TOP 45dBZ in km]
    ET45m=[]    # 1.7      [median height of ECHO TOP 45dBZ in km]
    ET15=[]     # 6.5      [maximum height of ECHO TOP 15dBZ in km]
    ET15m=[]    # 6.1      [median height of ECHO TOP 15dBZ in km]
    VIL=[]      # 2.0      [vertical integrated liquid water in kg/m2]
    maxH=[]     # 3.4      [maximum height of max. radar signal in km]
    maxHm=[]    # 2.2      [median height of max. radar signal in km]
    POH=[]      # 0        [always 0 placeholder at the moment (probability of hail in percent)]
    RANK=[]     # 0        [rank as integer, better use RANKr]
    Dvel_x=[]   # 0.58     [change of zonal velocity in km/h / 15min] 
    Dvel_y=[]   # 3.17     [change of meridional velocity in km/h / 15min]                                                            
    # cell_contour_lon_lat=[] # 6.3680; 46.9506;  6.3684; 46.9326;  6.3554; 46.9235;  6.3558; 46.9055;  6.3429; 46.8964;  6.3438; 46.8514; 
    #                         # [(lon, lat) pairs of the contour]   
    # def f(self):
    #    return 'hello world'


# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# input 
# file           [string]   input file name, e.g. "/srn/data/TRTC/CZC1510916200T.trt"
def readRdt(filename, **kwargs):

    from os import stat
    from numpy import zeros, rot90, floor
    from ApproxSwissProj import ApproxSwissProj 
    import numpy as np

    # array of names of all available cells 
    traj_IDs=[]  # 2015041916100004
    RANKs=[]
    cell_area=[]

    if 'cell' in kwargs.keys():
        print "... read specific cell (rdt.py): ", kwargs['cell']
        cell=kwargs['cell']
    else:
        print "... read all cells"

    if 'min_rank' in kwargs.keys():
        print "... read only cell with minimum rank (rdt.py): ", kwargs['min_rank']
        min_rank=kwargs['min_rank']
    else:
        print "... read all ranks"
        min_rank=int(-1)

    # dictionary of cells
    TRTcells = {}

    ni=640
    nj=710
    # empty masked array
    cell_mask = ma.asarray(zeros(shape=(ni,nj)))

    # read data from file
    if filename != "" and stat(filename).st_size != 0:

        file = open(filename,"r")

        # define class for reprojection to ccs4 grid 
        ccs4 = ApproxSwissProj()

        # read the data from the input file 
        for line in file:
            line2=line.strip()
            if len(line2) > 0:
                if line2.startswith("@"):
                    # print "version: "+ line2
                    pass
                elif line2.startswith("#"):
                    # print "comment: " + line2
                    pass
                else:
                    data = line2.split("; ")
                    # print "data: " + line2
                    if ('cell' not in kwargs.keys()) or ('cell' in kwargs.keys() and data[0] == cell):
                        if ('min_rank' not in kwargs.keys()) or ('min_rank' in kwargs.keys() and int(data[11]) >= min_rank):
                            traj_IDs.append(data[0])                          # awk '{print $1}' CZC1420416000T.trt
                            RANKs.append(int(data[11]))                       
                            cell_area.append(int(data[7]))
                            #print "create TRT cell: " + data[0]
                            TRTcells[data[0]] = TRTcell()
                            TRTcells[data[0]].date         = data[ 1]         # awk '{print $2}' CZC1420416000T.trt
                            TRTcells[data[0]].lon          = float(data[ 2])  # awk '{print $3}' CZC1420416000T.trt
                            TRTcells[data[0]].lat          = float(data[ 3])  # awk '{print $4}' CZC1420416000T.trt
                            TRTcells[data[0]].ell_L        = float(data[ 4])  # awk '{print $5}' CZC1420416000T.trt  
                            TRTcells[data[0]].ell_S        = float(data[ 5])  # awk '{print $6}' CZC1420416000T.trt 
                            TRTcells[data[0]].angle        = float(data[ 6])  # awk '{print $7}' CZC1420416000T.trt 
                            TRTcells[data[0]].area         = int  (data[ 7])  # awk '{print $8}' CZC1420416000T.trt
                            TRTcells[data[0]].vel_x        = float(data[ 8])    
                            TRTcells[data[0]].vel_y        = float(data[ 9]) 
                            TRTcells[data[0]].det          = int  (data[10]) 
                            TRTcells[data[0]].RANKr        = int  (data[11])  # awk '{print $10}' CZC1420416000T.trt
                            TRTcells[data[0]].CG_minus     = int  (data[12])
                            TRTcells[data[0]].CG_plus      = int  (data[13])
                            TRTcells[data[0]].CG           = int  (data[14])
                            TRTcells[data[0]].perc_CG_plus = int  (data[15])
                            TRTcells[data[0]].ET45         = float(data[16])
                            TRTcells[data[0]].ET45m        = float(data[17])
                            TRTcells[data[0]].ET15         = float(data[18])
                            TRTcells[data[0]].ET15m        = float(data[19]) 
                            TRTcells[data[0]].VIL          = float(data[20])
                            TRTcells[data[0]].maxH         = float(data[21])
                            TRTcells[data[0]].maxHm        = float(data[22]) 
                            TRTcells[data[0]].POH          = int  (data[23])
                            TRTcells[data[0]].RANK         = int  (data[24])
                            TRTcells[data[0]].Dvel_x       = float(data[25])
                            TRTcells[data[0]].Dvel_y       = float(data[26])

                            # position in km (x->[-160:480], y->[255,965]
                            xCH = ccs4.WGStoCHx(float(data[3]),float(data[2])) # and x in the vertical
                            yCH = ccs4.WGStoCHy(float(data[3]),float(data[2])) # y is in horizontal direction

                            TRTcells[data[0]].iCH = (480-floor(xCH/1000.))  # i is in the vertical direction
                            TRTcells[data[0]].jCH = floor(yCH/1000.)-255    # j is in the horizontal dir. 

                            dx=float(data[ 4])
                            dy=float(data[ 5])
                            orientation = float(data[ 6])
                            rank = 0.1 * int (data[11])
                            add_cell(cell_mask, TRTcells[data[0]].iCH, TRTcells[data[0]].jCH, dx, 'empty_ellipse', dy=dy, angle=orientation, rank=rank)
                            #add_cell(cell_mask, TRTcells[data[0]].iCH, TRTcells[data[0]].jCH, dx, 'ellipse', dy=dy, angle=orientation, rank=rank)

                        # cell_contour_lon_lat=traj_ID + {data[0]:data[
        file.close()

        # max RANK is  2014072316000004 16 44 511.0 466.0
        print "                 traj_ID        Rank Area  x0   y0"
 
        if len(RANKs) > 0:
            RANKs_np = np.array(RANKs)
            i_max = np.argmax(RANKs_np)
            max_rank = RANKs_np[i_max]
            traj=traj_IDs[i_max]
            print '    max RANK is ', traj, TRTcells[traj].RANKr, TRTcells[traj].area, TRTcells[traj].iCH, TRTcells[traj].jCH

            cell_area_np=np.array(cell_area)
            i_max = np.argmax(cell_area_np)
            max_area = cell_area_np[i_max]
            traj=traj_IDs[i_max]
            print '    max area is ', traj, TRTcells[traj].RANKr, TRTcells[traj].area, TRTcells[traj].iCH, TRTcells[traj].jCH

        #print traj_IDs
        #
        #print traj, TRTcells[traj].date, TRTcells[traj].lon, TRTcells[traj].lat
        #traj=traj_IDs[1]
        #print traj, TRTcells[traj].date, TRTcells[traj].lon, TRTcells[traj].lat

    else:
        print '*** Warning, empty lightning input file ', file

    return traj_IDs, TRTcells, cell_mask

# input 
# prop  -> array that contains the data field in Swiss coord.
# i     -> position of the cell
# j     -> position of the cell
# form  -> form of the cell in the array, e.g. 
#          'square', 'circle', or 'ellipse'
# dx    -> length of form
#          if form='square' -> edge length
#          if form='circle' -> radius
#          if form='ellipse' -> larger semi axis
# dy    -> if form='ellipse' -> smaller semi axis
# angle -> if form='ellipse' -> orientation of ellipse

def add_cell(prop, i, j, dx, form, dy=1, angle=0, rank=1):

    from numpy import arange, meshgrid

    # print "add lightning at", i, j
    if form == 'square':
        # mark a square 
        prop[i-dx:i+dx,j-dx:j+dx]+=1
    elif form == 'circle':
        # mark a circle 
        rr=dx  # radius
        x = arange(0, 710)  # 709-
        y = arange(0, 640)
        X, Y = meshgrid(x, y)
        #print X
        #print Y        
        #print X.shape
        #print Y.shape
        #print lightnings.shape
        interior = ((X-j)**2 + (Y-i)**2) < rr**2
        prop += rank*interior
    elif form == 'empty_circle':
        # mark a circle 
        rr=dx  # radius
        x = arange(0, 710)  # 709-
        y = arange(0, 640)
        X, Y = meshgrid(x, y)
        from numpy import logical_and
        interior = logical_and( (rr-0.5)**2 < ((X-j)**2 + (Y-i)**2), ((X-j)**2 + (Y-i)**2) < (rr+0.5)**2 )
        prop += rank*interior
    elif form == 'ellipse':
        # mark a circle 
        x = arange(0, 710)  # 709-
        y = arange(0, 640)
        X, Y = meshgrid(x, y)
        #print X
        #print Y        
        #print X.shape
        #print Y.shape
        #print lightnings.shape
        from math import sin, cos, pi
        a_rad = 0.5*pi + (angle / 360. * 2. * pi)

        interior = ( ( ((X-j)*cos(a_rad) + (Y-i)*sin(a_rad)) / dx)**2 + 
                     ( ((X-j)*sin(a_rad) - (Y-i)*cos(a_rad)) / dy)**2) < 1
        prop += rank*interior
    elif form == 'empty_ellipse':
        # mark a circle 
        x = arange(0, 710)  # 709-
        y = arange(0, 640)
        X, Y = meshgrid(x, y)
        #print X
        #print Y        
        #print X.shape
        #print Y.shape
        #print lightnings.shape
        from math import sin, cos, pi
        a_rad = 0.5*pi + (angle / 360. * 2. * pi)

        from numpy import logical_and
        interior =  logical_and( 0.82 < ( ( ((X-j)*cos(a_rad) + (Y-i)*sin(a_rad)) / dx)**2 + 
                                          ( ((X-j)*sin(a_rad) - (Y-i)*cos(a_rad)) / dy)**2),
                                          ( ( ((X-j)*cos(a_rad) + (Y-i)*sin(a_rad)) / dx)**2 + 
                                          ( ((X-j)*sin(a_rad) - (Y-i)*cos(a_rad)) / dy)**2) < 1.18 )
        prop += rank*interior
    else:
        print '*** ERROR in add_lightning (lightning.py)'
        print '    unknown form ', form
