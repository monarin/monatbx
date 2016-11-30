import sys
import numpy as np
import psana
from Detector.PyDetector import PyDetector
import argparse
from matplotlib import pyplot as plt
import RegDB.experiment_info

from xfel.cxi.cspad_ana import cspad_tbx
from libtbx import easy_pickle as ep
from scitbx.array_family import flex

plt.ion()

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run")
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--nevt", help="number of events", type=int)
args = parser.parse_args()
if not args.run:
    if len(sys.argv)>1:
        run=sys.argv[1]
    else:
        run='%s'%(RegDB.experiment_info.experiment_runs('XPP')[-1]['num'])
else:
    run=args.run
if not args.exp:
    expname=RegDB.experiment_info.active_experiment('XPP')[1]
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/ffb/xpp/%s/xtc:live'%expname
else:
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/psdm/xpp/%s/xtc:live'%expname

print 'getting images from dataset: ',dsname
ds = psana.DataSource(dsname)
det = PyDetector(psana.Source('rayonix'),ds.env())

print "WORKING..."

images=[]
for ievt,evt in enumerate(ds.events()):
    if det.raw(evt) is not None:
        image = det.raw(evt)
        image = flex.int(image.astype(np.intc))
        lambda_wavelength = 12398./(ds.env().epicsStore().value("SIOC:SYS0:ML00:AO627"))
        pixelSize = (170./3840.)*2. #standard 2x2 binning
        saturation=65535 #maybe thus the trustable value should be lower?

        eventTime=evt.get(psana.EventId).time()[0]+evt.get(psana.EventId).time()[1]*1e-9

        # Determine timestamp using cspad_tbx
        evt_time = cspad_tbx.evt_time(evt)
        evt_timestamp = cspad_tbx.evt_timestamp(evt_time)

        #this assumes that we will get the correct numbers for one position
        #then we will use the robot position to get relative numbers assuming
        #the beam centr stays about the same
        robX=ds.env().epicsStore().value("robot_x")
        robY=ds.env().epicsStore().value("robot_y")
        robZ=ds.env().epicsStore().value("robot_z")

        #beamX=602.015+(ds.env().epicsStore().value("robot_x"))
        #beamY=51.682+(ds.env().epicsStore().value("robot_y"))
        #beamDist=-563.083+60.+ (ds.env().epicsStore().value("robot_z"))

        beamX = 969
        beamY = 968
        beamDist = 300

        # Assemble data pack
        data = cspad_tbx.dpack(data=image,
                               distance=beamDist,
                               pixel_size=pixelSize,
                               wavelength=lambda_wavelength,
                               beam_center_x=beamX,
                               beam_center_y=beamY,
                               ccd_image_saturation=saturation,
                               saturated_value=saturation,
                               timestamp = evt_timestamp)

        # output image pickle
        fpklname='/reg/d/psdm/xpp/%s/ftc/Run%s/Run_%s_idx_{}.pickle'%(expname,run,run)
        ep.dump(fpklname.format(ievt), data)

        print 'Importing image {}, beamX = {}, beamY = {}, Z-offset = {}, wavelength = {} '.format(fplkname, beamX, beamY, beamDist, lambda_wavelength)

    if ievt>5:
        break

print 'is done '
