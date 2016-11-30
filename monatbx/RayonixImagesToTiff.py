import sys
import numpy as np
import psana
from Detector.PyDetector import PyDetector
import argparse
from matplotlib import pyplot as plt
import RegDB.experiment_info
import Image

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="run")
parser.add_argument("--exp", help="experiment name")
parser.add_argument("--nevt", help="number of events", type=int)
args = parser.parse_args()
if not args.run:
    if len(sys.argv)>1:
        run=sys.argv[1]
    else:
        run=raw_input('please enter a run#: ')
else:
    run=args.run
if not args.exp:
    expname=RegDB.experiment_info.active_experiment('XPP')[1]
    dsname='exp='+expname+':run='+run+':smd:dir=/reg/d/ffb/xpp/%s/xtc'%expname
else:
    expname=args.exp
    dsname='exp=%s:run=%s:smd:dir=/reg/d/psdm/xpp/%s/xtc'%(args.exp,run,args.exp)

print 'getting images from dataset: ',dsname
ds = psana.DataSource(dsname)
det = PyDetector(psana.Source('rayonix'),ds.env())

images=[]
for ievt,evt in enumerate(ds.events()):
    if det.raw(evt) is not None:
        image = det.raw(evt)
        # output image pickle
        tiffname='/reg/d/psdm/xpp/%s/ftc/Run%s/Run_%s_%i-%i-%i.tiff'%(expname,run,run,evt.get(psana.EventId).time()[0],evt.get(psana.EventId).time()[1],evt.get(psana.EventId).fiducials())
        im = Image.fromarray(image)
        im.save(tiffname)

    if ievt>5:
        break

print 'is done '
