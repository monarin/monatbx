import cPickle as pickle
import sys

pickle_filename = sys.argv[1]
mypickle = pickle.load(open(pickle_filename,"rb"))
for key in mypickle.keys():
  print key
