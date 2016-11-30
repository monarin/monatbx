import cPickle as pickle
import sys

sol_fname = sys.argv[1]
sol_pickle = pickle.load(open(sol_fname, "rb"))

ind_fname = sys.argv[2]
ind_pickle = pickle.load(open(ind_fname, "rb"))
cn_match = 0
for key in sol_pickle:
  if key in ind_pickle:
    if sol_pickle[key] == ind_pickle[key]:
      cn_match += 1
    else:
      print key, sol_pickle[key], ind_pickle[key]

print len(sol_pickle.keys()), cn_match
