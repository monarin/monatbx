import os
import sys
import shutil

file_a = open(sys.argv[1],'r')
data_a=file_a.read().split("\n")
print len(data_a)
destdir = sys.argv[2]

cn_file = 0
for line in data_a:
  line_arr = line.split('/')
  if len(line_arr) > 0:
    filename = line_arr[len(line_arr)-3]+'_'+line_arr[len(line_arr)-1]
    if True:
      shutil.copy(line, destdir+'/'+filename)
      print 'copy ', line, 'to', filename
      cn_file += 1
    #except IOError:
    #  print line
