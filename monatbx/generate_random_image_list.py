import os
import sys
import random

p = sys.argv[1]
n_images = int(sys.argv[2])

frame_files  = []
if os.path.isdir(p):
  for pickle_filename in os.listdir(p):
    if pickle_filename.endswith('.pickle'):
      frame_files.append(p+'/'+pickle_filename)

i_rand = random.sample(range(len(frame_files)),n_images)
frame_files_sel = [frame_files[i] for i in i_rand]

txt_out = ''
for frame in frame_files_sel:
  txt_out += frame + '\n'

f = open('frame_rand_'+str(n_images)+'.lst', 'w')
f.write(txt_out)
f.close()

