import sys

file_a = open(sys.argv[1],'r')
data_a=file_a.read().split("\n")

search_word = 'MEAN DENSITY'
cn_line = 0
for line in data_a:
  if line.find(search_word) > 0:
    break
  cn_line += 1
line_start = cn_line + 2

search_word = 'Summary of fit of model to map'
cn_line = 0
for line in data_a:
  if line.find(search_word) > 0:
    break
  cn_line += 1
line_end = cn_line - 4

aa_res_cc_dict = dict([('ALA', 0), ('ARG',0), ('ASN',0), ('ASP',0), ('ASX',0), ('CYS',0), ('GLU',0), ('GLN',0), \
  ('GLX',0), ('GLY',0), ('HIS',0), ('ILE',0), ('LEU',0), ('LYS',0), ('MET',0), ('PHE',0), ('PRO',0), ('SER',0), ('THR',0), \
  ('TRP',0), ('TYR',0), ('VAL',0), ('CA',0)])
aa_res_mean_den_all_dict = dict([('ALA', 0), ('ARG',0), ('ASN',0), ('ASP',0), ('ASX',0), ('CYS',0), ('GLU',0), ('GLN',0), \
  ('GLX',0), ('GLY',0), ('HIS',0), ('ILE',0), ('LEU',0), ('LYS',0), ('MET',0), ('PHE',0), ('PRO',0), ('SER',0), ('THR',0), \
  ('TRP',0), ('TYR',0), ('VAL',0), ('CA',0)])
aa_res_mean_den_main_dict = dict([('ALA', 0), ('ARG',0), ('ASN',0), ('ASP',0), ('ASX',0), ('CYS',0), ('GLU',0), ('GLN',0), \
  ('GLX',0), ('GLY',0), ('HIS',0), ('ILE',0), ('LEU',0), ('LYS',0), ('MET',0), ('PHE',0), ('PRO',0), ('SER',0), ('THR',0), \
  ('TRP',0), ('TYR',0), ('VAL',0), ('CA',0)])
aa_res_mean_den_side_dict = dict([('ALA', 0), ('ARG',0), ('ASN',0), ('ASP',0), ('ASX',0), ('CYS',0), ('GLU',0), ('GLN',0), \
  ('GLX',0), ('GLY',0), ('HIS',0), ('ILE',0), ('LEU',0), ('LYS',0), ('MET',0), ('PHE',0), ('PRO',0), ('SER',0), ('THR',0), \
  ('TRP',0), ('TYR',0), ('VAL',0), ('CA',0)])
aa_res_cn_dict = dict([('ALA', 0), ('ARG',0), ('ASN',0), ('ASP',0), ('ASX',0), ('CYS',0), ('GLU',0), ('GLN',0), \
  ('GLX',0), ('GLY',0), ('HIS',0), ('ILE',0), ('LEU',0), ('LYS',0), ('MET',0), ('PHE',0), ('PRO',0), ('SER',0), ('THR',0), \
  ('TRP',0), ('TYR',0), ('VAL',0), ('CA',0)])
for i in range(line_start, line_end):
  row = data_a[i].split()
  if row[0] in aa_res_cc_dict:
    aa_res_cc_dict[row[0]] += float(row[3])
    aa_res_mean_den_all_dict[row[0]] += float(row[4])
    aa_res_mean_den_main_dict[row[0]] += float(row[5])
    aa_res_mean_den_side_dict[row[0]] += float(row[6])
    aa_res_cn_dict[row[0]] += 1

for key in aa_res_cc_dict.keys():
  if aa_res_cn_dict[key] > 0:
    aa_res_cc_dict[key] = aa_res_cc_dict[key]/aa_res_cn_dict[key]
    aa_res_mean_den_all_dict[key] = aa_res_mean_den_all_dict[key]/aa_res_cn_dict[key]
    aa_res_mean_den_main_dict[key] = aa_res_mean_den_main_dict[key]/aa_res_cn_dict[key]
    aa_res_mean_den_side_dict[key] = aa_res_mean_den_side_dict[key]/aa_res_cn_dict[key]

  print key, aa_res_cn_dict[key], aa_res_cc_dict[key], aa_res_mean_den_all_dict[key], aa_res_mean_den_main_dict[key], aa_res_mean_den_side_dict[key]
