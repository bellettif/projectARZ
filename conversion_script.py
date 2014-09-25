import re
import string

import sys

if len(sys.argv) == 0:
    print 'Need at least one file name'
    raise Exception('No file name given')

file_list = sys.argv[1:]

for file_name in file_list:
    input_file_name = file_name
    output_file_name = ''.join(file_name.split('.')[:-1]) + '.csv'
    header = ['vehicule_ID',
              'frame_ID',
              'tot_frames',
              'time_since_epoch_ms',
              'local_x',
              'local_y',
              'global_x',
              'global_y',
              'veh_length',
              'veh_width',
              'veh_class',
              'veh_v',
              'veh_acc',
              'lane_id',
              'prec_veh',
              'follow_veh',
              'spacing',
              'headway']
    input_file = open(input_file_name, 'r')
    output_file = open(output_file_name, 'w')
    input_file.readline() # skip first line
    output_file.write(','.join(header) + '\n')
    for line in input_file:
        output_file.write(re.sub(' +', ',', line.strip(), re.UNICODE) + '\n')
    input_file.close()
    output_file.close()