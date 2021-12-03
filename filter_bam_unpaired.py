#!/usr/bin/env python
#coding=utf-8

import pysam
import sys

#### READ FILE ####
try:
    in_bam, out_bam = sys.argv[1:3]
except:
    print("Usage: {} <in_bam> <out_filtered_bam>".format(sys.argv[0]))
    sys.exit()


samfile = pysam.AlignmentFile(in_bam,'rb')
tmpfile = pysam.AlignmentFile(out_bam, "wb", template=samfile)

lineCount = 0
max_hit_num = 3
pre_read_id = ""
cur_read_list = list()
for read in samfile:
    cur_read_id = read.query_name
    if cur_read_id == pre_read_id:
        cur_read_list.append(read)
    else:
        if 1 < len(cur_read_list) < max_hit_num:
            for cur_read in cur_read_list:
                tmpfile.write(cur_read)
        cur_read_list.clear()
        cur_read_list.append(read)
        pre_read_id = cur_read_id
