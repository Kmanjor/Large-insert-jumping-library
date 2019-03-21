from collections import defaultdict
import pprint
import numpy as np
from Bio import SeqIO


amp_tag_start = defaultdict(list)
amp_tag_end = defaultdict(list)
ves_start = defaultdict(list)
ves_end = defaultdict(list)

f_blastn = open('./gu_t.blastn', 'r')
for line in f_blastn:
    line = line.strip('\n').split('\t')
    if line[2] == 'amp_tag':
        amp_tag_start[line[0]].append(int(line[4]))
        amp_tag_end[line[0]].append(int(line[5]))
    else:
        ves_start[line[0]].append(int(line[4]))
        ves_end[line[0]].append(int(line[5]))
f_blastn.close()

# pprint.pprint(amp_tag)
# pprint.pprint(ves)
index_position = {}
for i in amp_tag_start:
    amp_start = np.array(amp_tag_start[i]).min()
    amp_end = np.array(amp_tag_end[i]).max()
    if i in ves_start:
        left = []
        right = []
        for n in ves_end[i]:
            if n <= amp_start:
                left.append(n)
        for m in ves_start[i]:
            if m >= amp_end:
                right.append(m)
        if left != []:
            left_end = np.array(left).max()
        else:
            left_end = None
        if right != []:
            right_start = np.array(right).min()
        else:
            right_start = None
    else:
        left_end = None
        right_start = None
    index_position[i] = [left_end, amp_start, amp_end, right_start]
# print(len(index_position))
# print(len(amp_tag_start))
# print(len(ves_start))
# pprint.pprint(index_position)
seq_dict = SeqIO.to_dict(SeqIO.parse('/public/home/tongli/G_Y/CCS/gu_t/m54045_190112_025400.ccs.fasta', 'fasta'))

left_list = []
right_list = []
for i in index_position:
    record = seq_dict[i]
    per = index_position[i]
    if per[0] == None and per[3] == None:
        if per[1] >= 50 and len(record.seq) - per[2] >= 50:
            fosmid_left = record[:per[1]-1]
            fosmid_right = record[per[2]:]
            fosmid_left.id, fosmid_left.name, fosmid_left.description = fosmid_left.id + '/1', '', ''
            fosmid_right.id, fosmid_right.name, fosmid_right.description = fosmid_right.id + '/2', '', ''
            left_list.append(fosmid_left)
            right_list.append(fosmid_right)
    if per[0] != None and per[3] != None:
        if per[1] - per[0] >= 50 and per[3] - per[2] >= 50:
            fosmid_left = record[per[0]:per[1] - 1]
            fosmid_right = record[per[2]:per[3]-1]
            fosmid_left.id, fosmid_left.name, fosmid_left.description = fosmid_left.id + '/1', '', ''
            fosmid_right.id, fosmid_right.name, fosmid_right.description = fosmid_right.id + '/2', '', ''
            left_list.append(fosmid_left)
            right_list.append(fosmid_right)
    if per[0] == None and per[3] != None:
        if per[1] >= 50 and per[3] - per[2] >= 50:
            fosmid_left = record[:per[1] - 1]
            fosmid_right = record[per[2]:per[3] - 1]
            fosmid_left.id, fosmid_left.name, fosmid_left.description = fosmid_left.id + '/1', '', ''
            fosmid_right.id, fosmid_right.name, fosmid_right.description = fosmid_right.id + '/2', '', ''
            left_list.append(fosmid_left)
            right_list.append(fosmid_right)
    if per[0] != None and per[3] == None:
        if per[1] - per[0] >= 50 and len(record.seq) - per[2] >= 50:
            fosmid_left = record[per[0]:per[1] - 1]
            fosmid_right = record[per[2]:]
            fosmid_left.id, fosmid_left.name, fosmid_left.description = fosmid_left.id + '/1', '', ''
            fosmid_right.id, fosmid_right.name, fosmid_right.description = fosmid_right.id + '/2', '', ''
            left_list.append(fosmid_left)
            right_list.append(fosmid_right)

SeqIO.write(left_list, 'gu_t_left.fa', 'fasta')
SeqIO.write(right_list, 'gu_t_right.fa', 'fasta')
