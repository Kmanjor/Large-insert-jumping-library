import pprint
import pandas as pd
from Bio import SeqIO


chr2num = {}
# 染色体名字转化为数字
chrLen = open('./data_test/yugu1_ChrLen.txt', 'r')
n = 0
for line in chrLen:
    per = line.split(' ')[0]
    n += 1
    chr2num[per] = n
chrLen.close()
# pprint.pprint(chr2num)


# 输入 nucmer 的 mcoords 结果
region_max_dict = {}
f_mcoords = open('./data_test4/guzi_sspacek1.mcoords', 'r')
for line in f_mcoords:
    line = line.strip('\n').split('\t')
    if line[12] not in region_max_dict:
        region_max_dict[line[12]] = [int(line[0]), int(line[1]), int(line[2]), int(line[3]), line[11]]
    else:
        map_length_old = region_max_dict[line[12]][1] - region_max_dict[line[12]][0] + 1
        map_length_new = int(line[1]) - int(line[0]) + 1
        if map_length_new > map_length_old:
            region_max_dict[line[12]] = [int(line[0]), int(line[1]), int(line[2]), int(line[3]), line[11]]
f_mcoords.close()
print(len(region_max_dict))
# pprint.pprint(region_max_dict)

# 利用pandas排序
df = pd.DataFrame(columns=('scaff', 'chr', 'start', 'strand'))
for i in region_max_dict:
    per = region_max_dict[i]
    if per[2] < per[3]:
        s = pd.Series({'scaff': i, 'chr': int(chr2num[per[4]]), 'start': int(per[0]), 'strand': 'plus'})
        df = df.append(s, ignore_index=True)
    else:
        s = pd.Series({'scaff': i, 'chr': int(chr2num[per[4]]), 'start': int(per[0]), 'strand': 'minus'})
        df = df.append(s, ignore_index=True)
dfs = df.sort_values(by=['chr', 'start'])
# print(dfs)


# 输入原始contigs序列文件
scaff_seq = SeqIO.to_dict(SeqIO.parse('./data_test6/yeast_pb50.final.scaffolds.fasta', 'fasta'))
shuchu = []
n = 0
for per in dfs.values:
    oldseq = scaff_seq[per[0]]
    n += 1
    if per[3] == 'plus':
        newseq = oldseq
        newseq.id = 'contig%s' % n
        newseq.name = 'contig%s' % n
        newseq.description = 'contig%s' % n
        shuchu.append(newseq)
    else:
        newseq = oldseq.reverse_complement()
        newseq.id = 'contig%s' % n
        newseq.name = 'contig%s' % n
        newseq.description = 'contig%s' % n
        shuchu.append(newseq)
print(len(shuchu))
# 输出排序更名后的序列文件
SeqIO.write(shuchu, './data_test6/yeast_pb50_2_scaff_Ss.fa', 'fasta')