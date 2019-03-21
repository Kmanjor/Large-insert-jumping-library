from Bio import SeqIO
import pandas as pd


chr2num = {}
# 染色体名字转化为数字
chrLen = open('./data_test/s288c_ChrLen.txt', 'r')
n = 0
for line in chrLen:
    per = line.split(' ')[0]
    n += 1
    chr2num[per] = n
chrLen.close()


df = pd.DataFrame(columns=('scaff', 'chr', 'start', 'strand'))
# 输入blastn比对结果
blastn = open('./data_test7/yeast_30_20.blastn', 'r')
for line in blastn:
    line = line.strip('\n').split('\t')
    if line[0] not in df['scaff'].values:
        if line[7] == 'plus':
            s = pd.Series({'scaff': line[0], 'chr': int(chr2num[line[2]]), 'start': int(line[8]), 'strand': line[7]})
            df = df.append(s, ignore_index=True)
        else:
            s = pd.Series({'scaff': line[0], 'chr': int(chr2num[line[2]]), 'start': int(line[9]), 'strand': line[7]})
            df = df.append(s, ignore_index=True)
blastn.close()
dfs = df.sort_values(by=['chr', 'start'])


# 输入原始contigs序列文件
scaff_seq = SeqIO.to_dict(SeqIO.parse('./data_test7/yeast_30_20.final.scaffolds.fasta', 'fasta'))
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
SeqIO.write(shuchu, './data_test7/yeast_30_20_Ss.fa', 'fasta')
