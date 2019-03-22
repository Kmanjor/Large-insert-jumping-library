from collections import defaultdict
import pprint
import networkx as nx


def paired_map():
    return [0, 0, 0, None, 0, None]


def paired_pos_info():
    return [None, None, None, None]

# minimap2 两个末端分别比对
allreadInfo = defaultdict(paired_map)
f_left_sam = open('./data_test7/guzi_m_left_mmi_q60F2048_Abb.sam', 'r')
for line in f_left_sam:
    line = line.strip('\n').split('\t')
    allreadInfo['/'.join(line[0].split('/')[:-1])][0] = 1
    allreadInfo['/'.join(line[0].split('/')[:-1])][2] = line[1]
    allreadInfo['/'.join(line[0].split('/')[:-1])][3] = line[2]
f_left_sam.close()
f_right_sam = open('./data_test7/guzi_m_right_mmi_q60F2048_Abb.sam', 'r')
for line in f_right_sam:
    line = line.strip('\n').split('\t')
    allreadInfo['/'.join(line[0].split('/')[:-1])][1] = 1
    allreadInfo['/'.join(line[0].split('/')[:-1])][4] = line[1]
    allreadInfo['/'.join(line[0].split('/')[:-1])][5] = line[2]
f_right_sam.close()

unambi_sam = set()
n = 0
for i in allreadInfo:
    readInfo = allreadInfo[i]
    if readInfo[0] + readInfo[1] == 2:
        n += 1
        unambi_sam.add(i)
print(n)

# minimap2 两个末端分别比对
allbedInfo = defaultdict(paired_pos_info)
f_bed_left = open('./data_test7/guzi_m_left_mmi_q60F2048.bed', 'r')
for line in f_bed_left:
    line = line.strip('\n').split('\t')
    name = '/'.join(line[3].split('/')[:-1])
    allbedInfo[name][0] = line[1]
    allbedInfo[name][1] = line[2]
f_bed_left.close()
f_bed_right = open('./data_test7/guzi_m_right_mmi_q60F2048.bed', 'r')
for line in f_bed_right:
    line = line.strip('\n').split('\t')
    name = '/'.join(line[3].split('/')[:-1])
    allbedInfo[name][2] = line[1]
    allbedInfo[name][3] = line[2]
f_bed_right.close()

unambi_bed = set()
m = 0
for i in allbedInfo:
    bedInfo = allbedInfo[i]
    if None not in bedInfo:
        m += 1
        unambi_bed.add(i)
print(m)

print(unambi_bed == unambi_sam)

# 得到newform
f_newform = []
for i in unambi_sam:
    newline = [i] + allreadInfo[i][2:] + allbedInfo[i]
    f_newform.append(newline)
# pprint.pprint(f_newform)


# 各自方向去重复。正正、负负、正负、负正
def mvdup(newform, stand1, stand2, qstart, qend):
    strand2strand = defaultdict(list)
    for line in newform:
        strand1 = (int(line[1]) >> (5 - 1)) & 1
        strand2 = (int(line[3]) >> (5 - 1)) & 1
        if strand1 == stand1 and strand2 == stand2:
            chr1 = line[2]
            chr2 = line[4]
            chrComb = chr1 + '/' + chr2
            strand2strand[chrComb].append(line)
    cluster_graph = nx.Graph()
    for line in newform:
        strand1 = (int(line[1]) >> (5 - 1)) & 1
        strand2 = (int(line[3]) >> (5 - 1)) & 1
        if strand1 == stand1 and strand2 == stand2:
            chr1 = line[2]
            chr2 = line[4]
            chrComb = chr1 + '/' + chr2
            perChr = strand2strand[chrComb]
            q_pos_start = int(line[qstart])
            q_pos_end = int(line[qend])
            for per in perChr:
                s_pos_start = int(per[qstart])
                s_pos_end = int(per[qend])
                startGap = abs(q_pos_start - s_pos_start)
                endGap = abs(q_pos_end - s_pos_end)
                if startGap <= 20 and endGap <= 20:
                    cluster_graph.add_edge(tuple(line), tuple(per))
    newclusters = list(nx.connected_components(cluster_graph))
    f_result = []
    for cluster in newclusters:
        cluster = list(cluster)
        maxRead = cluster[0]
        comesfrom = []
        for read in cluster:
            comesfrom.append(read[0])
            maxLength = (int(maxRead[6]) - int(maxRead[5])) + (int(maxRead[8]) - int(maxRead[7]))
            readLength = (int(read[6]) - int(read[5])) + (int(read[8]) - int(read[7]))
            if readLength > maxLength:
                maxRead = read
        f_result.append(maxRead)
    return f_result

# 各自方向去重复。正正、负负、正负、负正
pp2pp = mvdup(f_newform, 0, 0, 5, 8)
mm2mm = mvdup(f_newform, 1, 1, 6, 7)
pp2mm = mvdup(f_newform, 0, 1, 5, 7)
mm2pp = mvdup(f_newform, 1, 0, 6, 8)
print(len(pp2pp), len(mm2mm), len(pp2mm), len(mm2pp))

# print(pp2pp)
# print(mm2mm)

# 正负链去重复，因为测序有两个方向。正正负负、正负正负、负正负正。
def mvstrand(newform1, newform2, sstart, send, qstart, qend):
    strand2strand_1 = defaultdict(list)
    for line in newform1:
        chr1 = line[2]
        chr2 = line[4]
        chrComb = chr1 + '/' + chr2
        strand2strand_1[chrComb].append(line)
    strand2strand_2 = defaultdict(list)
    for line in newform2:
        chr1 = line[2]
        chr2 = line[4]
        chrComb = chr1 + '/' + chr2
        strand2strand_2[chrComb].append(line)
    cluster_graph = nx.Graph()
    for line in newform2:
        chr1 = line[2]
        chr2 = line[4]
        chrComb = chr2 + '/' + chr1
        if strand2strand_1[chrComb]:
            perChr = strand2strand_1[chrComb]
            q_pos_start = int(line[qstart])
            q_pos_end = int(line[qend])
            for per in perChr:
                s_pos_start = int(per[sstart])
                s_pos_end = int(per[send])
                startGap = abs(q_pos_start - s_pos_end)
                endGap = abs(q_pos_end - s_pos_start)
                if startGap <= 20 and endGap <= 20:
                    cluster_graph.add_edge(tuple(line), tuple(per))
    for perChr in strand2strand_2:
        for per in strand2strand_2[perChr]:
            if tuple(per) not in cluster_graph.nodes:
                cluster_graph.add_edge(tuple(per), tuple(per))
    for perChr in strand2strand_1:
        for per in strand2strand_1[perChr]:
            if tuple(per) not in cluster_graph.nodes:
                cluster_graph.add_edge(tuple(per), tuple(per))
    newclusters = list(nx.connected_components(cluster_graph))
    f_result = []
    for cluster in newclusters:
        cluster = list(cluster)
        maxRead = cluster[0]
        comesfrom = []
        for read in cluster:
            comesfrom.append(read[0])
            maxLength = (int(maxRead[6]) - int(maxRead[5])) + (int(maxRead[8]) - int(maxRead[7]))
            readLength = (int(read[6]) - int(read[5])) + (int(read[8]) - int(read[7]))
            if readLength > maxLength:
                maxRead = read
        f_result.append(maxRead)
    return f_result

# 正负链去重复，因为测序有两个方向。正正负负、正负正负、负正负正。
mv_strand_ppmm = mvstrand(pp2pp, mm2mm, 5, 8, 6, 7)
mv_strand_pmpm = mvstrand(pp2mm, pp2mm, 5, 7, 5, 7)
mv_strand_mpmp = mvstrand(mm2pp, mm2pp, 6, 8, 6, 8)
print(len(mv_strand_ppmm))
print(len(mv_strand_pmpm))
print(len(mv_strand_mpmp))
# print(mv_strand_ppmm)
# print(mv_strand_pmpm)
# print(mv_strand_mpmp)


# 输出newform验证。
f_shuchu = open('./data_test7/guzi_m_mmi_q60F2048.newform', 'w')
for i in mv_strand_ppmm:
    waitWrite = '\t'.join(list(i)[0:9]) + '\n'
    f_shuchu.write(waitWrite)
for i in mv_strand_pmpm:
    waitWrite = '\t'.join(list(i)[0:9]) + '\n'
    f_shuchu.write(waitWrite)
for i in mv_strand_mpmp:
    waitWrite = '\t'.join(list(i)[0:9]) + '\n'
    f_shuchu.write(waitWrite)
f_shuchu.close()
