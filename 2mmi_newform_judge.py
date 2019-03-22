from collections import defaultdict
import pprint


def paired_map():
    return [0, 0, 0, None, 0, None]


def paired_pos_info():
    return [None, None, None, None]


# minimap2 两个末端分别比对
allreadInfo = defaultdict(paired_map)
f_left_sam = open('./data_test/yeast_1_left_mmi_Abb.sam', 'r')
for line in f_left_sam:
    line = line.strip('\n').split('\t')
    allreadInfo['/'.join(line[0].split('/')[:-1])][0] = 1
    allreadInfo['/'.join(line[0].split('/')[:-1])][2] = line[1]
    allreadInfo['/'.join(line[0].split('/')[:-1])][3] = line[2]
f_left_sam.close()
f_right_sam = open('./data_test/yeast_1_right_mmi_Abb.sam', 'r')
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
# print(n)
# print(len(unambi_sam))
print(len(allreadInfo), 'total pair')

# minimap2 两个末端分别比对
allbedInfo = defaultdict(paired_pos_info)
f_bed_left = open('./data_test/yeast_1_left_mmi.bed', 'r')
for line in f_bed_left:
    line = line.strip('\n').split('\t')
    name = '/'.join(line[3].split('/')[:-1])
    allbedInfo[name][0] = line[1]
    allbedInfo[name][1] = line[2]
f_bed_left.close()
f_bed_right = open('./data_test/yeast_1_right_mmi.bed', 'r')
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
print(m, '双端比对上')
# print(len(unambi_bed))
# print(len(allbedInfo))
# print(unambi_bed == unambi_sam)

# 得到newform
f_newform = []
for i in unambi_sam:
    newline = [i] + allreadInfo[i][2:] + allbedInfo[i]
    f_newform.append(newline)
# pprint.pprint(f_newform)

# 双端均未比对上的
unaligned_pair = set()
for line in f_newform:
    num1 = (int(line[1]) >> (3 - 1)) & 1
    num2 = (int(line[3]) >> (3 - 1)) & 1
    if num1 == 1 and num2 == 1:
        unaligned_pair.add(line[0])
print(len(unaligned_pair), '双端均未比对上')
# print(unaligned_pair == unambi_bed)
# print(unaligned_pair)

# 一端比对上的
onlyaligned = set()
for line in f_newform:
    num1 = (int(line[1]) >> (3 - 1)) & 1
    num2 = (int(line[3]) >> (3 - 1)) & 1
    if num1 + num2 == 1:
        onlyaligned.add(line[0])
print(len(onlyaligned), '一端比对上')

# 嵌合比对
chimeric = set()
for line in f_newform:
    num1 = (int(line[1]) >> (12 - 1)) & 1
    num2 = (int(line[3]) >> (12 - 1)) & 1
    if num1 == 1 or num2 == 1:
        chimeric.add(line[0])
print(len(chimeric), '嵌合比对')
# print(chimeric)

# 多位点比对
multiple = set()
for line in f_newform:
    num1 = (int(line[1]) >> (9 - 1)) & 1
    num2 = (int(line[3]) >> (9 - 1)) & 1
    if num1 == 1 or num2 == 1:
        multiple.add(line[0])
print(len(multiple), '多位点比对')