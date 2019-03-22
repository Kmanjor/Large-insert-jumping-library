from collections import defaultdict


def paired_map():
    return [0, 0, 0, None, 0, None]


def paired_pos_info():
    return [None, None, None, None]

# minimap2 两个末端分别比对
allreadInfo = defaultdict(paired_map)
f_left_sam = open('./data_test/yeast_1_left_mmi_q60F2048_Abb.sam', 'r')
for line in f_left_sam:
    line = line.strip('\n').split('\t')
    allreadInfo['/'.join(line[0].split('/')[:-1])][0] = 1
    allreadInfo['/'.join(line[0].split('/')[:-1])][2] = line[1]
    allreadInfo['/'.join(line[0].split('/')[:-1])][3] = line[2]
f_left_sam.close()
f_right_sam = open('./data_test/yeast_1_right_mmi_q60F2048_Abb.sam', 'r')
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
f_bed_left = open('./data_test/yeast_1_left_mmi_q60F2048.bed', 'r')
for line in f_bed_left:
    line = line.strip('\n').split('\t')
    name = '/'.join(line[3].split('/')[:-1])
    allbedInfo[name][0] = line[1]
    allbedInfo[name][1] = line[2]
f_bed_left.close()
f_bed_right = open('./data_test/yeast_1_right_mmi_q60F2048.bed', 'r')
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


link2Chr = set()
for line in f_newform:
    if line[2] != line[4]:
        link2Chr.add(line[0])
        # print(line[0])
print(len(link2Chr), 'link2chr')

inverted = set()
for line in f_newform:
    if line[0] not in link2Chr:
        num1 = (int(line[1]) >> (5 - 1)) & 1
        num2 = (int(line[3]) >> (5 - 1)) & 1
        if num1 != num2:
            inverted.add(line[0])
            # print(line[0])
print(len(inverted), 'inverted')


tandem = set()
for line in f_newform:
    if line[0] not in link2Chr:
        if line[0] not in inverted:
            num = (int(line[1]) >> (5 - 1)) & 1
            if num == 0:
                if int(line[5]) > int(line[7]):
                    tandem.add(line[0])
            if num == 1:
                if int(line[6]) < int(line[8]):
                    tandem.add(line[0])
print(len(tandem), 'tandem')


spacing100kb = set()
for line in f_newform:
    if line[0] not in link2Chr:
        if line[0] not in inverted:
            if line[0] not in tandem:
                num = (int(line[1]) >> (5 - 1)) & 1
                if num == 0:
                    if int(line[8]) - int(line[5])> 100000:
                        spacing100kb.add(line[0])
                        # print(line[0])
                if num == 1:
                    if int(line[6]) - int(line[7])> 100000:
                        spacing100kb.add(line[0])
print(len(spacing100kb), 'spacing 100kb')


more1less20more50kb = set()
for line in f_newform:
    if line[0] not in link2Chr:
        if line[0] not in inverted:
            if line[0] not in tandem:
                # if line[0] not in spacing100kb:
                num = (int(line[1]) >> (5 - 1)) & 1
                if num == 0:
                    temp = int(line[8]) - int(line[5])
                    if 1000 < temp < 20000 or 50000 < temp <= 100000:
                        more1less20more50kb.add(line[0])
                if num == 1:
                    temp = int(line[6]) - int(line[7])
                    if 1000 < temp < 20000 or 50000 < temp <= 100000:
                        more1less20more50kb.add(line[0])
print(len(more1less20more50kb), 'more1less20more50kb')
# print(more1less20more50kb)

less1kb = set()
for line in f_newform:
    if line[0] not in link2Chr:
        if line[0] not in inverted:
            if line[0] not in tandem:
                if line[0] not in spacing100kb:
                    num = (int(line[1]) >> (5 - 1)) & 1
                    if num == 0:
                        temp = int(line[8]) - int(line[5])
                        if 1000 >= temp:
                            less1kb.add(line[0])
                    if num == 1:
                        temp = int(line[6]) - int(line[7])
                        if 1000 >= temp:
                            less1kb.add(line[0])
print(len(less1kb), 'less1kb')


allRead = set()
rightPlace = set()
for line in f_newform:
    allRead.add(line[0])
    if line[2] == line[4]:
        orit1 = (int(line[1]) >> (5 - 1)) & 1
        orit2 = (int(line[3]) >> (5 - 1)) & 1
        if orit1 == orit2:
            if orit1 == 0:
                temp = int(line[8]) - int(line[5])
                if 20000 <= temp <= 50000:
                    rightPlace.add(line[0])
            if orit1 == 1:
                temp = int(line[6]) - int(line[7])
                if 20000 <= temp <= 50000:
                    rightPlace.add(line[0])
print(len(rightPlace), 'right placement')

consist = link2Chr | tandem | inverted | spacing100kb | less1kb | more1less20more50kb | rightPlace
print(consist == allRead)