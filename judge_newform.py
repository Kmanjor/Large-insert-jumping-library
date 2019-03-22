link2Chr = set()
f_newBed = open('./data_test/sindFosEnd_bwa_q60.newform', 'r')
for line in f_newBed:
    line = line.strip('\n').split('\t')
    if line[2] != line[4]:
        link2Chr.add(line[0])
        # print(line[0])
f_newBed.close()
print(len(link2Chr))


inverted = set()
f_newBed = open('./data_test/sindFosEnd_bwa_q60.newform', 'r')
for line in f_newBed:
    line = line.strip('\n').split('\t')
    if line[0] not in link2Chr:
        num1 = (int(line[1]) >> (5 - 1)) & 1
        num2 = (int(line[3]) >> (5 - 1)) & 1
        if num1 != num2:
            inverted.add(line[0])
            # print(line[0])
f_newBed.close()
print(len(inverted))


tandem = set()
f_newBed = open('./data_test/sindFosEnd_bwa_q60.newform', 'r')
for line in f_newBed:
    line = line.strip('\n').split('\t')
    if line[0] not in link2Chr:
        if line[0] not in inverted:
            num = (int(line[1]) >> (5 - 1)) & 1
            if num == 0:
                if int(line[5]) > int(line[7]):
                    tandem.add(line[0])
            if num == 1:
                if int(line[6]) < int(line[8]):
                    tandem.add(line[0])
f_newBed.close()
print(len(tandem))


spacing100kb = set()
f_newBed = open('./data_test/sindFosEnd_bwa_q60.newform', 'r')
for line in f_newBed:
    line = line.strip('\n').split('\t')
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
f_newBed.close()
print(len(spacing100kb))


more1less20more50kb = set()
f_newBed = open('./data_test/sindFosEnd_bwa_q60.newform', 'r')
for line in f_newBed:
    line = line.strip('\n').split('\t')
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
f_newBed.close()
print(len(more1less20more50kb))
print(more1less20more50kb)


less1kb = set()
f_newBed = open('./data_test/sindFosEnd_bwa_q60.newform', 'r')
for line in f_newBed:
    line = line.strip('\n').split('\t')
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
f_newBed.close()
print(len(less1kb))


allRead = set()
rightPlace = set()
f_newBed = open('./data_test/sindFosEnd_bwa_q60.newform', 'r')
for line in f_newBed:
    line = line.strip('\n').split('\t')
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
f_newBed.close()
print(len(rightPlace))

consist = link2Chr | tandem | inverted | spacing100kb | less1kb | more1less20more50kb | rightPlace
print(consist == allRead)
