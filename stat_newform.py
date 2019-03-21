import numpy as np


n = 0
inverted = []
spaceing_error = []
length_list = []
f_newform = open('./data_test7/gu_yu18_mmi_q60F2048.newform', 'r')
for line in f_newform:
    n += 1
    line = line.strip('\n').split('\t')
    if line[2] == line[4]:
        strand1 = (int(line[1]) >> (5 - 1)) & 1
        strand2 = (int(line[3]) >> (5 - 1)) & 1
        if strand1 != strand2:
            inverted.append(line[0])
        else:
            if strand1 == 0:
                length = int(line[8]) - int(line[5])
                if length < 26600 or length > 49400:
                    spaceing_error.append((line[0], length))
                else:
                    length_list.append(length)
            else:
                length = int(line[6]) - int(line[7])
                if length < 26600 or length > 49400:
                    spaceing_error.append((line[0], length))
                else:
                    length_list.append(length)
f_newform.close()
print(n)
print(inverted)
print(spaceing_error)
print(round(len(length_list)/n, 4)*100)
print(round(np.array(length_list).mean(), 2))
print(round(np.array(length_list).std(), 2))
