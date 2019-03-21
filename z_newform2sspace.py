f_newform = open('./data_test7/gu_yu18_mmi_q60F2048.newform', 'r')
f_tabsspace = open('./data_test7/gu_yu18_mmi_q60F2048_newform2Sspace.tab', 'w')
for line in f_newform:
    line = line.strip('\n').split('\t')
    strand1 = (int(line[1]) >> (5 - 1)) & 1
    strand2 = (int(line[3]) >> (5 - 1)) & 1
    if strand1 == 0 and strand2 == 0:
        newline = [line[2], line[5], line[6], line[4], line[7], line[8]]
        waitwrite = '\t'.join(newline) + '\n'
        f_tabsspace.write(waitwrite)
    if strand1 == 1 and strand2 == 1:
        newline = [line[2], line[6], line[5], line[4], line[8], line[7]]
        waitwrite = '\t'.join(newline) + '\n'
        f_tabsspace.write(waitwrite)
    if strand1 == 0 and strand2 == 1:
        newline = [line[2], line[5], line[6], line[4], line[8], line[7]]
        waitwrite = '\t'.join(newline) + '\n'
        f_tabsspace.write(waitwrite)
    if strand1 == 1 and strand2 == 0:
        newline = [line[2], line[6], line[5], line[4], line[7], line[8]]
        waitwrite = '\t'.join(newline) + '\n'
        f_tabsspace.write(waitwrite)
f_tabsspace.close()
f_newform.close()
