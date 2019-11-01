import numpy as np
import random


seq = 'GTCCATGTATCGTTTTATAATACTTTTTTAAGCACTTTCTTTAATAAATATTCCATTGAA' \
      'GTACTGTTACTGAAATGAGATGAACTGTTCAGAATGTAGAAATGGCGCCAGAAATCAATC'

def introduce_indelS(seq):
    seq = list(seq)
    length = len(seq)
    # generate errors
    insertErrors = abs(np.random.normal(0.00939, 0.0165))
    deleteErrors = abs(np.random.normal(0.00957, 0.01084))
    subsituteErrors = abs(np.random.normal(0.0374, 0.01585))
    insertNum = round(length * insertErrors)
    deleteNum = round(length * deleteErrors)
    subsituteNum = round(length * subsituteErrors)
    errorNum = insertNum + deleteNum + subsituteNum
    position = random.sample(range(length), errorNum)
    insertPosition = sorted(position[:insertNum])
    deletePosition = position[insertNum:insertNum + deleteNum]
    subsitutePosition = position[insertNum + deleteNum:]
    # introduce errors
    for i in subsitutePosition:
        subs_base = np.random.choice(['A', 'T', 'C', 'G'], p=[0.25, 0.25, 0.25, 0.25])
        seq[i] = subs_base
    for i in deletePosition:
        seq[i] = 'N'
    if insertPosition != []:
        afters = []
        afters += seq[0:insertPosition[0] + 1]
        inse_base = np.random.choice(['A', 'T', 'C', 'G'], p=[0.25, 0.25, 0.25, 0.25])
        afters += [inse_base]
        for i in range(len(insertPosition) - 1):
            afters += seq[insertPosition[i] + 1:insertPosition[i + 1] + 1]
            inse_base = np.random.choice(['A', 'T', 'C', 'G'], p=[0.25, 0.25, 0.25, 0.25])
            afters += [inse_base]
        afters += seq[insertPosition[-1] + 1:]
        finals = [i for i in afters if i != 'N']
    else:
        finals = [i for i in seq if i != 'N']
    # print(insertNum, deleteNum, subsituteNum)
    # print(insertPosition, deletePosition, subsitutePosition)
    return ''.join(finals)

print(introduce_indelS(seq))
print(seq)
