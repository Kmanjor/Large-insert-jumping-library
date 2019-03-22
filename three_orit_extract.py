import pandas as pd
import pprint
from Bio import SeqIO


total_df = pd.read_table('gu2B_ccs_minPase2_three_feature.blastn',
                         names=['qseqid', 'qlen', 'sseqid', 'slen', 'qstart', 'qend', 'pident', 'sstrand', 'sstart', 'send'])
total_dict = {i: total_df.loc[total_df.qseqid == i] for i in set(total_df.qseqid)}
seq_dict = SeqIO.to_dict(SeqIO.parse('m54045_180307_194900.ccs.fasta', 'fasta'))

left_list = []
right_list = []

the_true = []
for i in total_dict:
    df = total_dict[i]
    if len(df) == 3:
        amp_tag = df.loc[df.sseqid == 'amp_tag']
        IsceI_1 = df.loc[df.sseqid == 'IsceI_1']
        IsceI_2 = df.loc[df.sseqid == 'IsceI_2']
        if len(amp_tag) == 1 and len(IsceI_1) == 1 and len(IsceI_2) == 1:
            # print(str(IsceI_1.sstrand.values))
            # print(str(IsceI_2.sstrand.values))
            if IsceI_1.sstrand.values == IsceI_2.sstrand.values:
                if IsceI_1.sstrand.values == ['plus']:
                    I1S = int(IsceI_1.qstart.values[0])
                    I2E = int(IsceI_2.qend.values[0])
                    SeqL = int(IsceI_1.qlen.values[0])
                    temp = SeqL - I2E
                    # print('p', I1S, I2E, SeqL)
                    if I1S < 21 and temp < 21:
                        the_true.append(IsceI_1.qseqid.values[0])
                        amp_start = int(amp_tag.qstart.values[0])
                        amp_end = int(amp_tag.qend.values[0])
                        I1E = int(IsceI_1.qend.values[0])
                        I2S = int(IsceI_2.qstart.values[0])
                        record = seq_dict[IsceI_1.qseqid.values[0]]
                        if amp_start - I1E > 50 and I2S - amp_end > 50:
                            fosmid_left = record[I1E:amp_start-1]
                            fosmid_right = record[amp_end:I2S-1]
                            fosmid_left.id, fosmid_left.name, fosmid_left.description = fosmid_left.id + '/1', '', ''
                            fosmid_right.id, fosmid_right.name, fosmid_right.description = fosmid_right.id + '/2', '', ''
                            left_list.append(fosmid_left)
                            right_list.append(fosmid_right)
                if IsceI_1.sstrand.values == ['minus']:
                    I2S = int(IsceI_2.qstart.values[0])
                    I1E = int(IsceI_1.qend.values[0])
                    SeqL = int(IsceI_2.qlen.values[0])
                    temp = SeqL - I1E
                    # print('m', I2S, I1E, SeqL)
                    if I2S < 21 and temp < 21:
                        the_true.append(IsceI_1.qseqid.values[0])
                        amp_start = int(amp_tag.qstart.values[0])
                        amp_end = int(amp_tag.qend.values[0])
                        I2E = int(IsceI_2.qend.values[0])
                        I1S = int(IsceI_1.qstart.values[0])
                        record = seq_dict[IsceI_2.qseqid.values[0]]
                        if amp_start - I2E > 50 and I1S - amp_end > 50:
                            fosmid_left = record[I2E:amp_start-1]
                            fosmid_right = record[amp_end:I1S-1]
                            fosmid_left.id, fosmid_left.name, fosmid_left.description = fosmid_left.id + '/1', '', ''
                            fosmid_right.id, fosmid_right.name, fosmid_right.description = fosmid_right.id + '/2', '', ''
                            left_list.append(fosmid_left)
                            right_list.append(fosmid_right)

SeqIO.write(left_list, '3orit_gu2B-1.fa', 'fasta')
SeqIO.write(right_list, '3orit_gu2B-2.fa', 'fasta')

print(len(the_true))


