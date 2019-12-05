import os
import sys
import argparse
import gzip
import re
import random
import subprocess
import math
from subprocess import Popen, PIPE

def sample_to_alignment(fn, b, wl):
    
    n = 0
    with gzip.open(fn, 'rt') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            n = n + 1
    sample_name = fn.split(".")       
    bf = open(b,'r')
    barcode = []
    for l in bf:
        ss = []
        ss.append(l.split(':')[0])
        ss.append(l.split(':')[1][6:20])
        ss.append(l.split(':')[1][6:20])
        barcode.append(ss)

    newfq = []
    for m in barcode:
        newfq.append( lines[int(m[0])-2] )
        newfq.append( lines[int(m[0])-1] )
        newfq.append( lines[int(m[0])] )
        newfq.append( lines[int(m[0])+1] )
        
    fqMatrix = [[0 for x in range(6)] for y in range(len(barcode))]
    for l in range(len(newfq)):
        if l%4 == 0:
            fqMatrix[int(l/4)][0] = newfq[l]
        elif l%4 == 1:
            fqMatrix[int(l/4)][1] = newfq[l]
        elif l%4 == 2:
            fqMatrix[int(l/4)][2] = newfq[l]
        else:
            fqMatrix[int(l/4)][3] = newfq[l]
        fqMatrix[int(l/4)][4] = barcode[int(l/4)][1]
        fqMatrix[int(l/4)][5] = barcode[int(l/4)][2]

    whiteList = []
    w = open(wl,'r')
    for wli in w:
        whiteList.append(wli.rstrip())

    newList = []
    num = 0
    for ll in range(len(fqMatrix)):
        if fqMatrix[ll][5] in whiteList:
            newList.append(fqMatrix[ll])
    keyList = []
    d = {}
    for bm in range(len(newList)):
        if newList[bm][5] not in d.keys():
            a = []
            a.append(newList[bm][0])
            a.append(newList[bm][1])
            a.append(newList[bm][2])
            a.append(newList[bm][3])
            d[newList[bm][5]] = a
        else:
            d[newList[bm][5]].append(newList[bm][0])
            d[newList[bm][5]].append(newList[bm][1])
            d[newList[bm][5]].append(newList[bm][2])
            d[newList[bm][5]].append(newList[bm][3])

#
# grep the target regions
#   
    d_t = {}
    for i in d.keys():
        s_list = d[i]
        ns_list = []
        f = open(i + ".fa", "w")
        for counter, line in enumerate(s_list):
            if counter % 4 == 0:
                f.write(">"+i)
                f.write("\n")
            if counter % 4 == 1:
                m = re.search(i, str(line))
                start = m.start()
                subline = line[start:(len(line) - 1)]
                target = subline[0:-22]
                if len(target) != 0 or start == None:
                    ns_list.append(target)
                    f.write(target)
                    f.write("\n")
                else:
                    ns_list.append(i)
                    f.write(i)
                    f.write("\n")
        d_t[i] = ns_list
        f.close()

                    
#
# run needleall alignment
#
    os.system('sh needleall_alignment.sh' + ' ' + sample_name[0])
    
def Day_zero(bpath):
    
    table = []
    bigD = {}

    for file in os.listdir(bpath):
        if (".needleall" in file and ".needleall.fa" not in file):
            current_file = os.path.join(bpath,file)
            name = file.split(".")
            sam = open(current_file, 'r')
            count = 0
            
            for i in sam:
                count += 1
                ss = []
                if count > 2:
                    line = str(i).split("\t")
                    if line[2] not in bigD.keys():
                        ll = []
                        if line[5] != "74M" and line[5] not in ll:
                            ll.append(line[5])
                        bigD[name[0]] = ll
                    else:
                        if line[5] != "74M" and line[5] not in bigD[name[0]]:
                            bigD[name[0]].append(line[5])
    return bigD

#
# recombination filter
#

def cigar_dic(cigar):
    start = 0
    #a = '1D10M1I10M1D20M'
    listC = []
    for num1, i_or_d, num2, m in re.findall('(\d+)([ID])(\d+)?([A-Za-z])?', cigar):
        print(num1, i_or_d, start)
        line = []
        line.append(num1)
        line.append(i_or_d)
        line.append(start)
        listC.append(line)
        if num1:
            start += int(num1)
        if num2:
            start += int(num2)
    return listC
    
    
def main(sample, day0, barcodes, barcodes_day0, whitelist, dNumber):
    
#get sample alignment results
    sample_to_alignment(sample, barcodes, whitelist)
#get Day0 alignment results
    sample_to_alignment(day0, barcodes_day0, whitelist)
    sample_name = sample.split(".")
    
#
# Day 0 filter
# 
    spath = sample
    bpath = day0
    d_filter = {}
    A = blacklist(bpath)
    for file in os.listdir(spath):
        if (".needleall" in file and ".needleall.fa" not in file):
            current_file = os.path.join(spath,file)
            nm = file.split(".")
            design = ''
            sam = open(current_file, 'r')
            #fw = open(file, 'w')
            count = 0
            df_list = []
            unedit_count = 0           
            for i in sam:
                ss = []
                design = ""
                if count < 2:
                    fw.write(i)
                if count >= 2:
                    line = str(i).strip("\n").split("\t")
                    design = line[2]
                    cigar = line[5]
                    if nm[0] in A.keys():
                        bl = A[nm[0]]
                        if cigar not in bl:
                            df_list.append(i)
                    else:
                        df_list.append(i)
                count += 1
            d_filter[nm[0]] = df_list

#
# downsample
#
    dsp_list = {}
    for i in d_filter.keys():
        if len(d_filter[i]) >= dNumber:
            m_list = random.sample(d_filter[i], dNumber)
            dsp_list[i] = m_list

#
# entropy calculation
#
    table = []
    bigD = []
    entropyD = {}
    for dline in dsp_list.keys():
        clist = []
        d_p = {}
        unedit_count = 0
        design = ""
        for count, i in enumerate(dsp_list[dline]):
            dlines = str(i).strip("\n").split("\t")
            cigar = dlines[5]
            if cigar == "74M":
                unedit_count +=1
            design = dlines[2]
            if cigar not in d_p.keys():
                d_p[cigar] = 1
            else:
                d_p[cigar] += 1 
        if count > 0:
            total_eff = 1 - (unedit_count / count)
        else:
            total_eff = 0
        table.append(design)
        table.append(d_p)
        bigD.append(table)
        nd = {}
        name = design
        entropy = 0
        for j in d_p.keys():
            nd[j] = d_p[j]/count
            if nd[j] > 0 and j != 'None':
                    entropy += - nd[j]*math.log(nd[j], 2)
        lef = []
        lef.append(entropy)
        lef.append(total_eff)
        entropyD[design] = lef
        sorted_x = sorted(entropyD.items(), key=lambda kv: kv[1], reverse=True)
    
    with open(sample_name[0] + ".txt", 'w') as fe:
        for item in sorted_x:
            fe.write(str(item[0]) + "\t" + str(item[1][0]) + "\t" + str(item[1][1]))
            fe.write("\n")
    fe.close()

if __name__ == '__main__':
    main()
