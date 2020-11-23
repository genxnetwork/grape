#!/usr/bin/python
import argparse
import os.path


def mergeIntervals(intervals):
    # Given list/set 'intervals', merge any overlaps
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def mergeTwoLists(list1,list2):
    tmp = list1 + list2
    tmp.sort()
    return tmp


def checkInKeysAndAddToList(dict,ind1,ind2,chr,item):
    # for dict, check whether dict[ind1][ind2][chr] exists, and if not, add it, and then add item
    if not ind1 in dict.keys():
        dict[ind1] = {}
    if chr != '':
        if not ind2 in dict[ind1].keys():
            dict[ind1][ind2] = {}
        if not chr in dict[ind1][ind2].keys():
            dict[ind1][ind2][chr] = []
        if item != '':
            dict[ind1][ind2][chr].append(item)
    else: #used primarily for IBD_lengths creations
        if not ind2 in dict[ind1].keys():
            dict[ind1][ind2] = [0,0]


def overlap(self, other):
    # whether the interval self overlaps interval other, and whether their haplotypes are opposite

    if self[1] >= other[0] and self[0] <= other[1] and self[2][0] != other[2][0] and self[2][1] != other[2][1]: #ensure different homologues
        return [max(self[0],other[0]),min(self[1],other[1])]
    else:
        return 0


if __name__ == '__main__':

    parser=argparse.ArgumentParser(
        description='''merge_refined_ibd.py, 3-27-2018''',
        epilog="""Requires the output of three Refined IBD runs (using different seeds on same datasets) and a map file of all markers used in the Refined IBD runs.\nCreates output file *.ibd12 which contains columns sample1, sample2, IBD1 length/proportion, IBD2 length/proportion, and optionally creates a *.seg which contains sample1, sample2, chromosome, IBD type, region start (in cM), region end (in CM).""")
    parser.add_argument('-o', type=str, nargs=1, default=['out'], help='Output file prefix', metavar='out')
    parser.add_argument('-f', type=str, nargs=3, required=True, help='Three Refined IBD files', metavar=('file1','file2','file3'))
    parser.add_argument('-m', type=str, nargs=1, required=True, help='Map file (PLINK format), non-zero cM column required', metavar='file.map')
    parser.add_argument('-t', type=int, nargs=1, default=[0], help='Whether to print IBD1 & IBD2 lengths in cM (=0), or IBD1 & IBD2 proportions (=1)',metavar='0/1')
    parser.add_argument('-s', type=int, nargs=1, default=[0], help='Whether to print *.seg file which contains the discovered pairwise IBD regions and whether they are IBD1 or IBD2', metavar='0/1')

    args=parser.parse_args() #get user input

    print('\n')
    if os.path.isfile(args.o[0]+'.ibd12') or os.path.isfile(args.o[0]+'.seg'):
      k=1
      new = args.o[0] + '-' + str(k)
      while os.path.isfile(new+'.ibd12') or os.path.isfile(new+'.seg'):
        k = k + 1
        new = args.o[0] + '-' + str(k)
      print('Output file with prefix already exists, changing prefix to '+new+'\n')
      args.o[0] = new

    print("Input files: "+args.f[0]+' '+args.f[1]+' '+args.f[2])
    print("Map file: "+args.m[0])
    print("Output files: "+args.o[0]+'.ibd12, '+args.o[0]+'.seg')
    if args.t[0] == 1:
      print("Printing IBD1 & IBD2 proportions")
    else:
      print("Printing IBD1 & IBD2 lengths in cM")
    if args.s[0] == 1:
      print("Printing both *.ibd12 and *.seg files\n")
    else:
      print("Printing only *.ibd2 file\n")


    # Read in map file
    loc = {}
    for chr in range(1,23):
      loc[chr] = {}

    loc_data = open(args.m[0],'r')

    chr_ends = {}
    for chr in range(1,23):
        chr_ends[chr] = [1e100,0]

    total_genome = 0
    for line in loc_data:
      l = str.split(line.rstrip())
      chr = int(l[0])
      loc[chr][l[3]] = float(l[2])
      if loc[chr][l[3]] < chr_ends[chr][0]:
          chr_ends[chr][0] = loc[chr][l[3]]
      elif loc[chr][l[3]] > chr_ends[chr][1]:
          chr_ends[chr][1] = loc[chr][l[3]]

    loc_data.close()

    # get chromosome information
    total_genome = 0
    for chr in range(1,23):
      total_genome = total_genome + chr_ends[chr][1] - chr_ends[chr][0]


    print("Total genome length: "+str(total_genome)+"\n")


    IBD1_all = {} #collection of IBD1 segments from each of the three runs
    IBD2_all = {} #collection of IBD2 segments from each of the three runs
    IBDseg_all = {} #all IBD segments in three runs


    for run in range(0, 3):
        IBD2 = {}
        IBDseg = {}  # pairwise IBD segments: IBDseg[ind1][ind2][chromosome] where value of ind1 < value of ind2
        IBDinput = open(args.f[run], 'r')  # open file containing Refined IBD output
        for line in IBDinput:
            fields = str.split(line.rstrip())
            if fields[0] < fields[2]:
                ind1 = fields[0]
                ind2 = fields[2]
                h1 = int(fields[1]) - 1 #store homologue number for ind1 as 0/1
                h2 = int(fields[3]) - 1 #store homologue number for ind2 as 0/1
            else:
                ind1 = fields[2]
                ind2 = fields[0]
                h1 = int(fields[3]) - 1
                h2 = int(fields[1]) - 1
            if 'chr' in fields[4]:
                chr = int(str.split(fields[4], 'chr')[1]) #in case of older/new version differences
            else:
                chr = int(fields[4])
            checkInKeysAndAddToList(IBDseg,ind1,ind2,chr,[loc[chr][fields[5]],loc[chr][fields[6]],str(h1)+str(h2)]) #add IBD segment as Interval with homologue info (interval_start, interval end, h1h2) to set
        IBDinput.close()

        # Find IBD2 regions, add to IBD2 dict
        for ind1 in IBDseg.keys():
            for ind2 in IBDseg[ind1].keys():
                for chr in IBDseg[ind1][ind2].keys():
                    if len(IBDseg[ind1][ind2][chr]) > 1:  # if more than one segment on chromosome 'chr' between ind1 and ind2, check for IBD2
                        maxval = len(IBDseg[ind1][ind2][chr])
                        IBDseg[ind1][ind2][chr].sort()
                        IBDlist = IBDseg[ind1][ind2][chr]
                        k = 0
                        while k < maxval - 1:
                            kk = k + 1
                            while IBDlist[k][1] > IBDlist[kk][0]:  # while first segment ends after a following segment starts
                                ov = overlap(IBDlist[k],IBDlist[kk]) #check whether segment k overlaps segment kk, returns that overlap
                                if ov != 0:
                                    checkInKeysAndAddToList(IBD2, ind1, ind2, chr, ov)  # prepare IBD2 dict, add 'overlap' IBD2 segment
                                kk = kk + 1
                                if kk >= maxval:
                                    break
                            k = k + 1

        # add this run's IBD segments to master list of IBD segments from all runs (IBDseg_all)
        for ind1 in IBDseg.keys():
            for ind2 in IBDseg[ind1].keys():
                for chr in IBDseg[ind1][ind2].keys():
                    checkInKeysAndAddToList(IBDseg_all,ind1,ind2,chr,'')
                    IBDseg_all[ind1][ind2][chr] += IBDseg[ind1][ind2][chr]

        # add IBD2 segments from this run to IBD2 master list of IBD2 from all runs (IBD2_all)
        for ind1 in IBD2.keys():
            for ind2 in IBD2[ind1].keys():
                for chr in IBD2[ind1][ind2].keys():
                    checkInKeysAndAddToList(IBD2_all,ind1,ind2,chr,'')
                    IBD2_all[ind1][ind2][chr] += IBD2[ind1][ind2][chr]




    # Done collecting IBD info from all three runs
    # Now greedily calculate IBD2 regions by looking at IBD2 in all three runs

    for ind1 in IBD2_all.keys():
        for ind2 in IBD2_all[ind1].keys():
            for chr in IBD2_all[ind1][ind2].keys():
                IBD2_all[ind1][ind2][chr] = mergeIntervals(IBD2_all[ind1][ind2][chr])  # greedily merge our IBD2 intervals collected from the three runs, outputs a list

    IBD_lengths = {}

    for ind1 in IBDseg_all.keys():
        for ind2 in IBDseg_all[ind1].keys():
            for chr in IBDseg_all[ind1][ind2].keys():
                checkInKeysAndAddToList(IBD1_all, ind1, ind2, chr, '') #prepare IBD1_all dict
                if len(IBDseg_all[ind1][ind2][chr]) > 1: #IBDseg_all[ind1][ind2][chr] is not yet merged at this point
                    # find IBD1 by finding all regions that are IBD and removing those which are IBD2
                    IBDseg_all[ind1][ind2][chr].sort()
                    merged = mergeIntervals(IBDseg_all[ind1][ind2][chr])  # merge all IBD intervals collected from the three runs, outputs a list

                    # check for overlap between 'merged' IBD segments and IBD2 segments
                    # don't need to worry about homologue info at this point
                    if ind1 in IBD2_all.keys() and ind2 in IBD2_all[ind1].keys() and chr in IBD2_all[ind1][ind2].keys() and len(IBD2_all[ind1][ind2][chr]) > 0:
                        k_init = 0
                        for i in range(0, len(merged)):
                            if k_init >= len(IBD2_all[ind1][ind2][chr]):
                                break
                            seg = merged[i] #IBD segment; need to determine which part(s) of segments are IBD1 or IBD2
                            k = k_init
                            while k < len(IBD2_all[ind1][ind2][chr]) and IBD2_all[ind1][ind2][chr][k][1] < seg[0]:
                                k = k + 1
                            # while IBD2 segment is ending before merged segment starts, move to next IBD2 segment
                            if k >= len(IBD2_all[ind1][ind2][chr]):
                                IBD1_all[ind1][ind2][chr].append([seg[0], seg[1]])
                                break
                            while k < len(IBD2_all[ind1][ind2][chr]):
                                # check if the IBD segment ends before the next IBD2 segment
                                if seg[1] < IBD2_all[ind1][ind2][chr][k][0]:
                                    if k > 0:
                                        # we already added seg to IBD1 list if a previous IBD2 segment overlapped it
                                        if IBD2_all[ind1][ind2][chr][k - 1][1] < seg[0]: #end of previous IBD2 segment occurs before start of this IBD segment
                                            IBD1_all[ind1][ind2][chr].append([seg[0], seg[1]])
                                    else:
                                        IBD1_all[ind1][ind2][chr].append([seg[0], seg[1]])
                                    break
                                if k < len(IBD2_all[ind1][ind2][chr]) - 1:
                                    # if this IBD2_all and next IBD2_all end before merged segment, add space between this IBD2_all segment and next IBD2_all segment
                                    if IBD2_all[ind1][ind2][chr][k][0] >= seg[0] and IBD2_all[ind1][ind2][chr][k][1] < seg[1]:
                                        if IBD2_all[ind1][ind2][chr][k + 1][1] < seg[1]:
                                            IBD1_all[ind1][ind2][chr].append([IBD2_all[ind1][ind2][chr][k][1], IBD2_all[ind1][ind2][chr][k + 1][0]])
                                        # if this IBD2_all segment ends before merged segment but next IBD2_all segment does not
                                        else:
                                            IBD1_all[ind1][ind2][chr].append([IBD2_all[ind1][ind2][chr][k][1], min(seg[1], IBD2_all[ind1][ind2][chr][k + 1][0])])
                                else:
                                    # if this IBD2_all segment ends before merged segment but next IBD2_all segment does not
                                    if IBD2_all[ind1][ind2][chr][k][0] >= seg[0] and IBD2_all[ind1][ind2][chr][k][1] < seg[1]:
                                        IBD1_all[ind1][ind2][chr].append([IBD2_all[ind1][ind2][chr][k][1], seg[1]])
                                if k > 0:
                                    # if this IBD2_all starts after merged segment starts and previous IBD2_all segment did not
                                    if IBD2_all[ind1][ind2][chr][k][0] > seg[0] and IBD2_all[ind1][ind2][chr][k - 1][1] < seg[0]:
                                        IBD1_all[ind1][ind2][chr].append([seg[0], IBD2_all[ind1][ind2][chr][k][0]])
                                    # if this IBD2_all segment starts after merged segment and previous IBD2_all segment ends before merged segment
                                    if seg[1] < IBD2_all[ind1][ind2][chr][k][0] and IBD2_all[ind1][ind2][chr][k - 1][1] < seg[0]:
                                        IBD1_all[ind1][ind2][chr].append([seg[0], seg[1]])
                                # if this IBD2_all starts after merged segment starts and there is no previous IBD2_all segment
                                else:
                                    if IBD2_all[ind1][ind2][chr][k][0] > seg[0]:
                                        IBD1_all[ind1][ind2][chr].append([seg[0], IBD2_all[ind1][ind2][chr][k][0]])
                                if IBD2_all[ind1][ind2][chr][k][0] > seg[1]:
                                    break
                                k = k + 1
                            k_init = k
                        if i < len(merged):
                            for ii in range(i,len(merged)):
                                start = max(merged[ii][0], IBD2_all[ind1][ind2][chr][-1][1])
                                end = merged[ii][1]
                                if start != end:
                                    if not [start,end] in IBD1_all[ind1][ind2][chr]:
                                        if len(IBD1_all[ind1][ind2][chr]):
                                            IBD1_all[ind1][ind2][chr].append([start, end])
                                        else:
                                            IBD1_all[ind1][ind2][chr].append([start,end])
                    else: #if there's no IBD2 segments between these two individuals and on chromosome chr, all the IBD segments must be IBD1
                        IBD1_all[ind1][ind2][chr] += merged
                else:  # only one segment on chromosome, can assume it's IBD1
                    IBD1_all[ind1][ind2][chr] = IBDseg_all[ind1][ind2][chr]

                # Calculate total lengths of IBD1 & IBD2
                checkInKeysAndAddToList(IBD_lengths, ind1, ind2, '', '') #initiate dict list to [0,0] for IBD1, IBD2
                for ibdseg in IBD1_all[ind1][ind2][chr]:
                    #print("Adding "+str(ibdseg[0])+' '+str(ibdseg[1])+'\n')
                    IBD_lengths[ind1][ind2][0] += (ibdseg[1] - ibdseg[0]) #add segment length

                if ind1 in IBD2_all.keys() and ind2 in IBD2_all[ind1].keys() and chr in IBD2_all[ind1][ind2].keys():
                    for ibdseg in IBD2_all[ind1][ind2][chr]:
                        IBD_lengths[ind1][ind2][1] += (ibdseg[1] - ibdseg[0]) #add segment length

            if args.t[0] == 1:
                IBD_lengths[ind1][ind2][0] /= total_genome
                IBD_lengths[ind1][ind2][1] /= total_genome

    outfile_IBD = open(args.o[0] + '.ibd12', 'w')
    for ind1 in IBD_lengths.keys():
        for ind2 in IBD_lengths[ind1].keys():
            # outputs [ ind1    ind2    IBD1_between_ind1_ind2  IBD2_between_ind1_ind2 ]
            outfile_IBD.write(ind1 + '\t' + ind2 + '\t' + str(IBD_lengths[ind1][ind2][0]) + '\t' + str(IBD_lengths[ind1][ind2][1]) + '\n')

    outfile_IBD.close()

    if args.s[0] == 1:
        outfile_IBD = open(args.o[0] + '.seg', 'w')
        for ind1 in IBD1_all.keys():
            for ind2 in IBD1_all[ind1].keys():
                for chr in IBD1_all[ind1][ind2].keys():
                    for ibdseg in IBD1_all[ind1][ind2][chr]:
                        if ibdseg[0] != ibdseg[1]:
                            outfile_IBD.write(ind1 + '\t' + ind2 + '\t' + str(chr) + '\tIBD1\t' + str(ibdseg[0]) + '\t' + str(ibdseg[1]) + '\n')

        for ind1 in IBD2_all.keys():
            for ind2 in IBD2_all[ind1].keys():
                for chr in IBD2_all[ind1][ind2].keys():
                    for ibdseg in IBD2_all[ind1][ind2][chr]:
                        if ibdseg[0] != ibdseg[1]:
                            outfile_IBD.write(ind1 + '\t' + ind2 + '\t' + str(chr) + '\tIBD2\t' + str(ibdseg[0]) + '\t' + str(ibdseg[1]) + '\n')

        outfile_IBD.close()
