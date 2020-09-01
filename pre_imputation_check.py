import os,sys

ATGC = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
atgc = {'A':1, 'T':1, 'G':0, 'C':0}

def is_ambiguous(a1, a2):
    """if SNP is ambiguous or not""" 
    if atgc[a1] == atgc[a2]:
        return True
    return False

def match_ref(a1, a2, ref, alt):
    """Match a1, a2 with ref and alt
    None: match, 1: swap, 2: flip, 3: flip&swap, 4: exclude
    """
    if a1 == ref and a2 == alt:
        return
    elif a1 == alt and a2 == ref:
        return 1
    else:
        ambiguous = is_ambiguous(a1, a2)
        if ambiguous:
            return 4
        else:
            a1_f, a2_f = [ATGC[i] for i in [a1, a2]]
            if a1_f == ref and a2_f == alt:
                return 2
            elif a1_f == alt and a2_f == ref:
                return 3
            else:
                return 4

SITE_1000GENOME = "/media/1000genome/allele_info/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.only_rs.biallelic.tab"

def pre_imputation_check(params, fn='./data/hapmap20', rs=True, reference=SITE_1000GENOME):
    """Given a plink bim file
    1. Remove SNPs can not be matched
    2. Flip SNPs that match after flip
    3. Swap SNPs that match after swap or flip&swap
    4. Update chromosomes and positions for these SNPs
    """
    # read in bim first
    prefix = params
    fn_new = prefix #+ '_filter.bim'
    #filter_bim(fn + '.bim', fn_new)
    a1_a2 = {}
    for i in open(fn_new):
        items = i.split()
        a1_a2[items[1]] = (items[-2], items[-1])

    # files to update bim
    chr_file = prefix + '.chr'
    pos_file = prefix + '.pos'
    force_allele_file = prefix + '.force_allele'
    flip_file = prefix + '.flip'
    w1 = open(chr_file, 'w')
    w2 = open(pos_file, 'w')
    w3 = open(force_allele_file, 'w')
    w4 = open(flip_file, 'w')

    # write the files to be used
    in_ref = 0
    exclude = 0
    strand_flip = 0
    swap = 0
    flip_swap = 0

    for i in open(reference):
        items = i.split()
        chr_, pos_ = items[0], items[1]
        if rs:
            id_ = items[2].split(':')[0]
        else:
            id_ = "{}:{}".format(chr_, pos_)

        if id_ in a1_a2:
            a1, a2 = a1_a2[id_]
            ref, alt = items[3], items[4]
            matching = match_ref(a1, a2, ref, alt)
            in_ref += 1
            if matching == 4:
                exclude += 1
            else:
                w1.write("{}\t{}\n".format(id_, chr_))
                w2.write("{}\t{}\n".format(id_, pos_))
                # set alt as A1, because recode vcf will code A1 as alt later
                w3.write("{}\t{}\n".format(id_, alt))
                if matching == 2:
                    w4.write(id_ + "\n")
                    strand_flip += 1
                elif matching == 1:
                    swap += 1
                elif matching == 3:
                    w4.write(id_ + "\n")
                    flip_swap += 1
    w1.close()
    w2.close()
    w3.close()
    w4.close()
    print("{} ids in {}.bim, {} can be found in reference.".format(len(a1_a2), fn, in_ref))
    print("Exclude: {} Keep: {}".format(exclude, in_ref - exclude))
    print("Total flip: {}.".format(strand_flip))
    print("Total swap: {}.".format(swap))
    print("Total flip&swap: {}.".format(flip_swap))
    return flip_file, force_allele_file, chr_file, pos_file


if __name__ == "__main__":
    pre_imputation_check(sys.argv[1])
