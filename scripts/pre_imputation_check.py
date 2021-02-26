import logging


ATGC = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
atgc = {'A': 1, 'T': 1, 'G': 0, 'C': 0}


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


def get_id(a1_a2, id_element, chr_, pos_, ref_, alt_):
    # first, try to use rs
    id_ = id_element.split(':')[0]
    if id_ in a1_a2:
        return id_
    # then, try chr:pos
    id_ = "{}:{}".format(chr_, pos_)
    if id_ in a1_a2:
        return id_
    # try chr:pos:ref:alt
    id_ = "{}:{}:{}:{}".format(chr_, pos_, ref_, alt_)
    if id_ in a1_a2:
        return id_
    # try chr:pos:alt:ref
    id_ = "{}:{}:{}:{}".format(chr_, pos_, alt_, ref_)
    if id_ in a1_a2:
        return id_
    return None


def pre_imputation_check(params, reference, fn='./data/hapmap20', rs=True):
    """Given a plink bim file
    1. Remove SNPs can not be matched
    2. Flip SNPs that match after flip
    3. Swap SNPs that match after swap or flip&swap
    4. Update chromosomes and positions for these SNPs
    """
    # read in bim first
    prefix = params
    fn_new = prefix  # + '_filter.bim'
    # filter_bim(fn + '.bim', fn_new)
    a1_a2 = {}
    for i in open(fn_new):
        items = i.split()
        # id: (ref, alt)
        a1_a2[items[1]] = (items[-2], items[-1])

    # files to update bim
    chr_file = prefix + '.chr'
    pos_file = prefix + '.pos'
    force_allele_file = prefix + '.force_allele'
    flip_file = prefix + '.flip'
    with open(chr_file, 'w') as w1, open(pos_file, 'w') as w2, open(force_allele_file, 'w') as w3, open(flip_file, 'w') as w4:

        # write the files to be used
        in_ref = 0
        exclude = 0
        strand_flip = 0
        swap = 0
        flip_swap = 0

        for i in open(reference):
            items = i.split()
            chr_, pos_, ref_, alt_ = items[0], items[1], items[3], items[4]
            # first, try to use rs
            id_ = get_id(a1_a2, items[2], chr_, pos_, ref_, alt_)

            if id_ is not None:
                a1, a2 = a1_a2[id_]
                matching = match_ref(a1, a2, ref_, alt_)
                in_ref += 1
                if matching == 4:
                    exclude += 1
                else:
                    w1.write("{}\t{}\n".format(id_, chr_))
                    w2.write("{}\t{}\n".format(id_, pos_))
                    # set alt as A1, because recode vcf will code A1 as alt later
                    w3.write("{}\t{}\n".format(id_, alt_))
                    if matching == 2:
                        w4.write(id_ + "\n")
                        strand_flip += 1
                    elif matching == 1:
                        swap += 1
                    elif matching == 3:
                        w4.write(id_ + "\n")
                        flip_swap += 1

    logging.info("{} ids in {}.bim, {} can be found in reference.".format(len(a1_a2), fn, in_ref))
    logging.info("Exclude: {} Keep: {}".format(exclude, in_ref - exclude))
    logging.info("Total flip: {}.".format(strand_flip))
    logging.info("Total swap: {}.".format(swap))
    logging.info("Total flip&swap: {}.".format(flip_swap))
    return flip_file, force_allele_file, chr_file, pos_file


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    pre_imputation_check(snakemake.input[0], snakemake.params[0])
