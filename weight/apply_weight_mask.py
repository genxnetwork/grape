import argparse

from ibd_segments_weigher import IBDSegmentsWeigher


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Apply weight mask stored in JSON to the IBIS IBD segments file (.seg)'
    )
    parser.add_argument('--input-ibd-segments-file', help='File with IBD segments detected with IBIS (.seg file)')
    parser.add_argument('--mask-file', help='Weight mask, computed with IBDSegmentsWeigher')
    parser.add_argument('--output-ibd-segments-file', help='Weighed IBD segments file in the IBIS format (.seg file)')
    args = parser.parse_args()

    weigher = IBDSegmentsWeigher.from_json_mask_file(args.mask_file)
    weigher.apply_mask(args.input_ibd_segments_file, args.output_ibd_segments_file)
