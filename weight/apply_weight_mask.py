"""
TODO: ...
"""
import argparse

from ibd_segments_weigher import IBDSegmentsWeigher


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-ibd-segments-file', help='TODO')
    parser.add_argument('--mask-file', help='TODO')
    parser.add_argument('--output-ibd-segments-file', help='TODO')
    args = parser.parse_args()

    weigher = IBDSegmentsWeigher.from_json_mask_file(args.mask_file)
    weigher.apply_mask(args.input_ibd_segments_file, args.output_ibd_segments_file)
