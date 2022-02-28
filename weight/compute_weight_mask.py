"""
TODO: ...
"""
import argparse

from ibd_segments_weigher import IBDSegmentsWeigher


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ibd-segments-file', help='TODO')
    parser.add_argument('--output-mask-file', help='Result weight mask file, stored in JSON format')
    parser.add_argument('--output-mask-plot', help='Result weight mask plot')
    args = parser.parse_args()

    weigher = IBDSegmentsWeigher.from_ibd_segments_file(args.ibd_segments_file)
    weigher.to_json_mask_file(args.output_mask_file)
    weigher.visualize(args.output_mask_plot)
