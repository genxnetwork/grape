import csv
import json
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from collections import defaultdict
from sklearn.covariance import MinCovDet
from matplotlib import colors
from matplotlib.collections import BrokenBarHCollection
from matplotlib.patches import Rectangle


class IBDSegmentsWeigher:
    """
    TODO: ...
    """

    """
    TODO: ...
    """
    SUPPORT_FRACTION = 0.8

    """
    TODO: ...
    """
    ERSA_EXPECTED_RECOMBINATIONS_NUMBER = 35.2548101

    """
    TODO: ...
    """
    ERSA_MASKED_REGIONS = {
        1: [(150.31366473433604, 168.04559253633377)],
        2: [
            (114.77930315192536, 121.68311587384397),
            (153.57547065138377, 162.4860072200155),
            (207.60503014697173, 212.22642675353882)
        ],
        8: [(26.230374556132432, 35.022394672355695)],
        9: [(63.7961474341056, 83.1742089558682)],
        10: [(67.09773230002142, 74.14102890625604)],
        15: [
            (0.0047557671149254, 31.947180711331043),
            (37.913057336700184, 44.024287153230425)
        ],
        16: [(38.70343480460814, 44.576002743051056)],
        17: [
            (119.40062554029079, 124.78083800766407),
            (87.20365693007463, 93.12865726554318)
        ],
        21: [(3.63045785673081, 10.70749201500928)],
        22: [(0.0051839685346697995, 23.021402055048853)]
    }

    @staticmethod
    def overlap_segments_with_window(
        window: tuple[float, float],
        segments: tuple[np.ndarray, np.ndarray]
    ) -> int:
        return np.sum(np.minimum(window[1], segments[1]) - np.maximum(window[0], segments[0]))

    @staticmethod
    def load_segments(path: str) -> dict[dict[str, np.ndarray]]:
        """
        TODO:

        :param path: path to the IBD segments file
        """

        segments = defaultdict(lambda: {'start': [], 'end': []})
        with open(path) as segments_file:
            reader = csv.reader(segments_file, delimiter='\t')
            for row in reader:
                chrom, genetic_start_pos, genetic_end_pos = row[2], float(row[6]), float(row[7])
                segments[chrom]['start'].append(genetic_start_pos)
                segments[chrom]['end'].append(genetic_end_pos)

        segments = dict(segments)
        for chrom in segments:
            segments[chrom]['start'] = np.array(segments[chrom]['start'])
            segments[chrom]['end'] = np.array(segments[chrom]['end'])

        return segments

    @staticmethod
    def compute_overlaps(
        segments: dict[dict[str, np.ndarray]]
    ) -> dict[str, tuple[float, float, float]]:
        """
        TODO: ...
        """

        overlaps = defaultdict(lambda: [])
        width = 1 # We use 1cM as the window size
        for chrom in segments:
            # Define chromosome windows
            start_position = segments[chrom]['start'].min()
            end_position = segments[chrom]['end'].max()

            window_starts = np.arange(start_position, end_position, width) # last_window_end is not included
            window_ends = window_starts + width
            window_ends[-1] = end_position

            # Compute overlap length for each chromosome window
            for window_start, window_end in zip(window_starts, window_ends):
                overlapping = ~((window_start > segments[chrom]['end']) | (window_end < segments[chrom]['start']))

                if sum(overlapping) == 0:
                    overlaps[chrom].append((window_start, window_end, 0))
                else:
                    window = (window_start, window_end)
                    overlapping_segments = (segments[chrom]['start'][overlapping], segments[chrom]['end'][overlapping])
                    length = IBDSegmentsWeigher.overlap_segments_with_window(window, overlapping_segments)
                    overlaps[chrom].append((window_start, window_end, length))

        return dict(overlaps)

    @staticmethod
    def from_ibd_segments_file(path):
        """
        TODO: ...

        :param path: path to the IBD segments file to compute weight mask from
        """
        segments = IBDSegmentsWeigher.load_segments(path)
        overlaps = IBDSegmentsWeigher.compute_overlaps(segments)

        mask = {}
        for chrom in overlaps:
            lengths = np.array([length for _, _, length in overlaps[chrom]])
            X = lengths.reshape(-1, 1)
            try:
                cov = MinCovDet(random_state=0, support_fraction=IBDSegmentsWeigher.SUPPORT_FRACTION).fit(X)
                inliers = cov.reweight_covariance(X)[2]
                upper_boundary = max(lengths[inliers])
                weights = np.ones(len(lengths))
                sel = lengths > upper_boundary
                weights[sel] = weights[sel] * upper_boundary / lengths[sel]
                mask[chrom] = [(window[0], window[1], weight) for window, weight in zip(overlaps[chrom], weights)]
            except ValueError:
                mask[chrom] = [(window[0], window[1], 1) for window in overlaps[chrom]]

        return IBDSegmentsWeigher(mask)

    @staticmethod
    def from_json_mask_file(path):
        with open(path) as mask_file:
            mask = json.load(mask_file)

        return IBDSegmentsWeigher(mask)

    def __init__(self, mask):
        """
        TODO: ...
        """
        self.mask = {}
        for chrom in mask:
            self.mask[chrom] = {
                'start': np.fromiter((start for start, _, _ in mask[chrom]), float),
                'end': np.fromiter((end for _, end, _ in mask[chrom]), float),
                'weight': np.fromiter((weight for _, _, weight in mask[chrom]), float)
            }

        average_masked_length = []
        for chrom in self.mask:
            average_masked_length.append(
                sum((self.mask[chrom]['end'] - self.mask[chrom]['start']) * (1 - self.mask[chrom]['weight']))
            )

        self.adjusted_expected_recombination_number = \
            IBDSegmentsWeigher.ERSA_EXPECTED_RECOMBINATIONS_NUMBER - sum(average_masked_length) / 100

    def to_json_mask_file(self, path):
        """
        TODO: ...
        """

        mask = {}
        for chrom in self.mask:
            mask[chrom] = [window for window in
                zip(self.mask[chrom]['start'], self.mask[chrom]['end'], self.mask[chrom]['weight'])
            ]

        with open(path, 'w') as mask_file:
            json.dump(mask, mask_file)

    def visualize(self, path):
        plt.figure(figsize=(12, 6.5))
        ax = plt.gca()

        for chrom in self.mask:
            window_starts = self.mask[chrom]['start']
            window_ends = self.mask[chrom]['end']
            weights = self.mask[chrom]['weight']

            # Window shift of 0.1 is used to avoid white lines between horizontal bars
            xranges = list(zip(window_starts - 0.1, window_ends - window_starts + 0.2))
            yrange = (int(chrom) - 0.4, 0.8)

            normalize = colors.Normalize(vmin=0, vmax=1, clip=True)
            mapper = cm.ScalarMappable(norm=normalize, cmap=cm.jet)
            facecolors = list(map(lambda weight: mapper.to_rgba(weight), weights))
            barh_collection = BrokenBarHCollection(xranges, yrange, facecolors=facecolors)
            ax.add_collection(barh_collection)

        ax.autoscale()
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.invert_yaxis()
        ax.set_xlabel('cM', fontsize=12)

        # ERSA masked regions
        matplotlib.rcParams['hatch.linewidth'] = 0.5
        for chrom in IBDSegmentsWeigher.ERSA_MASKED_REGIONS:
            for start, end in IBDSegmentsWeigher.ERSA_MASKED_REGIONS[chrom]:
                rectangle = Rectangle(
                    (start, chrom - 0.5),
                    end - start,
                    1,
                    fill=False,
                    linewidth=0.5,
                    hatch='xxxxxxxx',
                    edgecolor='#ffb6c1'
                )
                ax.add_patch(rectangle)


        plt.colorbar(mapper, ax=ax, label='Weight')
        plt.yticks([int(x) for x in self.mask.keys()], map(lambda x: 'chr' + x, self.mask.keys()))
        plt.xticks(np.arange(0, 300, 50))
        plt.savefig(path)

    def compute_weighted_length(self, chrom, start, end, length):
        """
        TODO:
        """

        if chrom not in self.mask:
            return length

        overlapping = (self.mask[chrom]['end'] >= start) & (self.mask[chrom]['start'] <= end)

        if not np.any(overlapping):
            return length

        length_penalty = sum((
            np.minimum(end, self.mask[chrom]['end'][overlapping]) -
            np.maximum(start, self.mask[chrom]['start'][overlapping])
        ) * (1 - self.mask[chrom]['weight'][overlapping]))

        return length - length_penalty

    def apply_mask(self, input_segments_path, ibis_segments_path, output_segments_path):
        """
        TODO:
        """

        with (
            open(input_segments_path) as input_segments_file,
            open(ibis_segments_path) as ibis_segments_file,
            open(output_segments_path, 'w') as output_segments_file
        ):
            input_segments_reader = csv.reader(input_segments_file, delimiter='\t')
            ibis_segments_reader = csv.reader(ibis_segments_file, delimiter='\t')
            output_segments_writer = csv.writer(output_segments_file, delimiter='\t', quoting=csv.QUOTE_MINIMAL)

            for input_row, ibis_row in zip(input_segments_reader, ibis_segments_reader):
                chrom, length = input_row[2], float(input_row[6])
                start, end = float(ibis_row[6]), float(ibis_row[7])

                input_row[6] = self.compute_weighted_length(chrom, start, end, length)
                output_segments_writer.writerow(input_row)
