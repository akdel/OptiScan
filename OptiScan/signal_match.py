from scipy import signal, ndimage
from numpy import array, median, ones

"""
OptiScan crosscorr based aligner. This is different from the Photoasm nick label aligner.
"""


class Matcher:
    """
    Matches two signals with best possible positions. Returns the position info and the degree of match
    """
    def __init__(self, signalA, signalB):
        self.A = array(signalA, dtype=float)
        self.B = array(signalB, dtype=float)
        if len(self.A) >= len(self.B):
            self.primary = self.A
            self.primary_index = 1
            self.secondary = self.B
            self.secondary_index = 2
        else:
            self.primary = self.B
            self.primary_index = 2
            self.secondary = self.A
            self.secondary_index = 1
        self.closest_match = int()
        self.match_index = int()
        self.secondary_match_index = int()
        self.match = None
        self.score = None
        self.adjusted_match_score = None
        self.contrast = None

    def convolve_signals_and_get_match_info(self):
        convolution_product = array(signal.correlate(self.primary, self.secondary, mode="same"), dtype=int)
        self.score = ndimage.white_tophat(convolution_product, structure=ones(7))
        self.match_index = convolution_product.argmax()
        self.closest_match = self.score[self.match_index]
        secondary_half_length = len(self.secondary)/2
        if self.match_index < secondary_half_length:
            secondary_match = self.secondary[secondary_half_length - self.match_index:]
            start_prim = 0
            self.secondary_match_index = (len(self.secondary) - secondary_half_length - self.match_index)/2
        else:
            secondary_match = self.secondary
            start_prim = self.match_index - secondary_half_length
        primary_match = self.primary[start_prim:self.match_index + secondary_half_length]
        self.match = (primary_match, secondary_match)

    def get_contrast(self):
        tophat_score = ndimage.white_tophat(self.score, structure=ones(10))
        self.score = tophat_score
        self.closest_match = max(self.score)
        med = median(tophat_score)
        self.contrast = self.closest_match/med
        self.adjusted_match_score = self.contrast/len(self.match[0])
