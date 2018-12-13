# _*_ coding：utf-8 _*_
'''
Codes for signal processing
e.g. get train data, preprocessing...
'''

import logging
from sklearn import preprocessing
import numpy as np
from scipy.signal import stft
from scipy import interpolate
import matplotlib.pyplot as plt

from LoadData import EEGDataSet


class SigProcess:
    # @para dataset dataset with data and label
    # @para start start time in seconds
    # @para end end time in seconds
    def __init__(self, dataset_class, start, end):
        self.__sfreq = dataset_class.sampleRate
        dataset = dataset_class.dataset
        if start < 0 or start >= 5 or end < 0 or end > 5 or start >= end:
            print("session index not valid")
            raise ValueError
        self.__sessions = "session " + str(start) + " to " + "session " + str(end-1)
        self.__dataInfo = list()
        self.__dataInfo.append(self.__sessions)
        self.data = dict()
        for subject in dataset:
            self.data[subject] = dict()
            self.data[subject]['origin'] = list()
            self.data[subject]['label'] = list()
            for session in range(start, end):
                self.data[subject]['origin'].extend(dataset[subject][session]['data'])
                self.data[subject]['label'].extend(dataset[subject][session]['label'])
        self.__dataInfo.extend(['origin', 'label'])

    # read get training data from original data
    # @para start start time in seconds
    # @para end end time in seconds
    # @para slide window slide time in seconds
    # @para window segment length in seconds
    def get_train_data(self, start=4.5, end=6.5, slide=0.1, window=2):
        for subject in self.data:
            train_data = list()
            train_label = list()
            for idx in range(len(self.data[subject]['origin'])):
                segments = self.__get_segments(self.data[subject]['origin'][idx], self.data[subject]['label'][idx],
                                            start, end, slide, window)

                train_data.append(segments[0])
                train_label.append(segments[1])
            self.data[subject]['train_data'] = train_data
            self.data[subject]['train_label'] = train_label
        self.__dataInfo.extend(['train_data', 'train_label'])

    # get segments data from one trial
    # @return tuple with segments and labels
    def __get_segments(self, ori_data, label, start, end, slide, window):
        if slide == 0:
            logging.warning("parameter slide cannot set 0 in SigProcess.getSegments()")
            slide = 0.1
        segments = list()
        segments_label = list()
        start = int(start * self.__sfreq)
        end = int(end * self.__sfreq)
        slide = int(slide * self.__sfreq)
        window = int(window * self.__sfreq)
        left = start
        right = start + window
        while right <= end:
            # averaged
            segments.append(ori_data[left:right, :])
            left += slide
            right += slide
            # segments_label.append(label)
        return segments, label

    # train data standardization
    def data_standardization(self):
        for subject in self.data:
            self.data[subject]['standard_data'] = list()
            for segments in self.data[subject]['train_data']:
                standard_segments = list()
                for segment_data in segments:
                    #  Scaling data with zero mean and unit variance
                    standard_data = preprocessing.scale(segment_data, axis=0)
                    scale_data = preprocessing.minmax_scale(standard_data, (0, 1), axis=0)
                    standard_segments.append(scale_data)
                self.data[subject]['standard_data'].append(standard_segments)
        self.__dataInfo.append('standard_data')

    def get_data_info(self):
        return self.__dataInfo

    # get specified data/labels of specified subject
    def get_specified_data(self, subject, data_name):
        if subject not in self.data or data_name not in self.data[subject]:
            print("data information:" + str(self.get_data_info()))
            logging.warning("cannot find specified subject or data name")
        return self.data[subject][data_name]

    # find kth largest number in a 1-d array
    @staticmethod
    def find_kth_largest(arr, k):
        k = k - 1
        lo = 0
        hi = len(arr) - 1
        while lo < hi:
            arr[lo], arr[int((lo + hi) / 2)] = arr[int((lo + hi) / 2)], arr[lo]
            left = lo
            right = hi
            pivot = arr[lo]
            while left < right:
                while left < right and arr[right] <= pivot:
                    right = right - 1
                arr[left] = arr[right]
                while left < right and arr[left] >= pivot:
                    left = left + 1
                arr[right] = arr[left]
            arr[left] = pivot
            if k <= left:
                hi = left - 1
            if k >= left:
                lo = left + 1
        res = arr[k]
        return res

    # rescale data
    # @percentage: percentage of value considered to be artifact
    @staticmethod
    def recale(data, percentage):
        m, n = data.shape
        arr = data.flatten()
        min_val = np.min(arr)
        max_val = SigProcess.find_kth_largest(arr, int(m * n * percentage))
        for i in range(m):
            for j in range(n):
                if data[i][j] > max_val:
                    data[i][j] = 1
                else:
                    data[i][j] = (data[i][j] - min_val) / (max_val - min_val)
        return data

class FeatureExtraction:
    def __init__(self, data, label):
        self.data = dict()
        self.data['train_data'] = data
        self.data['train_label'] = label

    def stft_process(self, mu=(4, 14), beta=(16, 32), wlen=64, nfft=512, fs=250, hop=14, rescale=False):
        if len(mu) != 2 or len(beta) != 2:
            logging.warning("mu and beta must be tuple or list with 2 integer elements")
            raise ValueError
        self.data['stft'] = list()
        self.data['stft_channel'] = list()
        freq_falg = True  # reduce compute complexity
        for segments in self.data['train_data']:
            segment_image = list()
            segment_channel_image = list()
            for segment_data in segments:
                stft_image = None
                channel_image = list()
                for chn in range(len(segment_data[0])):
                    f, t, Fstft = stft(segment_data[:, chn], fs=fs, window='hamming', nperseg=wlen, noverlap=wlen - hop,
                                       nfft=nfft, return_onesided=True, boundary=None, padded=False)
                    if freq_falg:
                        freq_falg = False
                        mu_left = np.where(f >= mu[0])[0][0]
                        mu_right = np.where(f >= mu[1])[0][0]
                        beta_left = np.where(f >= beta[0])[0][0]
                        beta_right = np.where(f >= beta[1])[0][0]
                    mu_feature_matrix = np.abs(Fstft[mu_left : mu_right])
                    beta_feature_matrix = np.abs(Fstft[beta_left : beta_right])

                    # beta band cubic interpolation
                    beta_interp = interpolate.interp2d(t, f[beta_left: beta_right], beta_feature_matrix, kind='cubic')
                    interNum = len(mu_feature_matrix)
                    f_beta = np.arange(beta[0], beta[1], (beta[1] - beta[0]) / (interNum))
                    beta_feature_matrix = beta_interp(t, f_beta)

                    # Note: recale feature matrix
                    if rescale == True:
                        mu_feature_matrix = SigProcess.recale(mu_feature_matrix, 0.05)
                        beta_feature_matrix = SigProcess.recale(beta_feature_matrix, 0.05)

                    # plt.pcolormesh(t, f_beta, beta_feature_matrix, vmin=0)
                    # plt.show()
                    # pause = input("pause")
                    channel_image.append(np.append(mu_feature_matrix, beta_feature_matrix, axis=0))
                    if stft_image is None:
                        stft_image = np.append(mu_feature_matrix, beta_feature_matrix, axis=0)
                    else:
                        stft_image = np.append(stft_image, mu_feature_matrix, axis=0)
                        stft_image = np.append(stft_image, beta_feature_matrix, axis=0)
                segment_image.append(stft_image)
                segment_channel_image.append(channel_image)
            self.data['stft'].append(segment_image)
            self.data['stft_channel'].append(segment_channel_image)
        self.data['stft'] = np.array(self.data['stft'])
        self.data['stft_channel'] = np.array(self.data['stft_channel'])

    def stft_standard(self):
        if 'stft' not in self.data:
            logging.warning("must run FeatureExtraction.stft_process() first...")
            return None
        self.data['std_stft'] = list()
        for segment_images in self.data['stft']:
            scaled_segment_images = list()
            for image in segment_images:
                #  Scaling stft data with zero mean and unit variance
                standard_image = preprocessing.scale(image, axis=1)
                scaled_image = preprocessing.minmax_scale(standard_image, (0, 1), axis=1)
                # plt.pcolormesh(scaled_image, vmin=0)
                # plt.show()
                # pause = input("pause")
                scaled_segment_images.append(scaled_image)
            self.data['std_stft'].append(scaled_segment_images)



if __name__ == '__main__':
    dataSetDir = r"D:\BCI\EEG dataset\bci competition IV 2b\Graz_mat"
    dataset = EEGDataSet(dataSetDir)
    process = SigProcess(dataset, 0, 5)
    process.get_train_data()
    # print(process.get_data_info())
    # print(len(process.data['sub01']['train_data']))
    # label = process.data['sub01']['train_label']
    # class_num = [i for i in label if i == 1]
    process.data_standardization()
    print(process.get_data_info())
    fe = FeatureExtraction(process.data['sub01']['train_data'], process.data['sub01']['train_label'])
    fe.stft_process()
    print(fe.data['stft_channel'].shape)
    fe.stft_standard()
    print('Finish signal processing...')
