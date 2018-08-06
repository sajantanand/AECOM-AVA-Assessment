# -*- coding: utf-8 -*-
"""
@author: mcreng
"""

import copy
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import progressbar as pbar
import cv2

empty = np.load('500/empty_new.npy')
full = np.load('500/full_new.npy')
base = np.load('500/base_new.npy')
proposed = np.load('500/prop_new.npy')

class VecFieldSim:
    """
    This class reads wind field data and finds their similarity.

    Sample usage:
        >>> v = VecFieldSim(empty, base, proposed,
                            labels=['Empty Scheme', 'Base Scheme', 'Proposed Scheme'])
        >>> v.calc_dist() # returns distances of histograms
        >>> v.plot_hist() # plots histograms
        >>> v.test_conv(n=50) # run comparisons for 50 times to test convergence
    """

    def __init__(self, *data, **kwargs):
        """
        This initializes the class by supplying wind field data.

        Args:
            data1, data2, ..., datan: numpy arrays
                Wind field data, assumes the first input is the baseline.
            labels: list of strings
                Use to label the data, must be in same order of data1, data2, ..., datan

        Raises:
            TypeError: Incorrect type.
            ValueError: Count in data and labels mismatch.
        """
        for dat in data:
            if not isinstance(dat, np.ndarray):
                raise TypeError("*data must be numpy arrays.")
        labels = kwargs.pop('labels')
        if not (labels and isinstance(labels, tuple)):
            if not all(isinstance(elem, str) for elem in labels):
                raise TypeError("labels must be tuple of strings.")
        if len(data) != len(labels):
            raise ValueError("Count in data and labels mismatch.")
        self.data = data
        self.data_cnt = len(data)
        self.labels = labels
        self.arr_curv_dist = []
        self.arr_vel_dist = []
        self.hist_curv_dist = []
        self.hist_vel_dist = []

        # Setup plt
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')

    @staticmethod
    def __read_curv(dat):
        """
        This reads the data of a scheme and returns array of mean of curvature against distance.

        Args:
            dat: numpy array
                Wind field data.

        Returns:
            Array of mean of curvature against distance.
        """
        p_dist = np.empty((0, 2))
        for _, sub_dat in enumerate(dat):
            length = len(sub_dat)
            ran_dat_1 = sub_dat[np.random.choice(sub_dat.shape[0], int(length*.8))]
            ran_dat_2 = sub_dat[np.random.choice(sub_dat.shape[0], int(length*.8))]
            curv = np.mean([ran_dat_1[:, -1], ran_dat_2[:, -1]], axis=0)
            dist = np.sqrt(np.sum((ran_dat_1[:, 1:4]-ran_dat_2[:, 1:4])**2, axis=1))
            p_dist = np.append(p_dist, np.vstack((curv, dist)).T, axis=0)
        return p_dist

    @staticmethod
    def __read_vel(dat):
        """
        This reads the data of a scheme and returns array of mean of curvature against distance.

        Args:
            dat: numpy array
                Wind field data.

        Returns:
            Array of mean of curvature against distance.
        """
        p_dist = np.empty((0, 2))
        for _, sub_dat in enumerate(dat):
            length = len(sub_dat)
            ran_dat_1 = sub_dat[np.random.choice(sub_dat.shape[0], int(length))]
            ran_dat_2 = sub_dat[np.random.choice(sub_dat.shape[0], int(length))]
            vdot = np.einsum('ij,ij->i', ran_dat_1[:, 4:7], ran_dat_2[:, 4:7])
            dist = np.sqrt(np.sum((ran_dat_1[:, 1:4]-ran_dat_2[:, 1:4])**2, axis=1))
            p_dist = np.append(p_dist, np.vstack((vdot, dist)).T, axis=0)
        return p_dist

    def __update_arrays(self, mode=2, show_bar=True):
        """
        This updates the quantity against distance arrays depending on mode.

        Args:
            mode: 0, 1, 2
                0 - Updates mean of cuvature against distance only;
                1 - Updates velocity dot product against distance only;
                2 - Updates both arrays.
        """
        if mode == 0:
            self.arr_curv_dist = []
            widgets = ['UpdateArray-Curv: ', pbar.Percentage(), ' ',
                       pbar.Bar(marker='=', left='[', right=']'), ' ',
                       pbar.ETA()]
            if show_bar:
                pro_bar = pbar.ProgressBar(widgets=widgets, max_value=self.data_cnt)
                pro_bar.start()
            for i, dat in enumerate(self.data):
                if show_bar:
                    pro_bar.update(i)
                self.arr_curv_dist.append(self.__read_curv(dat))
            if show_bar:
                pro_bar.finish()
        elif mode == 1:
            self.arr_vel_dist = []
            widgets = ['UpdateArray-Vel : ', pbar.Percentage(), ' ',
                       pbar.Bar(marker='=', left='[', right=']'), ' ',
                       pbar.ETA()]
            if show_bar:
                pro_bar = pbar.ProgressBar(widgets=widgets, max_value=self.data_cnt)
                pro_bar.start()
            for i, dat in enumerate(self.data):
                if show_bar:
                    pro_bar.update(i)
                self.arr_vel_dist.append(self.__read_vel(dat))
            if show_bar:
                pro_bar.finish()
        else:
            self.__update_arrays(mode=0)
            self.__update_arrays(mode=1)

    def __calc_hist(self, mode=2, show_bar=True):
        """
        This calculates the 2D histograms of the quantity against distance arrays depending on mode.

        Args:
            mode: 0, 1, 2
                0 - Calculates the one for mean of cuvature against distance only;
                1 - Calculates the one for velocity dot product against distance only;
                2 - Calculates for both.
        """
        self.__update_arrays(mode, show_bar)
        if mode == 0:
            self.hist_curv_dist = []
            arr_p_min = []
            arr_p_max = []
            arr_d_max = []
            for _, dat in enumerate(self.arr_curv_dist):
                arr_p_min.append(np.percentile(dat[:, 0], 0))
                arr_p_max.append(np.percentile(dat[:, 0], 65))
                arr_d_max.append(np.percentile(dat[:, -1], 90))
            bin_n = (np.linspace(0, max(arr_d_max), 50),
                     np.linspace(min(arr_p_min), max(arr_p_max), 50))
            for _, dat in enumerate(self.arr_curv_dist):
                hist = list(np.histogram2d(dat[:, -1], dat[:, 0], bins=bin_n))
                hist[0] = (hist[0] / np.sum(hist[0])).astype(np.float32)
                self.hist_curv_dist.append(hist)
        elif mode == 1:
            self.hist_vel_dist = []
            arr_p_min = []
            arr_p_max = []
            arr_d_max = []
            for _, dat in enumerate(self.arr_vel_dist):
                arr_p_min.append(np.percentile(dat[:, 0], 10))
                arr_p_max.append(np.percentile(dat[:, 0], 90))
                arr_d_max.append(np.percentile(dat[:, -1], 90))
            bin_n = (np.linspace(0, max(arr_d_max), 50),
                     np.linspace(min(arr_p_min), max(arr_p_max), 50))
            for _, dat in enumerate(self.arr_vel_dist):
                hist = list(np.histogram2d(dat[:, -1], dat[:, 0], bins=bin_n))
                hist[0] = (hist[0] / np.sum(hist[0])).astype(np.float32)
                self.hist_vel_dist.append(hist)
        else:
            self.__calc_hist(mode=0)
            self.__calc_hist(mode=1)

    def calc_dist(self, mode=2, calc_hist=False, show_bar=True):
        """
        This calculates the distances of histograms depending on mode.

        Args:
            mode: 0, 1, 2
                0 - Calculates the one for mean of cuvature against distance only;
                1 - Calculates the one for velocity dot product against distance only;
                2 - Calculates for both.
            calc_hist: boolean
                True - Recalcuates the histograms before plotting;
                False - Plot histograms only. However, if no previous calculation
                has been done, this will be overwritten to True.

        Returns:
            Mode 0 - The n by 3 array of metrics from cuvature
            Mode 1 - The n by 3 array of metrics from velocity
            Mode 2 - A list of two n by 3 arrays of metrics from [cuvature, velocity]

            All arrays are in format of:
                    CHISQR   CHISQR_ALT   HELLINGER
            D(0,0)     *         *            *
            D(0,1)     *         *            *
            D(0,2)     *         *            *
             ...      ...       ...          ...
            D(0,n)     *         *            *

        Raises:
            ValueError: Incorrect mode.
        """
        if mode != 0 and mode != 1 and mode != 2:
            raise ValueError('Incorrect mode.')
        if (mode == 0 or mode == 2) and not self.hist_curv_dist:
            calc_hist = True
        if (mode == 1 or mode == 2) and not self.hist_vel_dist:
            calc_hist = True
        if calc_hist:
            self.__calc_hist(mode, show_bar)
        if mode == 0:
            metric = np.empty((0, 3))
            for _, hist in enumerate(self.hist_curv_dist):
                metric = np.append(metric, [
                    [cv2.compareHist(self.hist_curv_dist[0][0], hist[0], cv2.HISTCMP_CHISQR),
                     cv2.compareHist(self.hist_curv_dist[0][0], hist[0], cv2.HISTCMP_CHISQR_ALT)/4,
                     cv2.compareHist(self.hist_curv_dist[0][0], hist[0], cv2.HISTCMP_BHATTACHARYYA)]
                    ], axis=0)
            return metric
        elif mode == 1:
            metric = np.empty((0, 3))
            for _, hist in enumerate(self.hist_vel_dist):
                metric = np.append(metric, [
                    [cv2.compareHist(self.hist_vel_dist[0][0], hist[0], cv2.HISTCMP_CHISQR),
                     cv2.compareHist(self.hist_vel_dist[0][0], hist[0], cv2.HISTCMP_CHISQR_ALT)/4,
                     cv2.compareHist(self.hist_vel_dist[0][0], hist[0], cv2.HISTCMP_BHATTACHARYYA)]
                    ], axis=0)
            return metric
        return [self.calc_dist(mode=0), self.calc_dist(mode=1)]

    @staticmethod
    def __factor_int(number):
        """
        This factorize an integer into two factors, with the difference between minimized.

        Args:
            number: int
                An integer to be factorized.

        Returns:
            val1, val2: The required factors.
        """
        solution = False
        val1 = math.ceil(math.sqrt(number))
        while not solution:
            val2 = int(number/val1)
            if val2 * val1 == float(number):
                solution = True
            else:
                val1 -= 1
        return val1, val2

    def plot_hist(self, mode=2, calc_hist=False, save=False):
        """
        This plots the histograms depending on mode

        Args:
            mode: 0, 1, 2
                0 - Plots the one for mean of cuvature against distance only;
                1 - Plots the one for velocity dot product against distance only;
                2 - Plots both.
            calc_hist: boolean
                True - Recalcuates the histograms before plotting;
                False - Plot histograms only. However, if no previous calculation
                has been done, this will be overwritten to True.

        Raises:
            ValueError: Incorrect mode.
        """
        if mode != 0 and mode != 1 and mode != 2:
            raise ValueError('Incorrect mode.')
        if (mode == 0 or mode == 2) and not self.hist_curv_dist:
            calc_hist = True
        if (mode == 1 or mode == 2) and not self.hist_vel_dist:
            calc_hist = True
        if calc_hist:
            self.__calc_hist(mode)

        if mode == 2:
            self.plot_hist(mode=0, save=save)
            self.plot_hist(mode=1, save=save)
            return

        if mode == 0:
            param = 'Mean of Curvature'
        else:
            param = 'Velocity Dot Product'

        fig = plt.figure(param, figsize=(8, 4), dpi=700)
        size = self.__factor_int(self.data_cnt)
        grid = plt.GridSpec(size[0], size[1]*2, hspace=0.4, wspace=0.4)
        ax = []
        cur_x = 0
        cur_y = 0
        for i in range(self.data_cnt):
            if ax != []:
                ax.append(fig.add_subplot(grid[cur_x, cur_y:cur_y+2], sharex=ax[0], sharey=ax[0]))
            else:
                ax.append(fig.add_subplot(grid[cur_x, cur_y:cur_y+2]))
            ax[i].set_title(self.labels[i])
            cur_y += 2
            if cur_y == size[1]*2:
                cur_y = 0
                cur_x += 1

        t = np.arange(self.data_cnt).reshape(size)
        for i in t[:, 0]:
            ax[i].set_ylabel(param)
        for i in t[-1]:
            ax[i].set_xlabel('Euclidean Distance')

        my_cmap = copy.copy(plt.cm.get_cmap('afmhot'))
        my_cmap.set_bad((0, 0, 0))

        if mode == 0:
            sup = np.max([np.max(hist[0]) for hist in self.hist_curv_dist])
            for i, hist in enumerate(self.hist_curv_dist):
                cb = ax[i].imshow(hist[0].T[::-1],
                                  extent=[hist[1][0], hist[1][-1], hist[2][0], hist[2][-1]],
                                  aspect='auto', cmap=my_cmap, vmax=sup, norm=LogNorm())
            fig.colorbar(cb, ax=ax)
        elif mode == 1:
            sup = np.max([np.max(hist[0]) for hist in self.hist_vel_dist])
            for i, hist in enumerate(self.hist_vel_dist):
                cb = ax[i].imshow(hist[0].T[::-1],
                                  extent=[hist[1][0], hist[1][-1], hist[2][0], hist[2][-1]],
                                  aspect='auto', cmap=my_cmap, vmax=sup, norm=LogNorm())
            fig.colorbar(cb, ax=ax)

        plt.show()
        if save:
            plt.savefig(param+".png", bbox_inches='tight')

    def test_conv(self, n=50, mode=2, save=False):
        """
        This tests the convergence of metrics by running the comparisons multiple times.

        Args:
            n: int > 0
                The number of comparisons to make in the test.
            mode: 0, 1, 2
                0 - Plots the one for mean of cuvature against distance only;
                1 - Plots the one for velocity dot product against distance only;
                2 - Plots both.

        Returns:
            An array of means of metrics in format of
            
                        D(0,0) D(0,1) D(0,2) ... D(0,n)
            CHISQR        *      *      *          *
            CHISQR_ALT    *      *      *          *
            HELLINGER     *      *      *          *

        Raises:
            ValueError: Incorrect inputs.
        """
        if mode != 0 and mode != 1 and mode != 2:
            raise ValueError('Incorrect mode.')
        if mode == 2:
            return self.test_conv(n=n, mode=0, save=save), self.test_conv(n=n, mode=1, save=save)

        if mode == 0:
            param = 'Mean of Curvature'
        else:
            param = 'Velocity Dot Product'
        metric_names = [r'$D_1$ ($\chi^2$)', r'$D_2$ (Alternative $\chi^2$)', r'$D_3$ (Hellinger)']
        fig = plt.figure(param, figsize=(2, 3), dpi=700)
        grid = plt.GridSpec(3, 1, hspace=1)
        axs = []
        for i in range(3):
            axs.append(fig.add_subplot(grid[i, 0]))
            axs[i].get_yaxis().set_visible(False)
            axs[i].set_ylim(-0.3, (self.data_cnt-1)*0.5+0.3)

        axs[1].set_xlim(-0.03, 1.03)
        axs[2].set_xlim(-0.03, 1.03)

        colors = ['gx', 'rx', 'bx', 'yx', 'mx', 'cx', 'kx']
        
        means = np.zeros((3, self.data_cnt))

        widgets = [pbar.Percentage(), ' ',
                   pbar.Bar(marker='=', left='[', right=']'), ' ',
                   pbar.ETA()]
        pro_bar = pbar.ProgressBar(widgets=widgets, max_value=3*n*self.data_cnt)
        pro_bar.start()
        for k in range(n):
            pro_bar.update(k)
            values = self.calc_dist(mode, calc_hist=True, show_bar=False)
            means += values.T
            for i in range(3):
                for j in range(self.data_cnt):
                    axs[i].plot(values[j, i], (self.data_cnt-j-1)*0.5, colors[j%len(colors)])
                    axs[i].set_title(metric_names[i])

        means /= n

        pro_bar.finish()
        plt.show()
        if save:
            plt.savefig(param+".png", bbox_inches='tight')
            np.save(param+'.npy', means)
        return means
