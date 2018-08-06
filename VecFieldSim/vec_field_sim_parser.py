# -*- coding: utf-8 -*-
"""
@author: mcreng
"""

import warnings
import numpy as np
import progressbar as pbar

warnings.filterwarnings("ignore")

class VecFieldSimParser:
    """
    This class parses *.his files into numpy arrays for class VecFieldSim.

    Sample usage:
        >>> v = VecFieldSimParser('empty.his', 'base.his', 'prop.his', 'full.his')
        >>> v.load_files() # Load in the data
        >>> v.calc_curv() # Calculate curvature
        >>> v.keep_local() # Remove data lying outside local block
        >>> v.remove_len() # Remove data with length of entries < 3
        >>> v.save() # Save to **_new.npy
    """
    def __init__(self, *filenames):
        """
        This initializes the class by supplying wind field *.his files.

        Args:
            filename1, filename2, ...: filenames
                Wind field data in *.his.

        Raises:
            TypeError: Incorrect type.
            ValueError: Incorrect file type.
        """
        if not (filenames and isinstance(filenames, tuple)):
            if not all(isinstance(elem, str) for elem in filenames):
                raise TypeError("files must be tuple of strings.")
        for filename in filenames:
            if not filename.endswith('.his'):
                raise ValueError(".his files needed.")

        self.filenames = filenames
        self.data = []
        self.files_cnt = len(filenames)

    @staticmethod
    def __find_header(filename):
        """
        This function finds the actual line number where data starts in *.his.

        Args:
            filename: filename
            The filename of the *.his file.

        Returns:
            The starting line number.

        Raises:
            EOFError: Incorrect document format.
        """
        with open(filename) as file:
            for num, line in enumerate(file, 2):
                if '---------------------------------------------' in line:
                    return num
        raise EOFError('Incorrect document format.')

    @staticmethod
    def __group(dat):
        """
        This function groups the numpy arrays by ID.

        Args:
            dat: numpy array

        Returns:
            Grouped numpy array
        """
        _, numeric_id = np.unique(dat[:, 0], return_inverse=True)
        _, cut_idx = np.unique(numeric_id, return_index=True)
        return np.array(np.split(dat, cut_idx)[1:])

    def load_files(self):
        """
        This function loads the *.his files into numpy arrays.

        Original Format:
            0              1     2     3     4     5     6     7    8
            ResidenceTime, XPos, YPos, ZPos, XVel, YVel, ZVel, ID, (COLORBY)

        New Format:
            0   1     2     3     4     5     6
            ID, XPos, YPos, ZPos, XVel, YVel, ZVel
        """
        widgets = ['LoadFiles: ', pbar.Percentage(), ' ',
                   pbar.Bar(marker='=', left='[', right=']'), ' ',
                   pbar.ETA()]
        pro_bar = pbar.ProgressBar(widgets=widgets, max_value=self.files_cnt)
        pro_bar.start()
        for i, filename in enumerate(self.filenames):
            pro_bar.update(i)
            self.data.append(self.__group(np.loadtxt(filename, delimiter='\t',
                                                     skiprows=self.__find_header(filename),
                                                     usecols=[7, 1, 2, 3, 4, 5, 6])))
        pro_bar.finish()

    def remove_len(self, length=3):
        """
        This function deletes entries with number smaller than length.

        Args:
            length: int

        Raises:
            TypeError: length must be int.
        """
        if not isinstance(length, int):
            raise TypeError('length must be int.')
        for i, dat in enumerate(self.data):
            index = np.empty(0, dtype=np.int32)
            for j, sub_dat in enumerate(dat):
                if len(sub_dat) >= length:
                    index = np.append(index, j)
            self.data[i] = dat[index]

    def keep_local(self):
        """
        This function removes entries that lie outside the local block.
        """
        for i, dat in enumerate(self.data):
            for j, sub_dat in enumerate(dat):
                sub_dat = sub_dat[7*sub_dat[:, 2]-10*sub_dat[:, 1] <= 36280 + 98.5*7]
                sub_dat = sub_dat[7*sub_dat[:, 2]-10*sub_dat[:, 1] >= 34817 - 98.5*7]
                sub_dat = sub_dat[10*sub_dat[:, 2]+7*sub_dat[:, 1] <= 146550 + 98.5*10]
                sub_dat = sub_dat[10*sub_dat[:, 2]+7*sub_dat[:, 1] >= 144464 - 98.5*10]
                dat[j] = sub_dat
            self.data[i] = dat
        self.remove_len(length=1)

    def calc_curv(self):
        """
        This function calculates the curvatures in the numpy array.

        Original Format:
            0   1     2     3     4     5     6
            ID, XPos, YPos, ZPos, XVel, YVel, ZVel

        New Format:
            0   1     2     3     4     5     6     7
            ID, XPos, YPos, ZPos, XVel, YVel, ZVel, Curv
        """
        widgets = ['CalcCurv: ', pbar.Percentage(), ' ',
                   pbar.Bar(marker='=', left='[', right=']'), ' ',
                   pbar.ETA()]
        pro_bar = pbar.ProgressBar(widgets=widgets, max_value=self.files_cnt)
        pro_bar.start()
        for i, dat in enumerate(self.data):
            pro_bar.update(i)
            for j, sub_dat in enumerate(dat):
                dt1 = sub_dat[:, [4, 5, 6]]
                if len(dt1) < 3:
                    continue
                dist = sub_dat[:, [1, 2, 3]]
                dist = np.diff(dist, axis=0)
                dist = np.hypot(dist[:, 0], dist[:, 1])
                dist = np.cumsum(dist)
                dist = np.insert(dist, 0, 0)
                dt2 = np.gradient(dt1, dist, axis=0, edge_order=2)
                norm = np.linalg.norm(np.cross(dt1, dt2), axis=1)
                curv = norm / (np.linalg.norm(dt1, axis=1)**3)
                curv[curv < 0] = 0
                temp = np.zeros((sub_dat.shape[0], sub_dat.shape[1]+1))
                temp[:, :-1] = sub_dat
                temp[:, -1] = curv
                temp = temp[np.logical_and(np.logical_not(np.isnan(temp[:, -1])),
                                           np.logical_not(np.isinf(temp[:, -1])))]
                dat[j] = temp
            self.data[i] = dat
        pro_bar.finish()

    def save(self, suffix='_new'):
        """
        This function saves the numpy arrays, appending suffix in the filename.

        Args:
            suffix: str
                Suffix to be appended.

        Raises:
            TypeError: Suffix needs to be str.
        """
        if not isinstance(suffix, str):
            raise TypeError('Suffix needs to be str.')
        for i, filename in enumerate(self.filenames):
            np.save(filename[:-4]+suffix+'.npy', self.data[i])
        