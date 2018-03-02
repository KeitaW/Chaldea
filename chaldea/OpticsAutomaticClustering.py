# -*- coding: utf-8 -*-
import numpy as np
from itertools import cycle
import scipy
import os
import sys
from optics import *
import AutomaticClustering as AutoC

class OpticsAutomaticClustering(object):
        def __init__(self, smoothing=9, epsilon=30, MinPts=50):
                self._smoothing = smoothing
                self._epsilon = epsilon
                self._MinPts = MinPts
        def fit(self, distmat):
                self._set_of_objects = setOfObjects(distmat, self._smoothing)
                prep_optics(self._set_of_objects, self._epsilon, self._MinPts)
                build_optics(self._set_of_objects, self._epsilon, self._MinPts)
                RPlot = []
                RPoints = []
                for obj in self._set_of_objects._ordered_list:
                    RPlot.append(self._set_of_objects._reachability[obj])
                    RPoints.append(obj)
                self.RPlot = RPlot
                # hierarchacally cluster the data
                rootNode = AutoC.automaticCluster(RPlot, RPoints)
                # array of the TreeNode objects, position in the array is the TreeNodes's level in the tree
                array = AutoC.getArray(rootNode, 0, [0])
                # get only leaves of the tree
                leaves = AutoC.getLeaves(rootNode, [])
                # extract clusters
                self.core_samples = []
                for item in leaves:
                    node = []
                    for v in range(item.start, item.end):
                        node.append(RPoints[v])
                    node = np.array(node)
                    self.core_samples.append(node)




