import numpy as np
import json
import copy
import matplotlib.pyplot as plt
from ripser import ripser
from persim import plot_diagrams
from dessinercestgagner.src import stableRANK as sr
import os
dir = os.path.dirname(__file__)
filename = os.path.join(dir, '/relative/path/to/file/you/want')


class Shape:

    def __init__(self, object_type, i : int = 1):
        '''
        Fetch the ith shape of object_type from json files
        :param object_type: String, Object in the data folder
        '''
        self.object_type = object_type
        dir = os.path.dirname(__file__)
        print(os.path.join(dir, '../../dessinercestgagner/data/Shapes/' + object_type + '-' + str(i) + '.json'))
        self.points = np.array([list(it.values()) for it in
                                json.load(open(os.path.join(dir, '../../dessinercestgagner/data/Shapes/' + object_type + '-' + str(i) + '.json')))['points']
                                ])

    def perturb(self, magn = 0.01):
        '''
        Modify points coordinates with U([-magn, magn])/2
        :param magn: float, magnitude of the modification
        :return: New Shape modified
        '''
        pert = copy.copy(self)
        pert.points = pert.points + magn * ( np.random.rand(*self.points.shape) - 0.5 )
        return pert

    def render(self):
        '''
        Draw a matplotlib figure of the object
        :return: None
        '''
        plt.plot(self.points[:,0], self.points[:,1],'+')
        plt.show()
        return

    def plot_barcode(self):
        '''
        Plot barcode of the shape
        :return: ax object
        '''
        diagrams = ripser(self.points)['dgms']
        return plot_diagrams(diagrams, show=True)


