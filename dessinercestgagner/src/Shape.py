import numpy as np
import json
import copy
import matplotlib.pyplot as plt
from ripser import ripser
from persim import plot_diagrams
from dessinercestgagner.src import stableRANK as sr
import os
dir = os.path.dirname(__file__)


class Shape:

    def __init__(self, object_type, i : int = 1):
        '''
        Fetch the ith shape of object_type from json files
        :param object_type: String, Object in the data folder
        '''
        self.object_type = object_type
        dir = os.path.dirname(__file__)
        self.points = np.array([list(it.values()) for it in
                                json.load(open(os.path.join(dir, '../../dessinercestgagner/data/Shapes/' + object_type + '-' + str(i) + '.json')))['points']
                                ])

    def perturb(self, method = 'move', magn = 0.01):
        '''
        Method = 'move', modify points coordinates with U([-magn, magn])/2
        Method = 'noise', add some points with probability magn
        :param magn: float, magnitude of the modification
        :param method: str , in ('move', 'noise')
        :return: New Shape modified
        '''
        pert = copy.copy(self)
        if method == 'move':
            pert.points = pert.points + magn * ( np.random.rand(*self.points.shape) - 0.5 )
        if method == 'noise':
            pert.points = np.append(pert.points, np.random.rand(int(self.points.shape[0]*magn), 2), axis = 0)
        return pert

    def extend(self, n = 10, magn = 0.01):
        '''
        Perturbate the shape with random modification, return a list of n of those perturbations
        :param n: number of shapes returned
        :param magn: magnitude of the modification
        :return: list(Shape), original shape modified
        '''
        return [self.perturb(magn) for _ in range(n)]

    def render(self):
        '''
        Draw a matplotlib figure of the object
        :return: None
        '''
        plt.plot(self.points[:,0], self.points[:,1],'+')
        plt.show()
        return

    def plot_diagram(self):
        '''
        Plot birth-death diagram of the shape components (persistence diagram)
        :return: ax object
        '''
        diagrams = ripser(self.points)['dgms']
        return plot_diagrams(diagrams, show=True)

    def plot_barcode(self):
        '''
        Plot barcode of the shape
        :return:
        '''


if __name__ == "__main__":
    original_shape = Shape('apple').perturb(method = 'noise', magn = 0.05)
    shapes = original_shape.extend()
    for shape in shapes:
        shape.render()
        # shape.get_profile().display()
        shape.plot_diagram()


