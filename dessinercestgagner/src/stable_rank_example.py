import dessinercestgagner.src.stableRANK as sr
import numpy as np
from dessinercestgagner.src.Shape import Shape
inf=float("inf")
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import _pickle as pickle
from ripser import ripser
sample_shapes = {k: sr.euc_object(Shape(k).points) for k in ['apple', 'bat', 'beetle', 'bell']}
maxdim = 1
thresh = inf
coeff = 3
distance_matrix = False
do_cocycles = False
metric = "euclidean"
barcodes = {}
for k in sample_shapes.keys():
    barcodes[k] = sample_shapes[k].barcode(maxdim, coeff, metric)

plt.figure(figsize=(20, 20))
i = 1
for k in barcodes.keys():
    plt.subplot(2, 3, i)
    barcodes[k].plot("bar")
    i = i + 1
plt.show()