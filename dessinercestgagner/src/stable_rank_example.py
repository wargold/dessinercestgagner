import dessinercestgagner.src.stableRANK as sr
import numpy as np
from dessinercestgagner.src.Shape import Shape
import operator
import os
import midict
import matplotlib.pyplot as plt
import time
plt.style.use('ggplot')
inf = float("inf")
output_dir = os.path.join(os.path.dirname(__file__), '../../output/')
plt.style.use('ggplot')

method = 'move'
list_of_objects = ['apple', 'bat', 'bell']  # , 'bird'
colours_shape = {'apple': 'g', 'bat': 'k', 'bell': 'y'}
size_of_object = 21
perturb_magn = 0.3
perturb_number = 10

maxdim = 1
thresh = inf
coeff = 3
distance_matrix = False
do_cocycles = False
metric = "euclidean"
barcodes = {}

temp = None
contours = {"H0": [[[0], [1]], ("area", inf, inf)]}
border = 1
linewidth = 1.0
i = 0


def plot_signatures(signature_array, homology, title, colour):
    plt.figure(figsize=(20, 20))
    for k in signature_array:
        k[homology].plot(border, colour, linewidth)
    plt.title("Signature of " + title + ", " + homology + " homology")
    plt.show()
    return


def get_signatures(sample_shapes):
    stable_rank_max = {}
    signatures = {}
    for k in sample_shapes.keys():
        stable_rank_max[k] = sample_shapes[k]\
            .barcode(maxdim,
                     thresh,
                     coeff,
                     distance_matrix,
                     do_cocycles,
                     metric).stable_rank_max(contours)
    for shape_type in list_of_objects:
        signatures[shape_type] = {}
        sign_shape_H0 = []
        sign_shape_H1 = []
        sign_shape_H1_H0 = []
        for shape_number in range(1, size_of_object):
            for shape_perturb in range(perturb_number):
                temp = stable_rank_max[(shape_type, shape_number, shape_perturb)]
                sign_shape_H0.append(temp['H0'])
                sign_shape_H1.append(temp['H1'])
                sign_shape_H1_H0.append(temp['H1'] * temp['H0'].invert())
        signatures[shape_type]['H0'] = sign_shape_H0
        signatures[shape_type]['H1'] = sign_shape_H1
        signatures[shape_type]['H1/H0'] = sign_shape_H1_H0
    return  signatures

def Average_pcf(pcf_array, title, color, plot = False):
    a = 0
    for pcf in pcf_array:
        a += pcf
    a = a * (1 / len(pcf_array))
    if plot == True:
        plt.figure(figsize=(10, 10))
        border = 1
        linewidth = 1.0
        a.plot(border, color, linewidth)
        plt.title("Average PCF-" + title)
        plt.show()
    return a


def guess(nom, data_Type, types_mean, output):
    d = {}
    res = {fig_type: 0 for fig_type in types_mean.keys()}
    i = 0
    for pcf in nom:
        for type_name in types_mean.keys():
            d[type_name] = abs(pcf - types_mean[type_name]).integrate(0, 2)
        closest = min(d.items(), key=operator.itemgetter(1))[0]
        if output == True:
            plt.figure(figsize=(20,20))
            pcf.plot(1, color = 'b', linewidth = 1)
            for type_name in types_mean.keys():
                types_mean[type_name].plot(1, color = colours_shape[type_name], linewidth = 3)
            plt.title('real type:' + data_Type + ' guessed type:' + closest)
            plt.savefig(output_dir + 'guess_figure/' + data_Type + '_' + str(i) + '_' + closest + '.pdf')
            plt.show()
            plt.close()
        res[closest] = res[closest] + 1
    return res

def predict(signatures, homology, output = False):
    means = {}
    result = {}
    for object in list_of_objects:
        means[object] = Average_pcf(signatures[object][homology], object, colours_shape[object])
    for object in list_of_objects:
        print(object)
        result[object] = guess(signatures[object][homology], object, means, output)
    return result


def print_mean(normalized_sig_array, shape_name, colour):
    shape_mean = Average_pcf(normalized_sig_array[shape_name], shape_name, colour)
    plt.figure(figsize=(20, 20))
    for k in normalized_sig_array[shape_name]:
        k.plot(border, colour, linewidth)
    shape_mean.plot(border, 'y', linewidth + 10)
    plt.title("Signature of " + shape_name + ", normalized homology")
    plt.savefig(output_dir + 'signature_mean/' + shape_name + '_normalized_homology.pdf')
    plt.show()
    plt.close()
    return

if __name__ == '__main__':
    file = open(output_dir + 'result.txt', mode = 'w')
    start = time.time()
    for perturb_magn in np.arange(0,0.5,0.1):
        print(perturb_magn)
        print(time.time() - start)
        sample_shapes = {(k, y, i): sr.euc_object(Shape(k, y).perturb(method, perturb_magn).points) for k in list_of_objects
                                                                                             for y in range(1, size_of_object)                                                                    for i in range(perturb_number)
                     }
        print('Shapes created')
        print(time.time() - start)
        signatures = get_signatures(sample_shapes)
        print('Signatures computed')
        print(time.time()-start)
        output = False
        file.write(str(predict(signatures, homology = 'H1/H0', output = output)))
        del sample_shapes
        del signatures

    # Print mean and profile H1/H0
    # for object in list_of_objects:
    #     print_mean(object, 'b')