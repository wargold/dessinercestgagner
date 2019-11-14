# Dessinercestgagne

Dessinercestgagner aim at classifying 2D-shapes using TDA methods

## Installation

Simply clone the repo, make sure ripser is installed.

## Usage

Example is provided in classifier.py

## Output
### Objects
2D-representation of the objects to be classified, names are constructed as follow:
* {shape type e.g. 'apple'}\_{Transformation type e.g. 'move'}\_{modification magnitude e.g. 0.1}.pdf
### Guess_figure
H1/H0 signature of the object superposed with the mean signature of the three types.
File names are constructed as follow:
* {shape type e.g. 'apple'}\_{Shape number e.g. 1,2,3}\_{Classification result e.g. 'bell'}.pdf
### signature_mean
All the normalized H1/H0 signatures and the mean of these

## License
[MIT](https://choosealicense.com/licenses/mit/)
