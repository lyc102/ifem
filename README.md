# iFEM: an integrated finite element method package in MATLAB
iFEM is a MATLAB software package containing robust, efficient, and easy-following codes for the main building blocks of adaptive finite element methods on unstructured simplicial grids in both two and three dimensions. Besides the simplicity and readability, sparse matrixlization, an innovative programming style for MATLAB, is introduced to improve the efficiency. In this novel coding style, the sparse matrix and its operations are used extensively in the data structure and algorithms.


## Installation

Add the path to iFEM into the path library of MATLAB by either:

1. Graphical interface. Click
	File -> Set Path -> Add with Subfolders
   and chose the directory where the package iFEM is stored.
  
2. Command window. Go to the directory of iFEM and run 
`setpath`


## Help

1. `help funexample` displays a description of and syntax for the function `funexample`. For example,
`help mg` will show basic usage for `mg` function in the plain text.  
2. `ifem funexampledoc` show detailed description. For example, `ifem mgdoc` will explain the `mg` function step by step in html format. But not every function has a html documentation.


## Quick Start

1. Type `ifem introduction` to get an introduction on ifem.

2. Go through examples in `\example` directory.


## Feedback

If you like it, please send me an email lyc102@gmail.com. If you  feel it is helpful for your research, please acknowledge your use by citing:

> L. Chen. iFEM: an integrated finite element method package in MATLAB. Technical Report, University of California at Irvine, 2009.

```bibtex
@techreport{Chen:2008ifem,
	author = {Long Chen},
	journal = {Technical Report, University of California at Irvine},
	title = {{$i$FEM}: an integrated finite element methods package in {MATLAB}},
	url = {https://github.com/lyc102/ifem},
	year = {2009}}
```


## Acknowledgement

The author would like to thank Professor Michael Holst in University of California at San Diego and Professor Ludmil Zikatanov in Pennsylvania State University for many insightful discussion, and also Professor Chensong Zhang in Chinese Academy of Sciences for the effort in the development of AFEM@matlab, an early version of iFEM.

The author thanks students or postdocs Shuhao Cao, Ming Wang, Huayi Wei, Lin Zhong, and Jie Zhou for their contribution to iFEM in one way or another. Detailed credits can be found in the M-lint of several m files.

The author is also grateful to the NSF for the partial support over the years. 


Long Chen

--

Professor                 

Department of Mathematics
University of California at Irvine
http://math.uci.edu/~chenlong/

--