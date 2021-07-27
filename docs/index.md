---
layout: single
permalink: /
title: "iFEM: an Integrated Finite Element Methods Package in MATLAB"
excerpt: "A quick start guide."
sidebar:
    nav: docs
toc: true
toc_sticky: true
---

## Introduction to iFEM
*i*FEM is a [MATLAB](https://www.mathworks.com/) software package containing robust, efficient, and easy-following codes for the main building blocks of adaptive finite element methods and multigrid methods on unstructured simplicial grids in both two and three dimensions. Besides the simplicity and readability, sparse matrixlization, an innovative programming style for MATLAB, is introduced to improve the efficiency.

<!-- - [Basic Data Structure](https://www.math.uci.edu/~chenlong/ifemdoc/mesh/meshdoc.html)
- [Finite Element Methods](https://www.math.uci.edu/~chenlong/ifemdoc/fem/femcontent.html)
- [Adaptive Finite Element Methods](https://www.math.uci.edu/~chenlong/ifemdoc/afem/afemdoc.html)
- [Solvers of Linear Algebraic Equations](https://www.math.uci.edu/~chenlong/ifemdoc/solver/solverintroduction.html)
- [Mesh Generator and Smoothing/Optimization](https://www.math.uci.edu/~chenlong/ifemdoc/mesh/meshoptdoc.html)
- [Projects](https://www.math.uci.edu/~chenlong/ifemdoc/project/projectcontent.html) -->


## Installation

You can download the repository at [https://github.com/lyc102/ifem](https://github.com/lyc102/ifem), or alternatively you can get iFEM by using the following commands:

```bash
git clone https://github.com/lyc102/ifem.git
```

Then use MATLAB to run the `setpath.m` script in the root folder to add all the sub-folders to your MATLAB path. 

An [Octave](www.octave.org) version is also available at [https://github.com/lyc102/ifemOctave](https://github.com/lyc102/ifemOctave).

## Help, Guides, and Contributing

This documentation website will be constantly updated. If you have any questions, please feel free to [contact us](mailto:lyc102@gmail.com). If you like to suggest an additional equation to be implemented in iFEM, please go to the [GitHub repo submit an issue](https://github.com/lyc102/ifem/issues). If you like to contribute to the development of iFEM, please see our [Community hub page]({{ site.baseurl }}{% link _pages/research.md %}).

### Use MATLAB help/doc
1. `help funexample` displays a description of and syntax for the function `funexample`. For example, `help mg` will show basic usage for `mg` function in the plain text.  
2. `ifem funexampledoc` show detailed description. For example, `ifem mgdoc` will explain the `mg` function step by step in `html` format. But not every function has a html documentation.

## License and References

iFEM can be freely distributed under GPL 3.0. If you feel iFEM is helpful for your research, please acknowledge your use by citing:

> L. Chen. iFEM: an integrated finite element method package in MATLAB. Technical Report, University of California at Irvine, 2009.

```bibtex
@techreport{Chen:2008ifem,
author = {Long Chen},
journal = {Technical Report, University of California at Irvine},
title = {$i$FEM: an integrated finite element methods package in {MATLAB}},
url = {https://github.com/lyc102/ifem},
year = {2009}}
```



## Acknowledgement

The author would like to thank [Professor Michael Holst](http://cam.ucsd.edu/~mholst/) in University of California at San Diego and [Professor Ludmil Zikatanov](http://www.personal.psu.edu/ltz1/) in Pennsylvania State University for many insightful discussion, and also Professor [Chensong Zhang](http://lsec.cc.ac.cn/~zhangcs/) in Chinese Academy of Sciences for the effort in the development of AFEM@matlab, an early version of iFEM.

The author thanks students or postdocs: [Shuhao Cao](https://scaomath.github.io/), [Huayi Wei](weihuayi.github.io), Ming Wang, Lin Zhong, and Jie Zhou for their contribution to iFEM in one way or another. Detailed credits can be found in the M-lint of several `.m` files.

The author is also grateful to the NSF for the partial support over the years. 



<div style="width:400px" onclick="myhref('http://math.uci.edu/~chenlong/');"><hr/>
Long Chen
<br>
Professor                 
<br>
Department of Mathematics
<br>
University of California at Irvine
http://math.uci.edu/~chenlong/
</div>

<script type="text/javascript">
    function myhref(web){
      window.location.href = web;}
</script>