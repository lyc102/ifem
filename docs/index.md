---
layout: single
permalink: /
title: "iFEM: an Integrated Finite Element Methods Package in MATLAB"
excerpt: "A quick start guide."
sidebar:
    nav: docs
---

## Introduction to iFEM
*i*FEM is a MATLAB software package containing robust, efficient, and easy-following codes 
for the main building blocks of adaptive finite element methods and multigrid methods 
on unstructured simplicial grids in both two and three dimensions. Besides the simplicity and readability, sparse matrixlization, an innovative programming style for MATLAB, is introduced to improve the efficiency.

- [Basic Data Structure](../ifemdoc/mesh/meshdoc.ipynb)
- [Finite Element Methods](../ifemdoc/fem/femcontent.ipynb)
- [Adaptive Finite Element Methods](../ifemdoc/afem/afemdoc.ipynb)
- [Solvers of Linear Algebraic Equations](../ifemdoc/solver/solverintroduction.ipynb)
- [Mesh Generator and Smoothing/Optimization](../ifemdoc/mesh/meshoptdoc.ipynb)
- [Projects](../ifemdoc/project/projectcontent.ipynb)



## References

If you like to suggestion an additional equation to be implemented in iFEM, please go to the GitHub repo submit a feature request issue. If you feel it is helpful for your research, please acknowledge your use by citing:

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

The author would like to thank Professor Michael Holst in University of California at San Diego and Professor Ludmil Zikatanov in Pennsylvania State University for many insightful discussion, and also Professor Chensong Zhang in Chinese Academy of Sciences for the effort in the development of AFEM@matlab, an early version of iFEM.

The author thanks students or postdocs: Shuhao Cao, Huayi Wei, Ming Wang, Lin Zhong, and Jie Zhou for their contribution to iFEM in one way or another. Detailed credits can be found in the M-lint of several m files.

The author is also grateful to the NSF for the partial support over the years. 




Long Chen

--

Professor                 

Department of Mathematics
University of California at Irvine
http://math.uci.edu/~chenlong/

--