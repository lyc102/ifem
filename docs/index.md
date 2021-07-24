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

- [Basic Data Structure](./mesh/meshdoc.ipynb)
- [Finite Element Methods](./fem/femcontent.ipynb)
- [Adaptive Finite Element Methods](./afem/afemdoc.ipynb)
- [Solvers of Linear Algebraic Equations](./solver/solverintroduction.ipynb)
- [Mesh Generator and Smoothing/Optimization](./mesh/meshoptdoc.ipynb)
- [Projects](./project/projectcontent.ipynb)

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