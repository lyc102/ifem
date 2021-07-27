---
permalink: /projects/fmm/
title: "Project: Fast Multipole Methods"
sidebar:
    nav: project
---

The purpose of this project is to implement the tree code and fast multipole
methods for the N-body summation problem. Given $\boldsymbol x = (x_j)\in \mathbb R^N, \boldsymbol y = (y_i)\in \mathbb R^N$, let 

$$A_{i,j} = \frac{1}{\|x_j - y_i\|^2}.$$ 

For a given vector $\boldsymbol q \in \mathbb R^N$, we will compute the matrix-vector product

$$\boldsymbol u = A \boldsymbol q$$ 

in $\mathcal O(N\log N)$ and $\mathcal O(N)$ operations.

Reference: 

[Fast Multipole Methods](http://math.uci.edu/~chenlong/226/FMMsimple.pdf)

## Step 1: Direct Sum

Generate two random vectors `x, y` with length `N`. Although the direct sum
can be implemented in the double for loops, in MATLAB, it is better to
generate the matrix first and then compute the matrix-vector product for
another random vector `q`.

- Use double `for` loops to generate the matrix `A`
- Use vector product to generat `A` in one line

## Step 2: Compute the weight

- Loop over cells `for i=1:N` to compute the weight

- Try to remove the for loop using vectorization. 

The loop over levels is small (only $log N$ times) and thus can be kept. To store the weight in different levels, use `cell` structure. 

## Step 3: Evaluation procedure

- Find the interaction list. 
- Loop over each cell in a given level and compute the far field in the interaction list. 
- In the fines level, add the near field by direct sum or matrix-vector product using a small matrix.

## Step 4: Test

- Choose small `N` and `J = 1`. Make sure the code works for one level (only four intervals) first by comparing the result using tree algorithm with the result in Step 1.

- Test the performance for different `N` and plot the CPU time vs N for both direct method and your tree code.

## Step 5 (optional): Fast Multipole Methods

Modify the tree code to fast multipole methods.

- Compute the weight by restriction from the fine grid to coarse grid.

- Implement the `M2L`: multipole expansion to local expansion

- Change the evaluation of far field in the interaction list to the merge of coefficients `b` in the local expansion.

- Translate the local expansion using the prolongation operator.

- Evaluate in the finest level.

- Plot the CPU time vs N to confirm the `O(N)` complexity.
