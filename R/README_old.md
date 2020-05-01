# MCS in hyper-parameter optimization - in light of predicting Covid-19 cases

## Table of contents

* [Hyper-parameter optimization](#hyper-parameter-optimization)
* [The Model Confidence Set](#the-model-confidence-set)
* [The MN-SEIR model](#the-mn-seir-model)
* [Parameters, hyper-parameters](#parameters-,-hyper-parameters)
* [Figures](#figures)
* [Results](#results)
---

## [Hyper-parameter optimization]()

The ultimate objective of a typical learning algorithm $\mathcal{A}$ is to find a function g that minimizes some expected loss
$L_{g}(x)$ over independent and identically distributed samples $x$ that has a directly not observable distribution $G_{x}$.
$\mathcal{A}$ is a functional that maps a (train) data set $X$ (finite set of samples from $G_{x}$) to a function g. 
Very often a learning algorithm produces $g$ through the optimization of a training criterion with respect to a set of parameters
$\theta$. However, the learning algorithm itself has additional features called hyper-parameters $\lambda$,
so we are only able to produce results conditionally on it. As a result $g=g_{\lambda}=\mathcal{A}_{\lambda}(X)$ for an arbitrary
training set $X$. <br>
From now on it is a matter of taste how we choose $\lambda$, however we can define some rules of thumb. What we really need in
practice is a way to choose $\lambda$ so as to minimize generalization error $\mathbb{E}L_{\mathcal{A}_{\lambda}(X)}(x)$, 
where $X$ is a data set from the distribution $G_{x}$ and $x$ is a realization of such a random variable. From the previous example 
it is easily seen that for every possible $\lambda$ , $\mathcal{A}_{\lambda}$ might perform an inner optimization as well. 
The problem of finding the best fitting $\lambda$ is called the problem of hyper-parameter optimization.

---
## [The Model Confidence Set](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA5771)
aasdasdasd

---
## [The MN-SEIR model](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0020174)

asdadas

---
## [Parameters, hyper-parameters]()

asdasd

---

## [Figures]()

asdasda

---
## [Results]()

asda
