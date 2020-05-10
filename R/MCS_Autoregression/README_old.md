# MCS in hyper-parameter optimization - an encouraging example

## Table of contents

* [Hyper-parameter optimization](#hyper-parameter-optimization)
* [The Model Confidence Set](#the-model-confidence-set)
* [Autoregressive processes](#autoregressive-processes)
* [Modelling assumptions](#modelling-assumptions)
* [Maximum likelihood estimators](#maximum-likelihood-estimators)
* [Figures](#figures)
* [Results](#results)
---

## [Hyper-parameter optimization](https://en.wikipedia.org/wiki/Hyperparameter_optimization)

The ultimate objective of a typical learning algorithm $\mathcal{A}$ is to find a function g that minimizes some expected loss
$L_{g}(x)$ over independent and identically distributed samples $x$ that has a directly not observable distribution $G_{x}$.
<br>
$\mathcal{A}$ is a functional that maps a (train) data set $X$ (finite set of samples from $G_{x}$) to a function g.
Very often a learning algorithm produces $g$ through the optimization of a training criterion with respect to a set of parameters
$\theta$. However, the learning algorithm itself has additional bells and whistles called hyper-parameters $\lambda$,
so we are only able to produce results conditionally on it. As a result $g=g_{\lambda}=\mathcal{A}_{\lambda}(X)$ for an arbitrary
training set $X$.
<br>
From now on it is a matter of taste how we choose $\lambda$, however we can define some rules of thumb. What we really need in
practice is a way to choose $\lambda$ so as to minimize generalization error $\mathbb{E}L_{\mathcal{A}_{\lambda}(X)}(x)$,
where $X$ is a data set from the distribution $G_{x}$ and $x$ is a realization of such a random variable. From the previous example
it is easily seen that for every possible $\lambda$ , $\mathcal{A}_{\lambda}$ might perform an inner optimization as well.
The problem of finding the best fitting $\lambda$ is called the problem of hyper-parameter optimization.

---
## [The Model Confidence Set](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA5771)
From statistics courses, we can remember the term confidence set, when we examined how predicted coefficients vary in a range, presumably centered at the expectation of itself (if the estimator was unbiased). <br>
The confidence set of an estimated variable is a set that contains the exact (ground truth) value at $1-\alpha$ level of confidence. That is, doing the estimation again and again, the confidence set will contain the desired value $N\cdot(1-\alpha)$ times out of $N$ experiments in expectation. It can happen, that the ratio is slightly less or greater of course.
<br>
The Model Confidence Set reflects the very same phenomena. We have a bunch of models, which we use to predict the unknown values of the underlying process. How should we decide which to choose? Of course, we observe some losses on the validation dataset, and our assumption is that the loss is zero if we choose the appropriate model. Thus choosing the one that produces the minimum loss would be a great idea.
<br>
That is exactly what the grid search method uses. However, this does not tell anything to us about the confidence of our choice! For example, what happens when the observed losses form a flat function? Almost every entry is the same, so choosing a minimal among them does not make sense, truly.
<br>
As a result, we would like to obtain a set that reacts to this issue. Whenever the losses are flat, we are not quite confident about our choice: so why not choosing the whole set as a superior set. Whenever the loss is sharp, then we can be more confident by cutting the set towards a singleton set.
<br>
The Model Confidence Set was pioneered by Hansen et al (2011) (link the title that directs to the paper). It is commonly accepted as a forecast trimmer, yet no one really tried to apply the idea to hyper-parameter optimization. With a slight refactor, I formulated the underlying mathematics to be applicable to this type of process.
<br>
 All the mathematical details are in the paper of mine, click to the project's title to see.

---
## [Autoregressive processes](https://en.wikipedia.org/wiki/Autoregressive_model)

To shed light on the greatness of the MCS, I have chosen to model an autoregression. Such a process has very interesting properties, which are
discussed in the working paper.
An autoregression is a special type of a discrete stochastic model, i.e.:
1. 1. The process follows an autoregressive trend, with parameter $p\in\mathbb{N}$, meaning, $Y_{t}=\alpha+\phi_{1}Y_{t-1}+\dots+\phi_{p}Y_{t-p}+\epsilon_{t}$,where $\forall t\in\mathbb{N},\ \epsilon_{t}\sim\mathcal{N}(0,\sigma^{2})$ and $\alpha\in\mathbb{R},\ \phi\in\mathbb{R}^{p},\ \sigma^{2}\in\mathbb{R}^{+}$ are fixed constants.
2. To have stationary results suppose that the roots of the characteristic polynomial $H(z)=1-(\phi_{1}z+\dots+\phi_{p}z^{p})$ lie within the unit circle and so, there exists a stationary solution Y_{t} of the equation. If that holds, $\underset{t\rightarrow\infty}{lim}\mathbb{E}Y_{t}=\frac{\alpha}{1-\underset{i=1}{\overset{p}{\sum}}\phi_{i}}, \underset{t\rightarrow\infty}{lim}Var(Y_{t})\vcentcolon=D^{2}$.
3. 3. For every $i\neq j$, $\mathbb{E}\epsilon_{i}\epsilon_{j}=0$, so the process is homoskedastic and the error terms are uncorrelated.
---

## [Modelling assumptions]()
Our modeling assumption are the followings:
1. There is an initial set of orders to experiment on the grid, containing the true value $p$: $\mathcal{P}=(p,p_{1},\dots,p_{n-1})$. $\mathcal{P}$ is the set of hyper-parameters.
2. For every $i\in\mathcal{P}$,we model the the process assuming the order of the autoregressive process is exactly $i$. The $(\alpha,\phi_{1},\dots,\phi_{i},\sigma^{2})^{T}$ vector is found via an inner optimization, in this case by maximizing the likelihood, that $Y$ follows an autoregressive trend with order $i$.
3. The loss function is the $L_{2}$-norm.

---
## [Maximum likelihood estimators]()

In the working paper, I derive the maximum likelihood estimators for all the parameters of an autoregression.
I conclude the results here. <br>
Let $\theta_{i}=(\phi_{i},\sigma_{i}^{2})^{T}$ be the vector of parameters to find, where $\phi_{i}=[\alpha,\phi_{1},\dots,\phi_{i}]^{T}$.
Let $n$ be the cardinality of the set of hyperparameters. Let $x_{t,i}=[1,y_{t-1},\dots,y_{t-i}]^{T}$ and $X_{i}=[x_{T,i},x_{T-1,i},\dots,x_{n+1,i}]^{T}$. Furthermore let $y=[y_{T},y_{T-1},\dots,y_{n+1}]^{T}$.
<br>
In that case $\overset{\sim}{\phi_{i}}=(X_{i}^{T}X_{i})^{-1}X_{i}^{T}y$ and $\overset{\sim}{\sigma_{i}^{2}}=\frac{(y-X_{i}\overset{\sim}{\phi_{i}})^{T}(y-X_{i}\overset{\sim}{\phi_{i}})}{T-n-1}$.
---

## [Figures]()

See the working paper.

---
## [Results]()

Results are not yet available.
