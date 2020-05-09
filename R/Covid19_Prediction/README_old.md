# MCS in hyper-parameter optimization - in light of predicting Covid-19 cases

## Table of contents

* [Hyper-parameter optimization](#hyper-parameter-optimization)
* [The Model Confidence Set](#the-model-confidence-set)
* [The MN-SEIR model](#the-mn-seir-model)
* [MN-SEIR main equations](#mn-seir-main-equations)
* [Parameters, hyper-parameters](#parameters-,-hyper-parameters)
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
## [The MN-SEIR model](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0020174)

The most commonly used framework for epidemiological systems is the susceptible-infectious-recovered (SIR) class of models, in which the host population is categorized according to infection status as either susceptible, infectious or recovered, respectively (Kermack and McKendrick (1927)). Subsequent refinements of the model have incorporated an additional exposed (infected but not yet infectious) class so that the model became SEIR.
<br>
One of the fundamental mathematical assumptions in such models is that the rate of leaving the exposed and the infectious class is constant, irrespective of the period already spent in that class. While this makes the mathematics very easy to handle, this assumption gives rise to exponentially distributed latent and infectious periods, which is epidemiologically unrealistic for most infections (Sartwell et al. (1950), Bailey (1975), Wearing et al. (2005)).
<br>
Lloyd (2001) propose that a realistic or empirically provable distribution can be obtained by choosing $p(t)$ to be a gamma probability density function, with parameters \gamma and $n$ ($\sigma$ and $m$ for the exposed class). The expectation of such a variable is $\frac{1}{\gamma}(\frac{1}{\sigma})$, which fulfills our requirements described in the previous subsection. However, using the fact that the gamma distribution is a sum of exponentially distributed variables we introduce $n$ ($m$) different infectious (exposed) class $(I^{(1)},\dots,I^{(n)})$ ($(E^{(1)},\ldots,E^{(m)})$), with $n\gamma$ ($m\sigma$) being the rate of sequential progression through the sub classes. Equivalently, the time spent in each infectious class is exponentially distributed with $\frac{1}{n\gamma}$ ($\frac{1}{m\sigma}$) expectation, so that the whole infectious period is $\frac{1}{\gamma}$ ($\frac{1}{\sigma}$) in expectation, and it follows the proposed gamma distribution.


To incorporate contact tracing and isolation into the model, Wearing et al. (2005) proposed another extension of the model, while still dealing with gamma-distributed latent and infectious periods. In this model, isolation of newly infectious cases occurs at a daily rate $d_{I}$ after a delay of $\tau_{D}$ days, which represents a period when infected individuals are infectious but asymptomatic or undetectable ($I_{A}$).
<br>
A fraction of $q$ of those who had contact with an infectious and symptomatic person ($I_{S}$) (but did not contract the infection) are remove to the quarantined susceptible class. An identical fraction of newly exposed individuals are also quarantined.

---
## [MN-SEIR main equations]()
The following equations were introduced by Wearing et al. (2005), however I made a slight change in the formulation. To catch that people does not contact as often as before, as the epidemic increases, I assumed that the contact number decreases to a minimum of $k_{0}$. The minimal contact size is meaningful, because for example people has to do some shopping even in the middle of an epidemic. The rate how $k(t)$ shrinks is $\lambda$.

$\frac{dk}{dt}=-\frac{k(t)-k_{0}}{\lambda}$

$\frac{dS}{dt}=-\frac{(k(t)bI(t)+qk(t)(1-b)I_{S}(t))S(t)}{N}+\frac{qk(t)(1-b)S(t-\tau_{Q})I_{S}(t-\tau_{Q})}{N}$

$\frac{dS_{Q}}{dt}=\frac{qk(t)(1-b)S(t)I_{S}(t)}{N}-\frac{qk(t)(1-b)S(t-\tau_{Q})I_{S}(t-\tau_{Q})}{N}$

$\frac{dE_{1}}{dt}=\frac{k(t)b(I(t)-qI_{S}(t))S(t)}{N}-m\sigma E_{1}(t)$

$\frac{dE_{i}}{dt}=m\sigma E_{i-1}(t)-m\sigma E_{i}(t),\quad i=2,...,m$

$\frac{dI_{A,1}}{dt}=m\sigma E_{m}(t)-n\gamma I_{A,1}(t)-P_{I,1}(t)$

$\frac{dI_{A,i}}{dt}=n\gamma I_{A,i-1}(t)-n\gamma I_{A,i}(t)-P_{I,i}(t),\quad i=2,...,n$

$\frac{dI_{S,1}}{dt}=P_{I,1}(t)-(n\gamma+d_{I})I_{S,1}(t)$

$\frac{dI_{S,i}}{dt}=P_{I,i}(t)+n\gamma I_{S,i-1}(t)-(n\gamma+d_{I})I_{S,i}(t),\quad i=2,...,n$

$\frac{dQ}{dt}=\frac{qk(t)bS(t)I_{S}(t)}{N}+d_{I}I_{S}(t)$

$\frac{dD}{dt}=n\gamma\mu I_{S,n}(t)$

$\frac{dR}{dt}=n\gamma I_{A,n}(t)+n\gamma(1-\mu)I_{S,n}(t)$

where

$P_{I,i}(t)=m\sigma E_{m}(t-\tau_{D})e^{-n\gamma\tau_{D}}\frac{(n\gamma\tau_{D})^{i-1}}{(i-1)!}$

is the expected number of infectious individuals at time $t$ that are still in infectious class $i$ after a fixed delay of $\tau_{D}$ days.



---
## [Parameters, hyper-parameters]()

The parameters of the models (within each group $m,n$) are:
* $K_{0}$
* $\lambda$
* b
* q
* $\sigma$
* $\gamma$
* $\tau_{Q}$
* $\tau_{D}$
* $d_{I}$
* N
From this list, $N$ and $\tau_{Q}$ can be handled as constants: $N$ is the population size and $\tau_{Q}$ is the susceptible quarantine period which is fixed by the government. The only hyper-parameters of the system are $m$ and $n$.

---

## [Figures]()

See the working paper.

---
## [Results]()

Results are not yet available.
