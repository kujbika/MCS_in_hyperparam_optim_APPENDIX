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

The ultimate objective of a typical learning algorithm <img src="svgs/7651ba0e8e29ee7537841a819041a172.svg?invert_in_darkmode" align=middle width=13.12555859999999pt height=22.465723500000017pt/> is to find a function g that minimizes some expected loss
<img src="svgs/4ebc527357ad84ef102eae6f59aae243.svg?invert_in_darkmode" align=middle width=41.01552839999999pt height=24.65753399999998pt/> over independent and identically distributed samples <img src="svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> that has a directly not observable distribution <img src="svgs/db64baa2c720cc025301b17b18cd2f3b.svg?invert_in_darkmode" align=middle width=20.37901634999999pt height=22.465723500000017pt/>.
<br>
<img src="svgs/7651ba0e8e29ee7537841a819041a172.svg?invert_in_darkmode" align=middle width=13.12555859999999pt height=22.465723500000017pt/> is a functional that maps a (train) data set <img src="svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/> (finite set of samples from <img src="svgs/db64baa2c720cc025301b17b18cd2f3b.svg?invert_in_darkmode" align=middle width=20.37901634999999pt height=22.465723500000017pt/>) to a function g.
Very often a learning algorithm produces <img src="svgs/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode" align=middle width=8.430376349999989pt height=14.15524440000002pt/> through the optimization of a training criterion with respect to a set of parameters
<img src="svgs/27e556cf3caa0673ac49a8f0de3c73ca.svg?invert_in_darkmode" align=middle width=8.17352744999999pt height=22.831056599999986pt/>. However, the learning algorithm itself has additional bells and whistles called hyper-parameters <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/>,
so we are only able to produce results conditionally on it. As a result <img src="svgs/7ec65d83892a59b9f70290cdb4b7cccd.svg?invert_in_darkmode" align=middle width=118.16335739999998pt height=24.65753399999998pt/> for an arbitrary
training set <img src="svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/>.
<br>
From now on it is a matter of taste how we choose <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/>, however we can define some rules of thumb. What we really need in
practice is a way to choose <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> so as to minimize generalization error <img src="svgs/9202c931880bf0ea6b6d511883df3423.svg?invert_in_darkmode" align=middle width=85.38434354999998pt height=24.65753399999998pt/>,
where <img src="svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/> is a data set from the distribution <img src="svgs/db64baa2c720cc025301b17b18cd2f3b.svg?invert_in_darkmode" align=middle width=20.37901634999999pt height=22.465723500000017pt/> and <img src="svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> is a realization of such a random variable. From the previous example
it is easily seen that for every possible <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> , <img src="svgs/d27471ed880eb3004e845538e7587928.svg?invert_in_darkmode" align=middle width=20.922410849999988pt height=22.465723500000017pt/> might perform an inner optimization as well.
The problem of finding the best fitting <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> is called the problem of hyper-parameter optimization.

---
## [The Model Confidence Set](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA5771)
From statistics courses, we can remember the term confidence set, when we examined how predicted coefficients vary in a range, presumably centered at the expectation of itself (if the estimator was unbiased). <br>
The confidence set of an estimated variable is a set that contains the exact (ground truth) value at <img src="svgs/900d920909b4893be83c15f105c8ae1c.svg?invert_in_darkmode" align=middle width=38.88690464999999pt height=21.18721440000001pt/> level of confidence. That is, doing the estimation again and again, the confidence set will contain the desired value <img src="svgs/491e875b4a4bb3841311497d528e1332.svg?invert_in_darkmode" align=middle width=78.54428174999998pt height=24.65753399999998pt/> times out of <img src="svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> experiments in expectation. It can happen, that the ratio is slightly less or greater of course.
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
Lloyd (2001) propose that a realistic or empirically provable distribution can be obtained by choosing <img src="svgs/bb7a14d80e3cf63b2aa80d2c30c1687a.svg?invert_in_darkmode" align=middle width=26.99209754999999pt height=24.65753399999998pt/> to be a gamma probability density function, with parameters \gamma and <img src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687624999999pt height=14.15524440000002pt/> (<img src="svgs/8cda31ed38c6d59d14ebefa440099572.svg?invert_in_darkmode" align=middle width=9.98290094999999pt height=14.15524440000002pt/> and <img src="svgs/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode" align=middle width=14.433101099999991pt height=14.15524440000002pt/> for the exposed class). The expectation of such a variable is <img src="svgs/662cba0b94927c67c7bd63d950d18316.svg?invert_in_darkmode" align=middle width=34.346253149999995pt height=27.77565449999998pt/>, which fulfills our requirements described in the previous subsection. However, using the fact that the gamma distribution is a sum of exponentially distributed variables we introduce <img src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687624999999pt height=14.15524440000002pt/> (<img src="svgs/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode" align=middle width=14.433101099999991pt height=14.15524440000002pt/>) different infectious (exposed) class <img src="svgs/856f8555125e8ddd238cd12237ef4ceb.svg?invert_in_darkmode" align=middle width=103.21722344999999pt height=29.190975000000005pt/> (<img src="svgs/0036764651a93d9ff7dec55bbaf6ebfc.svg?invert_in_darkmode" align=middle width=115.88847269999998pt height=29.190975000000005pt/>), with <img src="svgs/fed7181182cb5881b6953d1ab10a3173.svg?invert_in_darkmode" align=middle width=19.290757199999987pt height=14.15524440000002pt/> (<img src="svgs/49345c53ba3b754e10d9dec135705e32.svg?invert_in_darkmode" align=middle width=24.41600204999999pt height=14.15524440000002pt/>) being the rate of sequential progression through the sub classes. Equivalently, the time spent in each infectious class is exponentially distributed with <img src="svgs/98dd3bf2c169083eb5ad9b4ae87d6d70.svg?invert_in_darkmode" align=middle width=15.738769199999997pt height=27.77565449999998pt/> (<img src="svgs/05aa60dc24533f565795beccabbc4b07.svg?invert_in_darkmode" align=middle width=19.695159pt height=27.77565449999998pt/>) expectation, so that the whole infectious period is <img src="svgs/52daf868b01428da25531e4b5f186c18.svg?invert_in_darkmode" align=middle width=7.612745249999997pt height=27.77565449999998pt/> (<img src="svgs/abf732cc96ea8814962ed6d287059a11.svg?invert_in_darkmode" align=middle width=8.030309099999997pt height=27.77565449999998pt/>) in expectation, and it follows the proposed gamma distribution.


To incorporate contact tracing and isolation into the model, Wearing et al. (2005) proposed another extension of the model, while still dealing with gamma-distributed latent and infectious periods. In this model, isolation of newly infectious cases occurs at a daily rate <img src="svgs/6fb3ad53b9c3e50139eac4c41d03389d.svg?invert_in_darkmode" align=middle width=15.27633194999999pt height=22.831056599999986pt/> after a delay of <img src="svgs/6834db793eb4d3eefb92e28742316579.svg?invert_in_darkmode" align=middle width=18.28822049999999pt height=14.15524440000002pt/> days, which represents a period when infected individuals are infectious but asymptomatic or undetectable (<img src="svgs/feb77772a5a378cdcd84637b317d36f3.svg?invert_in_darkmode" align=middle width=17.11196189999999pt height=22.465723500000017pt/>).
<br>
A fraction of <img src="svgs/d5c18a8ca1894fd3a7d25f242cbe8890.svg?invert_in_darkmode" align=middle width=7.928106449999989pt height=14.15524440000002pt/> of those who had contact with an infectious and symptomatic person (<img src="svgs/b0287aa6bdc305956f5b6eeb1d0ef2e5.svg?invert_in_darkmode" align=middle width=15.92700119999999pt height=22.465723500000017pt/>) (but did not contract the infection) are remove to the quarantined susceptible class. An identical fraction of newly exposed individuals are also quarantined.
<br>
<img src="svgs/MN_SEIR.png?invert_in_darkmode">

---
## [MN-SEIR main equations]()
The following equations were introduced by Wearing et al. (2005), however I made a slight change in the formulation. To catch that people does not contact as often as before, as the epidemic increases, I assumed that the contact number decreases to a minimum of <img src="svgs/575209ad9660c13e60763153293a4a53.svg?invert_in_darkmode" align=middle width=15.11042279999999pt height=22.831056599999986pt/>. The minimal contact size is meaningful, because for example people has to do some shopping even in the middle of an epidemic. The rate how <img src="svgs/dc1df253dd3e063c495a9d1f19c5ddbf.svg?invert_in_darkmode" align=middle width=27.796893299999986pt height=24.65753399999998pt/> shrinks is <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/>.

<img src="svgs/39402bed7e8ecd2fada4840fece66465.svg?invert_in_darkmode" align=middle width=98.9410521pt height=33.20539859999999pt/>

<img src="svgs/df9cc268a7dd2dd12f31129fa89b6e4f.svg?invert_in_darkmode" align=middle width=431.16313514999996pt height=33.95427420000001pt/>

<img src="svgs/e4e690ede85203967196f612296eaaf5.svg?invert_in_darkmode" align=middle width=356.85353715pt height=33.95427420000001pt/>

<img src="svgs/d6353363342f62758d8060bb667b308c.svg?invert_in_darkmode" align=middle width=263.99810579999996pt height=33.20539859999999pt/>

<img src="svgs/a5b0b35a0dc2e3a19f512f13a9a5b141.svg?invert_in_darkmode" align=middle width=306.28709595pt height=29.46111299999998pt/>

<img src="svgs/053f7e10f12ae0a73c090a9abcdad8f1.svg?invert_in_darkmode" align=middle width=276.87637889999996pt height=32.675741999999985pt/>

<img src="svgs/78ea738232a9601c999cf9d65456950c.svg?invert_in_darkmode" align=middle width=383.15515919999996pt height=32.675741999999985pt/>

<img src="svgs/03b5110b20b619a138f700ef0db1f4e0.svg?invert_in_darkmode" align=middle width=234.47583555pt height=32.675741999999985pt/>

<img src="svgs/f959f5eed91000c6017765da6236e6c2.svg?invert_in_darkmode" align=middle width=426.71209680000004pt height=32.675741999999985pt/>

<img src="svgs/2e2829344b640b819d38ae7612f2d3d5.svg?invert_in_darkmode" align=middle width=204.698109pt height=33.20539859999999pt/>

<img src="svgs/ea42945db7f51b368c4a97f122aeeaa3.svg?invert_in_darkmode" align=middle width=117.82386329999999pt height=28.92634470000001pt/>

<img src="svgs/e6fe29272496d9a4fca0aea9840392c3.svg?invert_in_darkmode" align=middle width=245.84684850000002pt height=28.92634470000001pt/>

where

<img src="svgs/bdf08f00ce9bff2f86a2e7629dffddd1.svg?invert_in_darkmode" align=middle width=289.01349675pt height=36.825666900000016pt/>

is the expected number of infectious individuals at time <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936097749999991pt height=20.221802699999984pt/> that are still in infectious class <img src="svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> after a fixed delay of <img src="svgs/6834db793eb4d3eefb92e28742316579.svg?invert_in_darkmode" align=middle width=18.28822049999999pt height=14.15524440000002pt/> days.



---
## [Parameters, hyper-parameters]()

The parameters of the models (within each group <img src="svgs/924c2a38ef139efbe6801016f51628cd.svg?invert_in_darkmode" align=middle width=31.605860549999992pt height=14.15524440000002pt/>) are:
* <img src="svgs/ce028a0a004b805d0985a960c34b60e4.svg?invert_in_darkmode" align=middle width=20.513758649999993pt height=22.465723500000017pt/>
* <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/>
* b
* q
* <img src="svgs/8cda31ed38c6d59d14ebefa440099572.svg?invert_in_darkmode" align=middle width=9.98290094999999pt height=14.15524440000002pt/>
* <img src="svgs/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode" align=middle width=9.423880949999988pt height=14.15524440000002pt/>
* <img src="svgs/917b0b5c3b69af0b1d9daa049c4443da.svg?invert_in_darkmode" align=middle width=17.53863374999999pt height=14.15524440000002pt/>
* <img src="svgs/6834db793eb4d3eefb92e28742316579.svg?invert_in_darkmode" align=middle width=18.28822049999999pt height=14.15524440000002pt/>
* <img src="svgs/6fb3ad53b9c3e50139eac4c41d03389d.svg?invert_in_darkmode" align=middle width=15.27633194999999pt height=22.831056599999986pt/>
* N
From this list, <img src="svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> and <img src="svgs/917b0b5c3b69af0b1d9daa049c4443da.svg?invert_in_darkmode" align=middle width=17.53863374999999pt height=14.15524440000002pt/> can be handled as constants: <img src="svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> is the population size and <img src="svgs/917b0b5c3b69af0b1d9daa049c4443da.svg?invert_in_darkmode" align=middle width=17.53863374999999pt height=14.15524440000002pt/> is the susceptible quarantine period which is fixed by the government. The only hyper-parameters of the system are <img src="svgs/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode" align=middle width=14.433101099999991pt height=14.15524440000002pt/> and <img src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687624999999pt height=14.15524440000002pt/>.

---

## [Figures]()

See the working paper.

---
## [Results]()

Results are not yet available.
