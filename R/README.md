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

The ultimate objective of a typical learning algorithm <img src="svgs/7651ba0e8e29ee7537841a819041a172.svg?invert_in_darkmode" align=middle width=13.12555859999999pt height=22.465723500000017pt/> is to find a function g that minimizes some expected loss
<img src="svgs/4ebc527357ad84ef102eae6f59aae243.svg?invert_in_darkmode" align=middle width=41.01552839999999pt height=24.65753399999998pt/> over independent and identically distributed samples <img src="svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> that has a directly not observable distribution <img src="svgs/db64baa2c720cc025301b17b18cd2f3b.svg?invert_in_darkmode" align=middle width=20.37901634999999pt height=22.465723500000017pt/>.
<img src="svgs/7651ba0e8e29ee7537841a819041a172.svg?invert_in_darkmode" align=middle width=13.12555859999999pt height=22.465723500000017pt/> is a functional that maps a (train) data set <img src="svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/> (finite set of samples from <img src="svgs/db64baa2c720cc025301b17b18cd2f3b.svg?invert_in_darkmode" align=middle width=20.37901634999999pt height=22.465723500000017pt/>) to a function g. 
Very often a learning algorithm produces <img src="svgs/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode" align=middle width=8.430376349999989pt height=14.15524440000002pt/> through the optimization of a training criterion with respect to a set of parameters
<img src="svgs/27e556cf3caa0673ac49a8f0de3c73ca.svg?invert_in_darkmode" align=middle width=8.17352744999999pt height=22.831056599999986pt/>. However, the learning algorithm itself has additional features called hyper-parameters <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/>,
so we are only able to produce results conditionally on it. As a result <img src="svgs/7ec65d83892a59b9f70290cdb4b7cccd.svg?invert_in_darkmode" align=middle width=118.16335739999998pt height=24.65753399999998pt/> for an arbitrary
training set <img src="svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/>. <br>
From now on it is a matter of taste how we choose <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/>, however we can define some rules of thumb. What we really need in
practice is a way to choose <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> so as to minimize generalization error <img src="svgs/9202c931880bf0ea6b6d511883df3423.svg?invert_in_darkmode" align=middle width=85.38434354999998pt height=24.65753399999998pt/>, 
where <img src="svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.908688849999992pt height=22.465723500000017pt/> is a data set from the distribution <img src="svgs/db64baa2c720cc025301b17b18cd2f3b.svg?invert_in_darkmode" align=middle width=20.37901634999999pt height=22.465723500000017pt/> and <img src="svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> is a realization of such a random variable. From the previous example 
it is easily seen that for every possible <img src="svgs/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> , <img src="svgs/d27471ed880eb3004e845538e7587928.svg?invert_in_darkmode" align=middle width=20.922410849999988pt height=22.465723500000017pt/> might perform an inner optimization as well. 
The problem of finding the best fitting <img src="s/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> is called the problem of hyper-parameter optimization.

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
