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
## [Autoregressive processes](https://en.wikipedia.org/wiki/Autoregressive_model)

To shed light on the greatness of the MCS, I have chosen to model an autoregression. Such a process has very interesting properties, which are
discussed in the working paper.
An autoregression is a special type of a discrete stochastic model, i.e.:
1. 1. The process follows an autoregressive trend, with parameter <img src="svgs/11e61a359d52e3bfa3ec5181155ae86a.svg?invert_in_darkmode" align=middle width=40.233883799999994pt height=22.648391699999998pt/>, meaning, <img src="svgs/d234ba5ba657a1cc39a47519278e0836.svg?invert_in_darkmode" align=middle width=258.10724939999994pt height=22.831056599999986pt/>,where <img src="svgs/0873212143ac1e5a70429ffaf2a064b3.svg?invert_in_darkmode" align=middle width=155.77262249999998pt height=26.76175259999998pt/> and <img src="svgs/f124b6dc8c466a019ebc0074531368e0.svg?invert_in_darkmode" align=middle width=176.87873609999997pt height=26.76175259999998pt/> are fixed constants.
2. To have stationary results suppose that the roots of the characteristic polynomial <img src="svgs/377c0952351dc741fc3b9e4ee16ab167.svg?invert_in_darkmode" align=middle width=217.42239419999996pt height=24.65753399999998pt/> lie within the unit circle and so, there exists a stationary solution Y_{t} of the equation. If that holds, <img src="svgs/b1aa42c312b171936c1c367a558ffd06.svg?invert_in_darkmode" align=middle width=270.63248684999996pt height=47.67709980000001pt/>.
3. 3. For every <img src="svgs/d1a177d007705974655aecd1ffe42d4c.svg?invert_in_darkmode" align=middle width=35.29127414999999pt height=22.831056599999986pt/>, <img src="svgs/2051930497302ad0342da3215a011b9b.svg?invert_in_darkmode" align=middle width=66.83976419999999pt height=22.648391699999998pt/>, so the process is homoskedastic and the error terms are uncorrelated.
---

## [Modelling assumptions]()
Our modeling assumption are the followings:
1. There is an initial set of orders to experiment on the grid, containing the true value <img src="svgs/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode" align=middle width=8.270567249999992pt height=14.15524440000002pt/>: <img src="svgs/a928308a1aa07e4be1b664615ae8405e.svg?invert_in_darkmode" align=middle width=149.2844496pt height=24.65753399999998pt/>. <img src="svgs/03c54486f0c18b8e265a9c922d83ad33.svg?invert_in_darkmode" align=middle width=12.78544904999999pt height=22.465723500000017pt/> is the set of hyper-parameters.
2. For every <img src="svgs/ec46dcb4a9298f51e7144dcb01eb8483.svg?invert_in_darkmode" align=middle width=38.53981394999999pt height=22.465723500000017pt/>,we model the the process assuming the order of the autoregressive process is exactly <img src="svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/>. The <img src="svgs/617ed95c38080f78ef0676365c36582f.svg?invert_in_darkmode" align=middle width=133.83051659999998pt height=27.6567522pt/> vector is found via an inner optimization, in this case by maximizing the likelihood, that <img src="svgs/91aac9730317276af725abd8cef04ca9.svg?invert_in_darkmode" align=middle width=13.19638649999999pt height=22.465723500000017pt/> follows an autoregressive trend with order <img src="svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/>.
3. The loss function is the <img src="svgs/394a677fe2755b93c575dfb3eeb8c276.svg?invert_in_darkmode" align=middle width=17.73978854999999pt height=22.465723500000017pt/>-norm.

---
## [Maximum likelihood estimators]()

In the working paper, I derive the maximum likelihood estimators for all the parameters of an autoregression.
I conclude the results here. <br>
Let <img src="svgs/7c731ec0dbeddc962abc9c743c64fbbf.svg?invert_in_darkmode" align=middle width=97.35703724999999pt height=27.6567522pt/> be the vector of parameters to find, where <img src="svgs/82416777084f3ab0ab11c1de26de5ea9.svg?invert_in_darkmode" align=middle width=142.69927814999997pt height=27.6567522pt/>.
Let <img src="svgs/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687624999999pt height=14.15524440000002pt/> be the cardinality of the set of hyperparameters. Let <img src="svgs/c412b371cf52c0a890baf7d01bfc506d.svg?invert_in_darkmode" align=middle width=175.82162399999999pt height=27.6567522pt/> and <img src="svgs/058501d883680bef9b6d430b1778e0aa.svg?invert_in_darkmode" align=middle width=219.78238094999995pt height=27.6567522pt/>. Furthermore let <img src="svgs/18af9f618583576c228a72d9712da005.svg?invert_in_darkmode" align=middle width=180.37616849999998pt height=27.6567522pt/>.
<br>
In that case <img src="svgs/d1bc722dfb2752bd7fd4aebcce4e2d22.svg?invert_in_darkmode" align=middle width=145.88814735pt height=38.245416pt/> and <img src="svgs/9a4564bc719354a8852c5724dc27df22.svg?invert_in_darkmode" align=middle width=163.7911341pt height=43.0457412pt/>.
---

## [Figures]()

See the working paper.

---
## [Results]()

Results are not yet available.
