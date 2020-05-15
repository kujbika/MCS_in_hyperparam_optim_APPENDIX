# Model Confidence Set in hyper-parameter optimization - online appendix
This is the code repository marked as the online appendix for the work: Marcell Kujbus - A Model Confidence Set extension of grid-search in hyper-parameter optimization
<br>

---
**Abstract**
<br>
Abstract Hyper-parameter optimization is a thoroughly investigated discipline in optimization theory, statistics and machine learning. While it is crucial to configure the hyper-parameters properly, the most widely used grid search technique has some hidden issues that impacts the efficiency of the learning algorithm. Such an issue is that the mentioned method does not take statistical significance of out-of-sample forecasting performance into consideration when choosing the hyper-parameter yielding the minimal loss. The following paper proposes an extension of this approach to make the modeling more robust to random fluctuation. The Model Confidence Set algorithm terminates at a superior set of hyper-parameters that are statistically indistinguishable at a given a significance level with respect to the loss their generated models produce. It is shown empirically and theoretically, that averaging over the predictions of such a superior set is more efficient with respect to out-of-sample loss than taking the prediction of a single realization of the ”minimum-yielding” hyper-parameter. The theoretical argument is based on showing a specific example, where the limitations of the grid search method are reached. The empirical study examines the predictability of the evolution of the novel Covid-19 infectious disease in New York state, U.S.A. Extending the famous SEIR dynamic system with gamma distributed latent and infectious periods gives rise to hyper-parameter optimization in disease prediction. 
<br>
<br>
**Keywords:** global optimization, hyper-parameter optimization, model selection, Model Confidence Set, autoregressive processes, maximum likelihood estimator, infectious disease models, SEIR model, gamma distribution
The repository is operated and maintained by the author (marcellkujbus@gmail.com).
<br><br>

---
<b> The working paper can be found here: </b><br>
https://marcellkujbus.wixsite.com/mysite/projects
<br>
For contribution, see the [contribution guideline](CONTRIBUTION.md)
