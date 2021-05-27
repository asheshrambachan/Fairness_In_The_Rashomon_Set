# Fairness In The Rashomon_Set (FaiRS)

Implementation of a reduction-based algorithm to characterize the range of predictive disparities and search for the absolute predictive disparity minimizing model over the set of good models (i.e., Rashomon Set). 

If you find this repository useful for your research, please consider citing our work:
```
@misc{coston2021characterizing,
      title         = {Characterizing Fairness Over the Set of Good Models Under       
                       Selective Labels}, 
      author        = {Amanda Coston and Ashesh Rambachan and Alexandra Chouldechova},
      year          = {2021},
      eprint        = {2101.00352},
      archivePrefix = {arXiv},
      primaryClass  = {cs.LG}
}
```
This paper will appear in ICML 2021. Here is the [arXiv link](https://arxiv.org/abs/2101.00352) to the paper.

### Requirements

We provide `R` code to implement the reduction-based algorithm. The following `R` packages are used:
- `tidyverse`
- `here` 
- `ranger` (if using the random forest learner)


### Datasets
We include the datasets:
- Compas
- Communities and Crime

We cannot publicly share the data used in the consumer lending experiments.

### Usage 
The `R` markdown scripts in the directories
- `Compas_Experiments/`
- `Communities_Experiments/`
can be used to reproduce the recidivism risk prediction and regression experiments in the paper. 

The functions in `Code/ExponentiatedGradient.R` implement the reduction-based algorithms. To search for the disparity minimizing model in the set of good models, see the function `run_expgrad_extremes`. To search for the absolute disparity minimizing in the set of good models, see the function `run_expgrad_minDisp`.
