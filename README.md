# Fairness In The Rashomon Set (FaiRS)

Implementation of a reduction-based algorithm to characterize the range of predictive disparities and search for the absolute predictive disparity minimizing model over the set of good models (i.e., Rashomon Set). 

If you find this repository useful for your research, please consider citing our work:
```

@InProceedings{pmlr-v139-coston21a,
  title = 	 {Characterizing Fairness Over the Set of Good Models Under Selective Labels},
  author =       {Coston, Amanda and Rambachan, Ashesh and Chouldechova, Alexandra},
  booktitle = 	 {Proceedings of the 38th International Conference on Machine Learning},
  pages = 	 {2144--2155},
  year = 	 {2021},
  editor = 	 {Meila, Marina and Zhang, Tong},
  volume = 	 {139},
  series = 	 {Proceedings of Machine Learning Research},
  month = 	 {18--24 Jul},
  publisher =    {PMLR}
}

```
Here are the [proceedings link](https://proceedings.mlr.press/v139/coston21a.html) and [arXiv link](https://arxiv.org/abs/2101.00352) to the paper.

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
