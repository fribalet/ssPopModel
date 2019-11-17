ssPopModel
==========
ssPopModel is a R package that uses size-structured matrix population model to estimate hourly division rates of microbial populations from SeaFlow data. These estimates are independent of variations in cell abundance and can be used to study physiologically-driven changes in population dynamics. 

The first version of the model is described in:
[Ribalet, F. et al. Light-driven synchrony of <i>Prochlorococcus</i> growth and mortality in the subtropical Pacific gyre. Proc. Natl. Acad. Sci. 112, 8008–8012 (2015)](http://www.pnas.org/lookup/doi/10.1073/pnas.1424279112)

_________
new release ssPopModel v2.0.0
(November 2019)

- Added respiration to the size-strucutured matrix population model.
- Added Huber loss as objective function and Covariance matrix adapting evolutionary strategy as optimization method (from cmaes package)
- Modified delta function by prohibiting cells less than twice the size of the smallest size to divide

Currently testing the model of the various SeaFlow cruises...
