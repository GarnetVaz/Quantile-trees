Quantile-trees
==============

R package implementation of decision trees using the absolute value loss function. (Supports LAD trees and other quantiles to extend linear Quantile Regression (http://en.wikipedia.org/wiki/Quantile_regression) to the tree framework)

This package provides an efficient implementation of LAD trees and the general case
of any quantile using the tilted absolute value function. This current implementation
has the same complexity as that of the 'tree' package in R.

Run times of tree vs qtree(tau=0.5) on CT-slices dataset (384 variables, default settings):

N:

 [1]    32    64   128   256   512  1024  2048  4096  8192 16384 32768

Average tree timings (secs):

 [1] 0.004567779 0.004790237 0.005605394 0.007774756 0.011044069 0.018317890
 
 [7] 0.034306114 0.068284928 0.137372828 0.290143928 0.618796518
 
Average qtree timings (secs):

 [1] 0.004271837 0.006552471 0.008399708 0.012150929 0.020915846 0.038780169
 
 [7] 0.074778188 0.156235678 0.308181793 0.629097475 1.299138459

The package relies on the Armadillo C++ library installed through 'RcppArmadillo'. The structure
is very similar to the 'tree' class and hence the same predict function is used.

Todo: 
* Categorical variables handling.
* Surrogate variables computation.


