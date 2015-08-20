# s4vd

This is a fork of the Sill-Kaiser implementation of sparse singular value biclustering
methods for R<sup>1</sup> by Lee, Shen, Huange, and Marron<sup>2</sup>.

We've tweaked the package for performance using, among other tricks, the irlba<sup>3</sup>
fast partial singular value decomposition method of Baglama and Reichel<sup>4</sup> and
the foreach package <sup>5</sup>.

More details on the changes are outlined in detail here
https://github.com/bwlewis/s4vd/blob/master/vignettes/s4vdp4.pdf


## References

1. http://cran.r-project.org/web/packages/s4vd
2. Biclustering via Sparse Singular Value Decomposition, M. Lee, H. Shen, J. Huang, J. S. Marron, Biometrics 66, pp. 1087-1095, December 2010.
3. http://cran.r-project.org/web/packages/irlba
4. Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005.
5. http://cran.r-project.org/web/packages/foreach


## Installation

```{r}
library(devtools)                  
install_github('mwsill/s4vd')
```