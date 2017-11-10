# HESS

HESS (Heritability Estimation from Summary Statistics) is a software package
for estimating and visualizing local SNP-heritability and genetic covariance
(correlation) from GWAS summary association data.

## Software requirement

### For running HESS

* [Python 2.7](https://www.python.org/download/releases/2.7/)
* [NumPy 1.9+](http://www.numpy.org/)
* [Pandas 0.18+](http://pandas.pydata.org/)
* [PySnpTools 0.3+](https://github.com/MicrosoftGenomics/PySnpTools)
* [SciPy 0.16+](https://www.scipy.org/)

### For running analysis and visualization tools

* [Python 2.7](https://www.python.org/download/releases/2.7/)
* [NumPy 1.9+](http://www.numpy.org/)
* [Pandas 0.18+](http://pandas.pydata.org/)
* [Matplotlib 1.4+](https://matplotlib.org/)
* [SciPy 0.16+](https://www.scipy.org/)

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: It is important to use the latest version of NumPy, which is
more numerically stable when performing eigendecomposition of the LD
matrix. Python 3.0+ may not work well ith HESS, which is written and tested
in Python 2.7.3. )
</div>

## References

* **Local SNP-heritability**: [Shi, Huwenbo, Gleb Kichaev, and Bogdan Pasaniuc. "Contrasting the genetic architecture of 30 complex traits from summary association data." The American Journal of Human Genetics 99, no. 1 (2016): 139-153.](http://www.sciencedirect.com/science/article/pii/S0002929716301483)
* **Local genetic covariance (correlation)**: [Shi, Huwenbo, Nicholas Mancuso, Sarah Spendlove, and Bogdan Pasaniuc. "Local genetic correlation gives insights into the shared genetic architecture of complex traits." The American Journal of Human Genetics 101, no. 5 (2017): 737-751.](http://www.sciencedirect.com/science/article/pii/S0002929717303919)

## Contact

* **Huwenbo Shi**: shihuwenbo_AT_ucla.edu
* **Bogdan Pasaniuc**: bpasaniuc_AT_mednet.ucla.edu 
