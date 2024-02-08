Cloud-Resolving Storm Simulator (CReSS) project

# Tentative announcement on retag! (Satoki Tsujino)
Ok, I messed up, and I pushed out an earlier version tagged as "v3.5.1m". I then fixed something, and retagged the fixed tree as "v3.5.1m" again. If you got the wrong tag, and want the new one, please delete the old one and fetch the new one by doing:
```
git tag -d v3.5.1m
git fetch origin tag v3.5.1m
```
to get my updated tag. You can test which tag you have by doing
```
git rev-parse v3.5.1m
```
which should return 3cb960025348.. if you have the new version. Sorry for the inconvenience.

# What is CReSS?
Cloud-Resolving Storm Simulator (CReSS) is a nonhydrostatic and regional atmosphere model for numerical simulations of tropical cyclones, thunderstorms, tornados, and other severe weather phenomena. 

The CReSS model is implemented in the Fortran 90 language.
For parallel computing, multithread parallelization is supprted by OpenMP, and process parallelization is supported by MPI.

[More information](http://www.rain.hyarc.nagoya-u.ac.jp/%7Etsuboki/kibanS2/src_eng/cress_synopsis_eng.html)

# How to use
Please see [Doc/0rig/readme_first.txt](https://cress-nagoya.github.io/CReSS/Doc/0rig/readme_first.txt). 

# Examples
You can find examples of the configuration and setting file for the CReSS model simulation in `Form/0rig/`. 

# Documents
* [Official English Document (Old version)](http://www.rain.hyarc.nagoya-u.ac.jp/~tsuboki/cress_html/src_cress/CReSS2223_users_guide_eng.pdf)
* [Japanese User's Guide](http://www.rain.hyarc.nagoya-u.ac.jp/~tsuboki/cress_html/from_kato/how_to_use_cress_20110413.pdf)

# References
1. Tsuboki, K., 2023: High-Resolution Simulations of Tropical Cyclones and Mesoscale Convective Systems Using the CReSS Model. Park, S.K. (Eds), _Numerical Weather Prediction: East Asian Perspectives._ Springer Atmospheric Sciences. Springer, Cham, 483-534. (https://doi.org/10.1007/978-3-031-40567-9_19)
2. Tsuboki, K., 2008: High-Resolution Simulations of High-Impact Weather Systems Using the Cloud-Resolving Model on the Earth Simulator. Hamilton, Kevin; Ohfuchi, Wataru (Eds.), _High Resolution Numerical Modelling of the Atmosphere and Ocean_, Springer New York, 141-156.
3. Tsuboki, K. and A. Sakakibara, 2007: Numerical prediction of high-impact weather systems. _The Textbook for Seventeenth IHP Training Course in 2007_, 281pp.
4. Tsuboki, K. and A. Sakakibara, 2002: Large scale parallel computing of Cloud Resolving Storm Simulator. H. P. Zima et al. (Eds), _High Performance Computing_, Springer, 243-259.
