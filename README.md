# pyrod
Pure python package to analyse CTR (Crystal  Truncation Rod) singal 

Overview
------------
The pyrod is a pure python package to analyse CTR (Crystal  Truncation Rod) singals
Perpare intial fitting model for phase recovery algorithm

- Convenient parameters setting 
- Organised data base of x-ray diffraction
- Expansibility
- Good fitting degree (R-square > 0.9)
- Provide initial data for further phase recovery algorithm

- Note: perovskite oxide heterojunction and superlattice structure only!

Directory structure
------------
- 00L.xlsx:         example CTR data: 00L rod of LaNiO3 (12 u.c.) / SrTiO3
- parameters.xlsx:  parameters setting example
- exceptions:	    exceptions list of api
- junction:         connect XRR and CTR fitting results, generate final result
- data (dir):       optimised XRR, CTR and connection data. Serialized by pickle
- base (dir):       data base of atomic mass, atomic scattering factor and atomic form factor
- read (dir):       initialize parameters and raw data. (parameters_xxxx.xlsx and xxL_xxxxx.xlsx)
- tool (dir):	    variables control methods and built-in algorithm
- apply (dir):      api of XRR and CTR fitting

Feedback
--------
- Your comments are always welcome! We would like to expand our api for your CTR system!
- Any comments or cooperation, please contact:
																						Han Xu xuhan@mail.ustc.edu.cn

Release History
---------------------------
- beta version
