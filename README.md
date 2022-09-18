# GEMS_python

> Developing grid algorithms for the Geostationary Environmental Monitoring Spectrometer

![ResultGIF](/title/Regridding_Result.gif)

GEMS python algorithm is a framework of regrid algorithm. 
README file is written by Dong Kyu Cho. 

For experiment, contact Dong Kyu Cho or *make experiment branch* instead. 

## Major Update (September)

From Satrec Initiative Team

- I updated .gitignore. Going to ignore all nc files. Do not upload any data-relevant things. 
- Now "qf_usage" works in more *efficient* way

  - Check L2converProcess.py and ql_algo method. 
  - For uncertainty quantification, you must **Modify ql_algo only**

- Grid parameter 
  
  - ea: Sparse Area 
  - kr: Dense Area (Korea penn.)

### Installation 

For now, we do not have any dependency problem. 
However we recommend you to install python > 3.9.5 

### Descrption 

CONFIG file determines major arguments. Check "grid" setting and File IO setting. Build folder directories as you want then change CONFIG. 

To run, use python Application.py YEAR MONTH DATE

```
python Application.py 2021 03 24
```
     



#### Caution 


