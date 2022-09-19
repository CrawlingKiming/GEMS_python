# GEMS_python

> Developing grid algorithms for the Geostationary Environmental Monitoring Spectrometer

![ResultGIF](/title/Regridding_Result.gif)

GEMS python algorithm is a framework of regrid algorithm. 
README file is written by Dong Kyu Cho. 

For experiment, contact Dong Kyu Cho or *make experiment branch* instead. 

## Major Update (September)

From Satrec Initiative Team

- I updated .gitignore. It will prevent uploading nc files. 
- Now "qf_usage" works in more *efficient* way

  - Check L2converProcess.py and ql_algo method. 
  - For uncertainty quantification, you must **Modify ql_algo only**
  - Check CONFIG "qf_usage"

- Grid parameter 
  
  - ea: Sparse Area 
  - kr: Dense Area (Korea penn.)

- Class Modified in more efficient way. Now, you should pass qf parameter. 

### Installation 

For now, we do not have any dependency problem. 
However we recommend you to install python > 3.9.5 

### Descrption 

CONFIG file determines major arguments. Check "grid" setting and File IO setting. Build folder directories as you want then change CONFIG. 

To run, use python Application.py YEAR MONTH DATE TIME

For example, to run a 202103240445 nc file, 
```
python Application.py 2021 03 24 04
```

You will have two outputs: kr and ea.       



#### Caution 


