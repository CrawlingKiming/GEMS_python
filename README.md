# GEMS_python

GEMS python algorithm is a framework of regrid algorithm. 
README file is written by Dong Kyu Cho. 

For experiment, contact Dong Kyu Cho or *make experiment branch* instead. 

### Installation 

For now, we do not have any dependency problem. 
However we recommend you to install python > 3.9.5 

### Descrption 

CONFIG file determines major arguments. Check "grid" setting and File IO setting. Build folder directories as you want then change CONFIG. 

To run, use python Application.py YEAR MONTH DATE


>>python Application.py 2021 03 24      



#### Caution 

There is a bug in Merge File.R 
Do not use the generated lat, lon file. 

Instead, use lat, lon from imported nc file. 
