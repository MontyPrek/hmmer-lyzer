hmmer-lyzer
===========

A python module with tools to organize, modify, and visualize HMMER data.

Current Usage
=============

Import the module
```
import hmmer-lyzer
```

HmmerAnalyze is the top level class. Instantiate an object using
```
obj_name hmmer-lyzer.HmmerAnalyze('path_to_sensor.csv', 'path_to_hxt.csv')
```

Enter which columns of the CSV you'd like to use
```
obj_name.setCSVVars(sensor_score_col, sensor_species_col, sensor_gene_col,
                    hxt_score_col, hxt_species_col, hxt_gene_col)
```

Pair up scores
```
obj_name.pairScores()
```

Generate plots
```
obj_name.showPlot('species_id_field_entry')
```
