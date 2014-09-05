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
obj_name = hmmer-lyzer.HmmerAnalyze('path_to_sensor.csv', 'path_to_hxt.csv')
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

Species field can now be truncated based on whitespace. Optional field `truncate_ident` selects how many fields. Default value of 0 selects all whitespace seperated 'words'. A value of 1 selects the first, 2 selects the first couple, and so on.
The below example will select only the genus, thus allowing plots by genus instead of species. Note that all truncated names always end in an underscore.
```
obj_name.pairScores(truncate_ident=1)
```
Setting the value to 2 will allow for better species classification as strains and so forth won't affect the label any more.

Generate plots
```
obj_name.showPlot('species_id_field_entry')
```
