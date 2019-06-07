#!/bin/bash

/opt/anaconda/bin/conda install -y --file /application/dependencies/R/packages.list

/opt/anaconda/bin/conda create --name cairo-env --file /application/dependencies/R/cairo-env.list
/opt/anaconda/bin/conda create --name cdo-env --file /application/dependencies/R/cdo-env.list
