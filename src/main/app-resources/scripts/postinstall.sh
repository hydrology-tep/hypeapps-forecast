#!/bin/bash

SUDO=/usr/bin/sudo 
YUM=/usr/bin/yum
#GREP=/bin/grep
#CUT=/bin/cut

${SUDO} ${YUM} erase -y proj
# Also removes some dependencies as gdal
${SUDO} ${YUM} install -y miniconda openjpeg2 proj-4.7.0-2.el6 libgfortran ghostscript-fonts urw-fonts gdal

/opt/anaconda/bin/conda install -y --file /application/dependencies/R/packages.list
/opt/anaconda/bin/conda create --name cairo-env --file /application/dependencies/R/cairo-env.list
