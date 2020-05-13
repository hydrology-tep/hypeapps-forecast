#!/bin/bash

SUDO=/usr/bin/sudo 
YUM=/usr/bin/yum
#GREP=/bin/grep
#CUT=/bin/cut

${SUDO} ${YUM} downgrade -y proj

/opt/anaconda/bin/conda install -y --file /application/dependencies/R/packages.list
/opt/anaconda/bin/conda create --name cairo-env --file /application/dependencies/R/cairo-env.list
