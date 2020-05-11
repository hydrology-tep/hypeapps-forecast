#!/bin/bash

#
# Add R libraries/packages and other necessary SW components to the develop and run-time environment.
# Prefereably, install the main part of the major SW components first via yum to reduce any
# incompability that may arise when conda upgrades/downgrades or installs additional library
# dependecies. And not all libraries/packages are available from the conda channels.
# If more SW components have been added, also add the installation dependency to pom.xml and the
# user guide for setting up a sandbox environment.
#

# yum list installed
# miniconda.x86_64                     4.6.14-1.el6               @ciop-stable
# openjpeg2.x86_64                     2.3.0-6.el6                @epel
# proj.x86_64                          4.7.0-2.el6                @ciop-stable
# libgfortran.x86_64                   4.4.7-23.el6               @base
# ghostscript-fonts.noarch             5.50-23.2.el6              @base
# urw-fonts.noarch                     2.4-11.el6                 @base

export PATH=/opt/anaconda/bin/:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:$PATH

logfile=/home/tobiaslagander/log-mvn.txt

SUDO=/usr/bin/sudo 
YUM=/usr/bin/yum
GREP=/bin/grep
CUT=/bin/cut

${SUDO} ${YUM} downgrade -y proj
echo $? > $logfile

echo $PATH >> $logfile

# exp_proj_ver="4.7.0-2.el6"
# #exp_proj_ver="4.8.0-2.rhel6"

# proj_ver=$(${YUM} list installed | ${GREP} -i proj | ${CUT} -d ' ' -f 27)

# if [ -n "${proj_ver}" ]; then

#     if [ "x${proj_ver}" != "x${exp_proj_ver}" ]; then
# 	#${SUDO} ${YUM} install -y proj-${exp_proj_ver}
# 	${SUDO} ${YUM} downgrade -y proj
#     fi
# fi

/opt/anaconda/bin/conda install -y --file /application/dependencies/R/packages.list
/opt/anaconda/bin/conda create --name cairo-env --file /application/dependencies/R/cairo-env.list
