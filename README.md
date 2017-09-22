## Developer Cloud Sandbox hypeapps-forecast application  

H-TEP Hydrological Modelling Thematic Application  - 10 day forecast with the Niger-HYPE model

## Quick link
 
* [Getting Started](#getting-started)
* [Installation](#installation)
* [Submitting the workflow](#submit)
* [Community and Documentation](#community)
* [Authors](#authors)
* [Questions, bugs, and suggestions](#questions)
* [License](#license)
* [Funding](#funding)

### <a name="getting-started"></a>Getting Started 

To run this application you will need a Developer Cloud Sandbox that can be requested at support (at) terradue.com

A Developer Cloud Sandbox provides Earth Sciences data access services, and helper tools for a user to implement, test and validate a scalable data processing application. It offers a dedicated virtual machine and a Cloud Computing environment.
The virtual machine runs in two different lifecycle modes: Sandbox mode and Cluster mode. 
Used in Sandbox mode (single virtual machine), it supports cluster simulation and user assistance functions in building the distributed application.
Used in Cluster mode (a set of master and slave nodes), it supports the deployment and execution of the application with the power of distributed computing for data processing over large datasets (leveraging the Hadoop Streaming MapReduce technology). 

### <a name="installation"></a>Installation

#### Pre-requisites

Requirements of the application in terms of software packages:

* miniconda
* gfortran compiler
* R and the rciop package, and the following R packages (conda and terradue compilations)

r-essentials
r-sp
r-rgdal
r-rgeos
r-raster
r-data.table
r-lmomco

To install these packages, run the simple steps below on the Developer Cloud Sandbox shell:

To install anaconda:

```bash
sudo yum install -y miniconda
```

Add anaconda binaries to PATH:

```bash
export PATH=/opt/anaconda/bin/:$PATH
```

To install R and rciop:

```bash
sudo conda install -c r  r-essentials
sudo conda install -c r -c terradue r-rciop
```

To install an R package from default channel or from Terradue channel:

```bash
sudo conda install -c r <package>
sudo conda install -c r -c terradue <package>
```

* Separate R installation with Cairo package to enable plotting (Cairo provides a graphic device, which is otherwise missing on the sandbox... no X11 there. Unfortunately, Cairo and rgdal have conflicting dependancies to jpeglib version 8 and 9, so the current workaround is to run the plotting in a separate R environment with only Cairo (and sp and data.table) installed).


#### HYPE model application

To run the application a HYPE model application is needed. The Niger-HYPE model was used in the development of the application.
However, the Niger-HYPE model is not open access, and the use of Niger-HYPE must be discussed with SMHI.



The model input files should be stored in a model subfolder, with a folder name identical to the R variable model.name:

hypeapps-forecast/src/main/app-resources/model/[model.name]



##### Example:

hypeapps-forecast/src/main/app-resources/model/my-hype

The R variable model.name is set in the util/hypeapps-enviroment.R script, under "Model settings".



##### HYPE model input file requirements

Appart from the standard HYPE input files located in the model folder defined above, the following input files are required:


1. INFO.TXT TEMPLATES: there should be two templates for info.txt, one for hindcast and one for forecast:

   hypeapps-forecast/src/main/app-resources/model/my-hype/info-hindcast.txt
   hypeapps-forecast/src/main/app-resources/model/my-hype/info-forecast.txt
   

2. SHAPEFILE: A shapefile with the model sub-basin polygons should also be included in the model folder:

   hypeapps-forecast/src/main/app-resources/model/my-hype/shapefiles/my-hype.shp


For map plots, it is furthermore necessary to prepare a Rdata file with the shapefile imported as a SpatialPolygonsDataFrame.

   hypeapps-forecast/src/main/app-resources/model/my-hype/shapefiles/my-hype.shp.Rdata
   
   (there is a small utility r script that can be used, hypeapps-forecast/src/main/app-resources/util/R/hypeapps-make-model.shp.rdata.R)



3. HYPE2CSV: A text file with information needed to transform HYPE model outputs to the TEP time-series csv format
             should be located in a sub-folder called hype2csv:

   hypeapps-forecast/src/main/app-resources/model/my-hype/hype2csv/my-hype2csv.txt


The textfile [model.name]2csv.txt should have the following structure with one entry for each sub-basin in the model:

SUBID	CENTERX	CENTERY	POURX	POURY	LON	LAT	ELEV
12388	-0.16042	17.525	-0.125	17.5125	-0.145408	17.526332	307
..


4. RETURN PERIOD LEVELS: For generation of river discharge forecast warning levels, a default file with corresponding return period levels should be provided as follows:

  hypeapps-forecast/src/main/app-resources/model/my-hype/returnperiods/my-hype-returnlevels-cout.txt
  
  
The textfile should have 4 columns, one for the subbasin identifier (SUBID), and three for selected flow levels, from which the warning levels 1-3 should be identified:

SUBID RP2  RP5  RP30
12388 0.02 0.13 0.86
..


The number in column name (RP2, RP5, etc) will be interpreted as the return period in years associated with the return level (aka warning level).



#### Meteorological forcing data

Data sources for meteorological data is selected by the R variable forcing.data.source in util/hypeapps-enviroment.R.

Currently, two data sources are available: 

forcing.data.source = "local"   will use the local file  hypeapps-forecast/src/main/app-resources/gfd/[model.name]/hindcast/archive.zip

and

forcing.data.source = "hydro-smhi" will request archive.zip from the TEP data catalogue.

However, currently there is only data for the Niger-HYPE model at the TEP data catalogue.



#### HYPE model executable

The HYPE model source code is included in the package, but the executable must be built.

After installing the package (see further below), build the HYPE executables with the following commands:

```bash
cd hypeapps-forecast/src/main/app-resources/util/fortran
./build-hype.4.12.0.sh
./build-HypeDataAssimilation.sh
```

Two executable files should now have been built and copied to the following location:

hypeapps-forecast/src/main/app-resources/util/bin/hype-4.12.0
hypeapps-forecast/src/main/app-resources/util/bin/hype_assimilation

Please also note that after installation, the file permission must be manually set on these files:

```bash
cd hypeapps-forecast
mvn clean
mvn install
chmod +x /application/util/bin/*
```

##### Using the releases

Log on the Developer Cloud Sandbox.

Install the package by running this command in a shell:

```bash
sudo yum -y install hypeapps-forecast
```

> At this stage there are no releases yet

#### Using the development version

Install the pre-requisites as instructed above.

Log on the Developer Cloud Sandbox and run these commands in a shell:

```bash
git clone https://github.com/hydrology-tep/hypeapps-forecast
cd hypeapps-forecast
mvn install
```

Please note the information above regarding the HYPE model applicaton and how to build and install the HYPE model executable files.

### <a name="submit"></a>Submitting the workflow

To submit the application with its default parameters, run the command below in the Developer Cloud Sandbox shell:

```bash
ciop-run
```
Or invoke the Web Processing Service via the Sandbox dashboard.

### <a name="community"></a>Community and Documentation

To learn more and find information go to 

* [Developer Cloud Sandbox](http://docs.terradue.com/developer-sandbox/)  

### <a name="authors"></a>Authors (alphabetically)

* David Gustafsson  ( david [dot] gustafsson [at] smhi [dot] se )
* Jafet Andersson   ( jafet [dot] andersson [at] smhi [dot] se )
* Frida Gyllensvärd ( frida [dot] gyllensvard [at] smhi [dot] se )

### <a name="questions"></a>Questions, bugs, and suggestions

Please file any bugs or questions as [issues](<app-url>) or send in a pull request if you corrected any.

### <a name="license"></a>License

Licensed under the GNU Lesser General Public License, Version 3.0: https://www.gnu.org/licenses/lgpl-3.0.en.html

### <a name="funding"></a>Funding

Put here any information about the funding.
