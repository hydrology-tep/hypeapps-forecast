# build hype v4.12.0 and copy to ../bin-file directory
cd HYPE-5.x.0
make -fmakefile_assimilation
cp ./hype_assimilation ../../bin/hype_assimilation-5.x.0.exe
make -fmakefile_assimilation clean
rm hype_assimilation
cd ..
