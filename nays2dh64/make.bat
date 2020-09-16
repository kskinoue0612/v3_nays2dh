ifort /MD /c /Qopenmp /O3 ../src/iricmi.f90
ifort /MD /c /Qopenmp /O3 ../src/Nays2DH.f90
ifort *.obj ../libs/cgnsdll.lib ../libs/iriclib.lib ../libs/iricmi_f.lib -o nays2dh.exe
del *.obj 
del *.mod 

