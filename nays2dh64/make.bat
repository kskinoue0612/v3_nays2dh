ifort ..\src\iric.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort ..\src\Nays2DH.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort *.obj ..\libs\iriclib.lib -o nays2dh.exe
