call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019

mkdir INSTALL

rem build nays2dh.exe

ifort src\iric.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort src\Nays2DH.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort *.obj libs\iriclib.lib -o nays2dh.exe

rem copy nays2dh.exe to INSTALL folder

copy nays2dh.exe INSTALL\nays2dh.exe

rem copy files to bundle to installer
rem ---------------------------------

copy nays2dh64\definition.xml INSTALL\definition.xml
copy nays2dh64\translation_ja_JP.ts INSTALL\translation_ja_JP.ts
