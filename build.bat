mkdir bin

rem build nays2dh.exe

ifort src\iric.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort src\Nays2DH.f90 /Qopenmp /nostandard-realloc-lhs /MD /c
ifort *.obj libs\iriclib.lib -o nays2dh.exe

rem copy nays2dh.exe to bin folder

copy nays2dh.exe bin\nays2dh.exe

rem copy files to bundle to installer
rem ---------------------------------

copy nays2dh64\definition.xml bin\definition.xml
copy nays2dh64\translation_ja_JP.ts bin\translation_ja_JP.ts
