ifort 	..\libs\cgnsdll.lib^
	..\libs\iriclib.lib ^
	..\src\Nays2DH.f90 /Qopenmp /nostandard-realloc-lhs /MD -o nays2dh.exe
	

del *.obj 
del *.mod 

