mkdir build
cd build
cl /c ..\src\*.c /I..\include
lib *.obj /out:tope.lib
cl ../exe/main.c tope.lib /Fepoly.exe /I..\include
del *.obj
