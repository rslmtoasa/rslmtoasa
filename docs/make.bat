@REM Sphinx batch file for Windows
@REM
@REM This will build the Sphinx documentation for Windows.

@SETLOCAL

set SPHINXBUILD=sphinx-build
set SOURCEDIR=source
set BUILDDIR=_build

if "%1"=="" goto help

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%
goto end

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%

:end
@ENDLOCAL
