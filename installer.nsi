outFile "installer.exe"
 
InstallDir C:\dzlab-tools
!define UNINSTALLER $INSTDIR\uninstaller.exe

# execute unintaller if it already exists
Function .onInit
Exec ${UNINSTALLER}
FunctionEnd
 
section
setOutPath $INSTDIR
file /r /x utilities /x .* /x .git /x *.exe /x *.nsi /x tmp *
writeUninstaller ${UNINSTALLER}
sectionEnd
 
section "Uninstall"
delete ${UNINSTALLER}
rmdir /r $INSTDIR\*
sectionEnd

