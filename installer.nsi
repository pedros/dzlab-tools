outFile "dzlab-tools-0.0.16.exe"
 
InstallDir C:\dzlab-tools
!define UNINSTALLER $INSTDIR\uninstaller.exe

# execute unintaller if it already exists
# http://nsis.sourceforge.net/When_I_use_ExecWait_uninstaller.exe_it_doesn't_wait_for_the_uninstaller
Function .onInit
ExecWait '"${UNINSTALLER}" _?=$INSTDIR'
FunctionEnd
 
section
setOutPath $INSTDIR
file /r /x utilities /x .* /x .git /x *.exe /x *.nsi /x t /x tmp /x *.gff *
writeUninstaller ${UNINSTALLER}
CreateDirectory "$SMPROGRAMS\Zilberman Lab"
CreateShortCut "$SMPROGRAMS\Zilberman Lab\Help.lnk" "$INSTDIR\help\index.html"
sectionEnd
 
section "Uninstall"
SetAutoClose true
Delete "$SMPROGRAMS\Zilberman Lab\Help.lnk"
delete ${UNINSTALLER}
rmdir /r $INSTDIR\*
sectionEnd

