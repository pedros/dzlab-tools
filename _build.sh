
rm -rf help
pod2projdocs -o help -l . -except doc/index.pod -except '^tmp'
makensis installer.nsi
