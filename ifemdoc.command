#!/bin/bash
rsync -av --exclude ".git/" ~/Dropbox/Math/Programming/GitHub/ifem/ ~/Dropbox/Math/Programming/ifem/
cd ~/Dropbox/Math/Programming/ifem/
./hgifemdoc.command
rsync -av --exclude ".hg/" ~/Dropbox/Math/Programming/ifem/ifemdoc/ ~/Dropbox/Math/Programming/GitHub/ifem/ifemdoc/
cd ~/Dropbox/Math/Programming/GitHub/ifem/
git add ifemdoc/;
git commit;
git pull;
git push;
rsync -av —delete -e ssh --progress ~/Dropbox/Math/Programming/ifem/ifemdoc/ /Users/Shared/Dropbox/Sites/public_html/ifemdoc/
rsync -av --delete -e ssh --progress /Users/Shared/Dropbox/Sites/public_html/ chenlong@home.ps.uci.edu:public_html/
