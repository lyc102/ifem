#!/bin/bash
rsync -av --exclude ".git/" ~/Dropbox/Math/Programming/GitHub/ifem/ ~/Dropbox/Math/Programming/ifem/
cd ~/Dropbox/Math/Programming/ifem/
./hgifemdoc.command
rsync -av --exclude ".hg/" ~/Dropbox/Math/Programming/ifem/ifemdoc/ ~/Dropbox/Math/Programming/GitHub/ifem/ifemdoc/
cd ~/Dropbox/Math/Programming/GitHub/ifem/
git add ifemdoc/;
git rm -r ifemdoc/afem/.ipynb_checkpoints/
git rm -r ifemdoc/fem/.ipynb_checkpoints/
git rm -r ifemdoc/mesh/.ipynb_checkpoints/
git rm -r ifemdoc/project/.ipynb_checkpoints/
git rm -r ifemdoc/solver/.ipynb_checkpoints/
git commit;
git pull;
rsync -av —delete -e ssh --progress ~/Dropbox/Math/Programming/ifem/ifemdoc/ ~/Dropbox/Sites/public_html/ifemdoc/
rsync -av --delete -e ssh --progress ~/Dropbox/Sites/public_html/ chenlong@home.ps.uci.edu:public_html/
