#!/bin/bash
git add ifemdoc/;
git commit -m "jupyter note book files";
python ifemdoc/myconvert.py;
git add ifemdoc/;
git commit;
git push;
rsync -av --exclude ".git/" ~/Dropbox/Math/Programming/GitHub/ifem/ ~/Dropbox/Math/Programming/ifem/
