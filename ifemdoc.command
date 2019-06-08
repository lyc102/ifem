#!/bin/bash
git add ifemdoc/;
git commit -m "nb_files";
python ifemdoc/myconvert.py;
git add ifemdoc/;
git commit;
git push;
