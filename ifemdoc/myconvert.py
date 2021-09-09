import hglib
import os
import nbconvert

client = hglib.open('.')
l = client.status(change='tip')
files = []
for t, f in l:
    if (f.find(b'checkpoints') < 0) and f.find(b'.ipynb') >= 0:
        files.append(f.decode("utf-8"))

print("Update the following files:\n", files)
for f in files:
    os.system('jupyter nbconvert --template ifemdoc/basic-linkcss.tpl ' +f)
    os.system('jupyter nbconvert --to markdown ' +f)
