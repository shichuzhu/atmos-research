##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

from ipyparallel import Client
import os, socket

rc = Client()

rc.ids

rc[:].apply_sync(os.getpid)

rc[:].apply_sync(socket.gethostname)


#%%

import pickle, os


#%%

with open(os.environ['HOME']+'/Downloads/hiwc.p','rb') as f:
    hehe = pickle.load(f)


#%%

hehe.to_hdf(os.environ['HOME']+'/Downloads/hiwc.h5','niubi')

