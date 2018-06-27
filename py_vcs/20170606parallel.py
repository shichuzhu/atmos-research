##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

import os
import ipyparallel as ipp

rc = ipp.Client()
ar = rc[:].apply_async(os.getpid)
pid_map = ar.get_dict()


#%%

def dull_printing():
    import os
    import time
    pid = os.getpid()
    for _ in range(10):
        print(str(pid) + " reporting ...")
        time.sleep(1)

import os
import ipyparallel as ipp

rc = ipp.Client()
ar = rc[:].apply_async(dull_printing)
pid_map = ar.get_dict()
print(pid_map)


#%%

ar.stdout


#%%

ar.get_dict()


#%%

dull_printing()

