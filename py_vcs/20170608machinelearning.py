##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

get_ipython().run_line_magic('pylab', 'inline')
from sklearn import datasets


#%%

iris = datasets.load_iris()
digits = datasets.load_digits()


#%%

digits.target


#%%

plt.imshow(digits.images[0],cmap=plt.cm.gray_r)


#%%

from sklearn import svm
clf = svm.SVC(gamma=0.001, C=100.)


#%%

clf.fit(digits.data[:-1], digits.target[:-1])


#%%

clf.predict(digits.data[-1:])

