##%%writefile jupyter_notebook_exporter_python_mycfg.tpl
{% extends 'python.tpl'%}

## remove markdown cells
{% block markdowncell -%}
{% endblock markdowncell %}

## change the appearance of execution count
## Rule based on: https://nbconvert.readthedocs.io/en/latest/customizing.html#Custom-Templates
## Format based on: https://github.com/DonJayamanne/vscodeJupyter#sample-python-usage
{% block in_prompt %}
#%%
{% endblock in_prompt %}
