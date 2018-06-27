##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

import numpy as np

from bokeh.layouts import layout
from bokeh.models import CustomJS, Slider, ColumnDataSource, WidgetBox
from bokeh.plotting import figure, output_file, show
from bokeh.io import output_notebook

# output_file('dashboard.html')
output_notebook()


#%%

tools = 'pan'


def bollinger():
    # Define Bollinger Bands.
    upperband = np.random.random_integers(100, 150, size=100)
    lowerband = upperband - 100
    x_data = np.arange(1, 101)

    # Bollinger shading glyph:
    band_x = np.append(x_data, x_data[::-1])
    band_y = np.append(lowerband, upperband[::-1])

    p = figure(x_axis_type='datetime', tools=tools)
    p.patch(band_x, band_y, color='#7570B3', fill_alpha=0.2)

    p.title.text = 'Bollinger Bands'
    p.title_location = 'left'
    p.title.align = 'left'
    p.plot_height = 600
    p.plot_width = 800
    p.grid.grid_line_alpha = 0.4
    return [p]


def slider():
    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    source = ColumnDataSource(data=dict(x=x, y=y))

    plot = figure(
        y_range=(-10, 10), tools='', toolbar_location=None,
        title="Sliders example")
    plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)

    callback = CustomJS(args=dict(source=source), code="""
        var data = source.data;
        var A = amp.value;
        var k = freq.value;
        var phi = phase.value;
        var B = offset.value;
        x = data['x']
        y = data['y']
        for (i = 0; i < x.length; i++) {
            y[i] = B + A*Math.sin(k*x[i]+phi);
        }
        source.change.emit();
    """)

    amp_slider = Slider(start=0.1, end=10, value=1, step=.1, title="Amplitude", callback=callback, callback_policy='mouseup')
    callback.args["amp"] = amp_slider

    freq_slider = Slider(start=0.1, end=10, value=1, step=.1, title="Frequency", callback=callback)
    callback.args["freq"] = freq_slider

    phase_slider = Slider(start=0, end=6.4, value=0, step=.1, title="Phase", callback=callback)
    callback.args["phase"] = phase_slider

    offset_slider = Slider(start=-5, end=5, value=0, step=.1, title="Offset", callback=callback)
    callback.args["offset"] = offset_slider

    widgets = WidgetBox(amp_slider, freq_slider, phase_slider, offset_slider)
    return [widgets, plot]


def linked_panning():
    N = 100
    x = np.linspace(0, 4 * np.pi, N)
    y1 = np.sin(x)
    y2 = np.cos(x)
    y3 = np.sin(x) + np.cos(x)

    s1 = figure(tools=tools)
    s1.circle(x, y1, color="navy", size=8, alpha=0.5)
    s2 = figure(tools=tools, x_range=s1.x_range, y_range=s1.y_range)
    s2.circle(x, y2, color="firebrick", size=8, alpha=0.5)
    s3 = figure(tools='pan, box_select', x_range=s1.x_range)
    s3.circle(x, y3, color="olive", size=8, alpha=0.5)
    return [s1, s2, s3]

l = layout([
    bollinger(),
    slider(),
    linked_panning(),
], sizing_mod='fixed')
show(l)


#%%

# my test
from bokeh.layouts import column
from bokeh.models import CustomJS, ColumnDataSource, HoverTool
from bokeh.plotting import Figure, output_notebook, show

output_notebook()

x0 = [x*0.005 for x in range(0, 200)]
y0 = x0
sB = ColumnDataSource(data=dict(x=x0, y=y0))
sA = ColumnDataSource(data = dict(x=[0],y=[0]))

callback = CustomJS(args=dict(sA=sA,sB=sB), code="""
    var geometry = cb_data['geometry'];
    var a = geometry.x; // current mouse x position in plot coordinates
    var b = geometry.y; // current mouse y position in plot coordinates
    var as = sA.get('data')['x'];
    var bs = sA.get('data')['y'];
    as[0] = a;
    bs[0] = b;
    sA.trigger('change');

    var data = sB.data;
    x = data['x']
    y = data['y']
    for (i = 0; i < x.length; i++) {
        y[i] = Math.pow(x[i], a+b)
    }
    sB.trigger('change');
""")

hover_tool = HoverTool(callback=callback)
A = Figure(x_range=(0,1), y_range=(0,1), tools= [hover_tool,
                        "crosshair,box_zoom,wheel_zoom,pan,reset"])
A.circle(x='x',y='y',source=sA)
B = Figure(x_range=(0,1), y_range=(0,1), tools= ["crosshair,box_zoom,wheel_zoom,pan,reset"])
B.line(x='x',y='y',source=sB)
layout = column(A,B)
show(layout)


#%%

import xarray as xr
import numpy as np
campaign = 'darwin'
sdsrc = xr.open_dataset('tmp/'+campaign+'/'+campaign+'_sync_gz.nc',group='/lamp')
bin_div = sdsrc.bin_div
bin_diff = np.diff(bin_div)
bin_mid = (bin_div[1:] + bin_div[:-1])/2

bin_mid=list(bin_mid)


#%%

# my test
from bokeh.layouts import column, row
from bokeh.models import CustomJS, ColumnDataSource, HoverTool
from bokeh.plotting import Figure, output_notebook, show

output_notebook()

sB = ColumnDataSource(data=dict(x=bin_mid, y=bin_mid))
sA = ColumnDataSource(data = dict(x=[0],y=[0]))

def callback(sA=sA,sB=sB,window=None):
    geometry = cb_data['geometry'];
    a = geometry.x;
    b = geometry.y;
    as0 = sA.get('data')['x'];
    bs0 = sA.get('data')['y'];
    as0[0] = a;
    bs0[0] = b;
    sA.trigger('change');

    data = sB.data;
    x = data['x']
    y = data['y']
    tmp = map( lambda d: d**a * 2.7182818**(-b*d) , x )
    tot = sum(tmp)
    for i in range(len(x)):
        y[i] = tmp[i]/tot
    sB.trigger('change');

hover_tool = HoverTool(callback=CustomJS.from_py_func(callback))
A = Figure(x_range=(-5,2), y_range=(0,0.01), tools= [hover_tool,
                        "crosshair,box_zoom,wheel_zoom,pan,reset"])
A.circle(x='x',y='y',source=sA)
B = Figure(x_range=(10,13000), y_range=(1e-19,1), x_axis_type="log", y_axis_type="log", 
           tools= ["crosshair,box_zoom,wheel_zoom,pan,reset"])
B.line(x='x',y='y',source=sB)
layout = column(A,B)
show(layout)

