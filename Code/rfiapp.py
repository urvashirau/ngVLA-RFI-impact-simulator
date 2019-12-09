#!/usr/bin/env python

"""rfiapp.py: Create a dash/plotly GUI for an RFI impact estimator."""

__author__      = "Urvashi R.V."
__email__ = "rurvashi@nrao.edu"

import dash
from dash.dependencies import Input, Output
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

#execfile('rfisim.py')

from rfisim import *

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

df = get_rfi_df(rfi)
#print "Initial RFI df" , df

app.layout = html.Div(
    children=[
                    html.H3(children='RFI Impact Simulator for the Next Generation VLA'),
        html.Div( [
            html.Div( 
                children=[ 
                    html.H6(children='RFI Characteristics'),
                    dash_table.DataTable(
                        id='rfitable',
                        #columns=[{"name": i, "id": i} for i in df],
                        columns = [ { 'id':'name', 'name': 'Type of RFI', 'type':'text'}, 
                                    { 'id':'arrayfrac', 'name': 'Array Visibility', 'type':'text'}, 
                                    { 'id':'timefrac', 'name': 'Time Fraction', 'type':'text'}, 
                                    { 'id':'timeres', 'name': 'Sig Len (s)', 'type':'text'}, 
                                    { 'id':'timegap', 'name': 'Sig Gap (s)', 'type':'text'}, 
                                    { 'id':'freqres', 'name': 'Chan Width (kHz)', 'type':'text'}, 
                                    { 'id':'freqrange', 'name': 'Freq Band [start, end] (GHz)', 'type':'text'}, 
                                    { 'id':'types', 'name': 'RFI Sources', 'type':'text'}],
                        data=df.to_dict('records'),
                        style_cell={'textAlign': 'right'},
                        style_data={ 'border': '1px solid black' },
                        style_header={ 'border': '1px solid black' },
                        editable=True,
                        style_table={'overflowX' : 'scroll'}, #, 'overflowY' : 'scroll' },
                    ), 
                ], style={'width':'50%', 'height':'50%','display':'inline-block', 'vertical-align':'top'}) , 
            html.Div( 
                children=[ 
                    
                    html.H6(children='RFI mitigation options'),
                    dcc.Checklist(
                        id='algorithm-picker',
                        options=[
                            {'label': 'Post-Processing Flagging', 'value': 'D', 'disabled':True},
                            {'label': 'Antenna-based Real Time Flagging', 'value': 'A'},
                            {'label': 'Baseline-based High time resolution Flagging (in-correlator)', 'value': 'B'},
                            {'label': 'RFI modeling and subtraction at high time resolution', 'value': 'C'}
                        ],
                        value=['D'],
                        #labelStyle={'display': 'inline-block'}
                    ),
                    
                    html.Div(
                        
                        children = [
                            
                            html.Div([
                                html.H6(children='RFI Decorrelation'),
                                dcc.RadioItems(
                                    id='decorr-picker',
                                    options=[
                                        {'label': 'None (ignore from calculations)', 'value': 'no_decorr'},
                                        {'label': 'RFI at 20deg from phase-center (practical estimate)', 'value': 'prac_decorr'},
                                        {'label': 'RFI at 90deg from phase-center (maximal decorrelation)', 'value': 'max_decorr'}
                                    ],
                                    value='prac_decorr',
                                    #labelStyle={'display': 'inline-block'}
                                )
                            ], style={'width':'50%', 'height':'50%','display':'inline-block', 'vertical-align':'top'} ),
                            
                            html.Div([
                                html.H6(children='Attenuation threshold'),
                                dcc.RadioItems(
                                    id='dynr-picker',
                                    options=[
                                        {'label': '20 dB', 'value': '1e-02'},
                                        {'label': '40 dB', 'value': '1e-04'},
                                        {'label': '60 dB', 'value': '1e-06'},
                                    ],
                                    value='1e-04',
                                    labelStyle={'display': 'inline-block'}
                                )
                            ],style={'width':'50%', 'height':'50%','display':'inline-block', 'vertical-align':'top'} )
                            
                        ] 
                    )
                ], style={'width':'45%', 'height':'50%','display':'inline-block', 'vertical-align':'top', 'marginLeft':'5%'}) ] ), 
        
        html.Div( [
            html.Div( 
                children=[

                    #                html.H6(children='Fraction of data lost to RFI'),
                    dcc.Graph(id='rfi_loss_frac')


                ], style={'width': '50%', 'display': 'inline-block','vertical-align':'top'}),
            html.Div( 
                children=[ 
                    dcc.Graph(id='extra_observing')
                ], style={'width':'50%', 'height':'50%','display':'inline-block', 'vertical-align':'top'}) , 
#            html.Div( 
#                children=[ 
#                    dcc.Graph(id='compute_intensity')
#                ], style={'width':'33%', 'height':'50%','display':'inline-block', 'vertical-align':'top'})  
        ] )
    ]
)

@app.callback(
    [Output('rfi_loss_frac', 'figure'),
     Output('extra_observing', 'figure')],
#     Output('compute_intensity','figure')],
    [Input('rfitable', 'data'),
     Input('rfitable', 'columns'),
     Input('algorithm-picker','value'),
     Input('decorr-picker','value'),
     Input('dynr-picker', 'value')
 ])
def update_figure(rfidata, rficolumns, algolist, decorr, dynr):

    df = pd.DataFrame(rfidata) #, columns=[c['name'] for c in rficolumns])
    rfidict = get_rfi_dict( df )

    ### Make plot of RFI loss fraction
    totalx, totals, procpars = get_data_loss_fraction(rfichar = rfidict, sol=algolist, decorr=decorr, dynr=dynr)
    extratime_min, extratime_max, fracloss_min, fracloss_max = get_observing_time_ratio(totalx, totals)
    
#    compute_intensity = get_compute_intensity( totalx, procpars )

    ##################################################################

    traces1=[]
    for atom in totals.keys():
        traces1.append( go.Scatter(x=totalx, y=totals[atom], fill='tozeroy',mode='none',name=atom))
    
    # Add band labels
    bandcols = ['green','cyan','green','cyan','green','cyan']
    bcnt=1
    for band in ngvlabands:
        traces1.append(go.Scatter(
            x=[band[0], (band[0]+band[1])/2.0, band[1]],
            y=[0.95, 0.95, 0.95],
#            x=[band[0], band[1]],
#            y=[0.95, 0.95],
            mode="lines+text",
            line_color=bandcols[bcnt-1],
            line_width=3,
            name="Band"+str(bcnt),
            text=["", "Band "+str(bcnt),""],
#            text=["    Band "+str(bcnt),""],
            textposition="bottom center",
            showlegend=False
        ))
        bcnt = bcnt+1

    # Add the oxygen band
    traces1.append(go.Scatter( x=[50.0, 70.0], y=[1.0,1.0], fill='tozeroy', mode='none', 
                               name='Atm Oxygen' , fillcolor='rgba(80,80,80,0.3)' ) ) 

    ##################################################################

    ### Make plot of extra observing time required.
    traces2 = []
    traces2.append( go.Scatter(x=totalx, y=extratime_min, name='RFI from multiple\nsources overlap\nin the data', 
                               showlegend=True,
                               opacity=0.5))
    traces2.append( go.Scatter(x=totalx, y=extratime_max, name='RFI from multiple\nsources affect different\nsubsets of the data', 
                               showlegend=True,
                               opacity=0.5))
    
    ymax = np.min( [np.max(extratime_max)+0.5, 12.0] )
    # Add band labels
    bcnt=1
    for band in ngvlabands:
        traces2.append(go.Scatter(
            x=[band[0], (band[0]+band[1])/2.0, band[1]],
#            y=[0.95+0.007*bcnt, 0.95+0.007*bcnt, 0.95+0.007*bcnt],
            y=[0.95*ymax, 0.95*ymax, 0.95*ymax],
            mode="lines+text",
            line_color=bandcols[bcnt-1],
            line_width=3,
            name="Band"+str(bcnt),
            text=["", "Band "+str(bcnt),""],
            textposition="bottom center",
            showlegend=False
        ))
        bcnt = bcnt+1
    # Add the oxygen band
    traces2.append(go.Scatter( x=[50.0, 70.0], y=[ymax, ymax], 
                               fill='tozeroy', mode='none', 
                               name='Atm Oxygen' , fillcolor='rgba(80,80,80,0.3)' , showlegend=False) ) 
    ##################################################################

    ### Make plot of compute load
#    traces3=[]
#    maxx = 0.0
#    for atom in compute_intensity.keys():
#        traces3.append( go.Scatter(x=totalx, y=compute_intensity[atom], fill='tozeroy',mode='none', name=atom))
#        if np.max(compute_intensity[atom]) > maxx:
#            maxx = np.max(compute_intensity[atom])

#    ymax = maxx
#    # Add band labels
#    bcnt=1
#    for band in ngvlabands:
#        traces3.append(go.Scatter(
#            x=[band[0], (band[0]+band[1])/2.0, band[1]],
#            y=[0.95*ymax, 0.95*ymax, 0.95*ymax],
#            mode="lines+text",
#            line_color=bandcols[bcnt-1],
#            line_width=3,
#            name="Band"+str(bcnt),
#            text=["", "Band "+str(bcnt),""],
#            textposition="bottom center",
#            showlegend=False
#        ))
#        bcnt = bcnt+1
#    # Add the oxygen band
#    traces3.append(go.Scatter( x=[50.0, 70.0], y=[ymax, ymax], 
#                               fill='tozeroy', mode='none', 
#                               name='Atm Oxygen' , fillcolor='rgba(80,80,80,0.3)' , showlegend=False) ) 
    ##################################################################



    return [{
        'data': traces1,
        'layout': go.Layout(
            #            xaxis={'title': 'Frequency (GHz)','range':[1.0,120.0]},
            xaxis={'title': 'Frequency (GHz)','type':'log', 'range':[0,pl.log(np.max(totalx))/2.3]},
            yaxis={'title': 'Fraction of data loss','range':[0,1]},
            #xaxis_type="log",
            title="Fraction of data loss   [Average :  %2.1f%% to %2.1f%% ]"%(fracloss_min*100.0, fracloss_max*100.0),
            #            width = 1000, height = 500,
            autosize = True,
            legend_orientation = 'h',
            legend=dict(x=0.05, y=1.15)
        ) 
    },
            {
                'data': traces2,
                'layout': go.Layout(
                    #            xaxis={'title': 'Frequency','range':[1.0,120.0]},
                    xaxis={'title': 'Frequency','type':'log', 'range':[0,pl.log(np.max(totalx))/2.3]},
                    yaxis={'title': 'Scale factor for observing time','range':[0.9, np.min( [np.max(extratime_max)+0.5 , 12] ) ]},
                    #xaxis_type="log",
                    title="Observing time required to reach target sensitivity (scale factor = 1 for no RFI)",
                    #            width = 1000, height = 500,
                    autosize = True,
                    legend_orientation = 'h',
                    legend=dict(x=0.05, y=1.15)
                ) 
            },
#            {
#                'data': traces3,
#                'layout': go.Layout(
#                    #            xaxis={'title': 'Frequency','range':[1.0,120.0]},
#                    xaxis={'title': 'Frequency','type':'log', 'range':[0,pl.log(np.max(totalx))/2.0]},
#                    yaxis={'title': 'Compute Intensity','type':'log','range':[0,12]},
#                    #xaxis_type="log",
#                    title="Compute Intensity = Data rate  x  Number of operations per datum",
#                    #            width = 1000, height = 500,
#                    autosize = True,
#                    legend_orientation = 'h',
#                    legend=dict(x=0.05, y=1.15)
#                ) 
#            }
        ]
    

if __name__ == '__main__':
    app.run_server(debug=True)
