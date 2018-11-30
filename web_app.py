import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import sys
import os.path
import flask
import cowell as cw
import numpy as np
import math as mp
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
from plotly import tools

pd.set_option("display.precision", 10)

app = dash.Dash()
app.title = 'Auth Propagation Model'
app.config['suppress_callback_exceptions']=True

propagatorObj = cw.propagator()

app.layout = html.Div([

    html.H1('AUTH Propagation Model'),

    html.Div(children=['Satellite Orbit Determination using the propagation model AUTH'], style={'margin': '20px', 'font-size':'16pt'}),

    html.Div(children=[

        html.Div(children=[
            html.H3('Initial Parameters'),
            html.Label('State Vector (m - m/s)', className="main"),
            # dcc.Textarea(value='+6465819, -1429000, +14827, -206.5, -865.7, +7707.7', id='state-vector', style={"margin-bottom": "5px", "min-height": "150px"}),
            dcc.Textarea(value='4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03', id='state-vector', style={"margin-bottom": "5px", "min-height": "150px"}),

            html.Label('Initial Time (sec)', className="main"),
            dcc.Input(value='0', type='text', id='init-time', style={"margin-bottom": "5px"}),

            html.Label('Final Time (sec)', className="main"),
            dcc.Input(value='26400', type='text', id='final-time', style={"margin-bottom": "5px"}),

            html.Label('Step Size (sec)', className="main"),
            dcc.Input(value='100', type='text', id='step-size', style={"margin-bottom": "5px"}),

            html.Label('Satellite Mass (kg)', className="main"),
            dcc.Input(value='720', type='text', id='sat-mass', style={"margin-bottom": "5px"}),
        ], style={'margin': '20px', }),

        html.Div(children=[
            html.H3('Force Model'),

            html.Label("Select Earth's Gravity Model", style={'margin-top': '10px'}, className="main"),
            dcc.Dropdown(
                id="earth-model",
                options=[
                    {'label': 'Keplerian', 'value': 1},
                    {'label': 'Geodynamic Model', 'value': 2},
                ],
                value=1
            ),
            html.Br(),

            html.Div(id="geodynamic-parameters", style={'display': 'none'}, children=[
                dcc.Dropdown(
                    id="geo-model-chose",
                    options=[
                        {'label': 'EGM08', 'value': 1},
                        {'label': 'GOCO03s', 'value': 2},
                        {'label': 'RWI TOPO 2012', 'value': 3},
                    ],
                    value=1
                ),

                html.Label('Order (n)', className="second"),
                dcc.Input(value='20', type='text', id='order-geo', style={"margin-bottom": "5px"}),
                html.Label('Degree (m)', className="second"),
                dcc.Input(value='20', type='text', id='degree-geo', style={"margin-bottom": "5px"}),

                html.Br(),

                html.Button(id="calculate-geo-button", children="Import Geodynamical Model", n_clicks=0, style={"margin-top": "5px"}),

                html.Div(id="output-geo-coeff"),

            ]),

            html.Label("Select other forces you want to include to the Model", className="main"),
            dcc.Checklist(
                id="other-forces",
                options=[
                    {'label': "Sun's and Moon's Gravitational Force", 'value': 1},
                    {'label': 'Solar Radiation', 'value': 2},
                    {'label': 'Atmospheric Drag', 'value': 3},
                ],
                values=[],
                style={"text-align": "left"}),
            html.Div(id="radiation-parameters", style={'display': 'none'}, children=[
                html.Label('Solar Radiation Reflectivity - Absorption (e)', className="main"),
                dcc.Input(value='0.8', type='text', id='solar-epsilon', style={"margin-bottom": "5px"}),
            ]),
            html.Div(id="drag-parameters", style={'display': 'none'}, children=[
                html.Label('Atmospheric Drag Coefficient (Cd)', className="main"),
                dcc.Input(value='2.1', type='text', id='drag-coeff', style={"margin-bottom": "5px"}),
            ]),
            html.Div(id="drag-model", style={'display': 'none'}, children=[
                html.Label('Atmospheric Drag Model', className="main"),
                dcc.Dropdown(
                    id="drag-model-choose",
                    options=[
                        {'label': 'Harris-Priester', 'value': 1},
                        {'label': 'Doorbnos Thermospheric Density Model (Only for GOCE)', 'value': 2},
                    ],
                    value=1
                ),
                html.Div(id="goce", style={"display": "none"}, children=[
                    html.Label('Date (yyyy-mm-dd)', className="main"),
                    dcc.Input(value='2012-11-20', type='text', id='goce-date', style={"margin-bottom": "5px", "text-align": "center"}),
                    html.Label('Time (hh:mm-ss)', className="main"),
                    dcc.Input(value='17:40:00', type='text', id='goce-time', style={"margin-bottom": "5px", "text-align": "center"}),
                ]),
            ])

        ], style={'margin': '20px'}),
    ], style={'display': 'inline-flex'}),

    html.Div(children=[

        html.Label('Please first import the Geodynamic Model before clicking Determine Satellite Orbit', id="import-error", className="info-error", style={"display": "none"}),
        html.Button(children="Determine Satellite Orbit", id="calculate-orbit", n_clicks=0)

    ], style={'margin': '40px'}),

    html.Div(id='output-state', style={'margin': '40px'})


], style={'columnCount': 1, "text-align": "center"})


@app.callback(Output('output-state', 'children'),
                [Input('calculate-orbit', 'n_clicks')],
                [State('state-vector', 'value'),
                 State('init-time', 'value'),
                 State('final-time', 'value'),
                 State('step-size', 'value'),
                 State('earth-model', 'value'),
                 State('other-forces', 'values'),
                 State('solar-epsilon', 'value'),
                 State('drag-coeff', 'value'),
                State('sat-mass', 'value'),
                State('goce-date', 'value'),
                State('goce-time', 'value'),
                State('drag-model-choose', 'value'),
                ])
def update_output(n_clicks, input1, input2, input3, input4, input5, input6, input7, input8, input9, input10, input11, input12):

    if n_clicks == 0:
        pass
    else:
        y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        index = 0
        parseInput1 = input1.split(",")
        for i in parseInput1:
            y[index] = float(i)
            index += 1

        if input2 == "" or input3 == "" or input4 == "":
            pass
        else:
            t0 = float(input2)
            tf = float(input3)
            step = float(input4)

            kepOrGeo = int(input5)
            other = list(input6)

            sunAndMoon = False
            solar = False
            drag = False

            if 1 in other:
                sunAndMoon = True

            if 2 in other:
                solar = True

            if 3 in other:
                drag = True

            loopTf = step
            loopIndex = int((tf - t0) / step)

            y = propagatorObj.rk4(y, t0, tf, step, kepOrGeo, solar, float(input7), sunAndMoon, drag, float(input8), input12, propagatorObj.C, propagatorObj.S, float(input9), input10, input11)

            altitude = np.sqrt(np.asarray(propagatorObj.keepStateVectors)[:, 0]**2 + np.asarray(propagatorObj.keepStateVectors)[:, 1]**2 + np.asarray(propagatorObj.keepStateVectors)[:, 2]**2)
            beforeDataFrame = np.zeros((loopIndex, 7))
            beforeDataFrame[:, 1:7] = propagatorObj.keepStateVectors[:]
            beforeDataFrame[:, 0] = propagatorObj.epochs[:]

            returnTable = pd.DataFrame(data=beforeDataFrame[:, :],
                                       index=range(0, loopIndex),
                                       columns=["Time(sec)", "rx(m)", "ry(m)", "rz(m)", "vx(m/s)", "vy(m/s)", "vz(m/s)"])

            propagatorObj.accelerations_graph(solar, sunAndMoon, drag)
            beforeAccelTable = np.zeros((len(propagatorObj.tickLabels), 4))
            for i in range(0, len(propagatorObj.tickLabels)):
                beforeAccelTable[i, 0] = propagatorObj.keepAbsoluteAccelerations[i][0]
                beforeAccelTable[i, 1] = propagatorObj.keepAbsoluteAccelerations[i][1]
                beforeAccelTable[i, 2] = propagatorObj.keepAbsoluteAccelerations[i][2]
                beforeAccelTable[i, 3] = propagatorObj.keepAbsoluteAccelerations[i][3]

            returnAcceTable = pd.DataFrame(
                data=beforeAccelTable,
                index=propagatorObj.tickLabels,
                columns=["Max(m^2/s)", "Min(m^2/s)", "Average(m^2/s)", "STD(m^2/s)"]
            )

            returnAcceTable.insert(0, "Source", returnAcceTable.index)

            if len(propagatorObj.tickLabels) == 1:
                depends = "accel-sizes-hide"
                depends1 = "accel-sizes-hide1"
            else:
                depends = "accel-sizes"
                depends1 = "accel-sizes1"

            # Keplerian Elements Graph
            returnStaticKepTable, kepElements = propagatorObj.keplerian_elements_graph()

            trace1 = go.Scatter(
                x=propagatorObj.epochs,
                y=kepElements[:, 0],
                opacity=0.7,
            )
            trace2 = go.Scatter(
                x=propagatorObj.epochs,
                y=kepElements[:, 1],
                opacity=0.7,
            )
            trace3 = go.Scatter(
                x=propagatorObj.epochs,
                y=kepElements[:, 2],
                opacity=0.7,
            )
            trace4 = go.Scatter(
                x=propagatorObj.epochs,
                y=kepElements[:, 3],
                opacity=0.7,
            )
            trace5 = go.Scatter(
                x=propagatorObj.epochs,
                y=kepElements[:, 4],
                opacity=0.7,
            )
            trace6 = go.Scatter(
                x=propagatorObj.epochs,
                y=kepElements[:, 5],
                opacity=0.7,
            )
            fig = tools.make_subplots(rows=3, cols=2, specs=[[{}, {}], [{}, {}], [{}, {}]],
                                      subplot_titles=('Semi Major Axis(km)', 'Eccentricity(float)',
                                                    'Inclination(Degrees)', 'Argument of Perigee(degrees)',
                                                    "Longitude of the ascending node(degrees)", "True anomaly(degrees)"),
                                      shared_xaxes=True, vertical_spacing=0.1)
            fig.append_trace(trace1, 1, 1)
            fig.append_trace(trace2, 1, 2)
            fig.append_trace(trace3, 2, 1)
            fig.append_trace(trace4, 2, 2)
            fig.append_trace(trace5, 3, 1)
            fig.append_trace(trace6, 3, 2)

            fig['layout'].update(height=1000, title='Keplerian Elements', showlegend=False)

            # Log Accel Graph
            allTraceLogs = []
            for i in range(0, len(propagatorObj.keepWantedAcceleration)):
                tracelog = go.Scatter(
                    x=propagatorObj.epochs,
                    y=propagatorObj.keepWantedAcceleration[i],
                    name=propagatorObj.keepWantedTickLabels[i],
                )
                allTraceLogs.append(tracelog)

            layoutlog = go.Layout(
                yaxis=dict(
                    type='log',
                    exponentformat='e',
                )
            )

            figlog = go.Figure(layout=layoutlog)
            figlog['layout'].update(title='Logarithmic Accelerations Graph (m^2/sec)')
            figlog.add_traces(allTraceLogs)

            return html.Div(children=
                        [html.Table(
                            # Header
                            [html.Tr([html.Th(col) for col in returnTable.columns])] +

                            # Body
                            [html.Tr([
                                html.Td(returnTable.iloc[i][col]) for col in returnTable.columns
                            ]) for i in range(min(len(returnTable), loopIndex))]
                        )], style={"padding": "20px 17%", "width": "100%"}),\
                    html.Div(id="info-download-state"), \
                    html.Button(children="Download State Vectors", id="download-state-vectors", n_clicks=0), \
                    html.Br(),\
                    dcc.Graph(
                        id='alitutde-graph',
                        figure={
                            'data': [
                                go.Scatter(
                                    x=propagatorObj.epochs,
                                    y=altitude/1000,
                                    opacity=0.7,
                                )
                            ], 'layout': {'title': 'Altitude (km)', 'exponentformat': "'e'"}
                        }
                    ), \
                    \
                    \
                   html.Div(children=
                       [html.Table(
                        # Header
                        [html.Tr([html.Th(col) for col in returnStaticKepTable.columns])] +

                        # Body
                        [html.Tr([
                            html.Td(returnStaticKepTable.iloc[i][col]) for col in returnStaticKepTable.columns
                        ]) for i in range(min(len(returnStaticKepTable), loopIndex))]
                        , id="kep-table")], style={"padding": "20px 17%", "margin": "25px 0px", "width": "100%"}), \
                   html.Div(id="info-download-kep"), \
                   html.Button(children="Download Keplerian Elements", id="download-kep-elements", n_clicks=0), \
                   html.Br(), \
                   dcc.Graph(
                        id="kep-graph",
                        figure=fig
                    ), \
                    \
                    \
                    html.Div(children=
                        [html.Table(
                            # Header
                            [html.Tr([html.Th(col) for col in returnAcceTable.columns])] +

                            # Body
                            [html.Tr([
                                html.Td(returnAcceTable.iloc[i][col]) for col in returnAcceTable.columns
                            ]) for i in range(min(len(returnAcceTable), loopIndex))]
                            , id="accel-table")], style={"padding": "20px 25%", "margin": "25px 0px"}), \
                   html.Div(id="info-download-accel"), \
                   html.Button(children="Download Acceleration Vectors", id="download-accel-vectors", n_clicks=0), \
                    html.Br(), \
                    dcc.Graph(
                        id=depends,
                        figure={
                            'data': [
                                {'x': propagatorObj.tickLabels,
                                 'y': propagatorObj.accelerationsGraph,
                                 'type': 'bar'}
                            ], 'layout': {'title': 'Accelerations Order of Magnitude'}
                        }
                    ), \
                   html.Br(), \
                   dcc.Graph(
                       id=depends1,
                       figure=figlog,
                   ),




@app.callback(Output('output-geo-coeff', 'children'),
                [Input("calculate-geo-button", "n_clicks")],
                [State('geo-model-chose', 'value'),
                State('order-geo', 'value'),
                 State('degree-geo', 'value'),
                ])
def choose_geo_model(n_clicks, input1, input2, input3):

    if n_clicks != 0:
        input2 = int(input2)
        input3 = int(input3)
        if input1 == 1:
            data = propagatorObj.read_geopotential_coeffs("restruct_EGM2008.gfc", False)
            propagatorObj.createNumpyArrayForCoeffs(data, input2, input3)
        elif input1 == 2:
            data = propagatorObj.read_geopotential_coeffs("restruct_GOCO03s.gfc", False)
            propagatorObj.createNumpyArrayForCoeffs(data, input2, input3)
        elif input1 == 3:
            data = propagatorObj.read_geopotential_coeffs("restruct_RWI_TOPO_2012_plusGRS80.gfc", False)
            propagatorObj.createNumpyArrayForCoeffs(data, input2, input3)

        n_clicks = 0


@app.callback(Output('info-download-state', 'children'),
                [Input("download-state-vectors", "n_clicks")])
def download_state_vectors(n_clicks):

    if n_clicks == 1:
        propagatorObj.download_state_vectors()
        return html.Label('The State Vectors csv has been created', className="info-message")
    elif n_clicks == 0:
        pass
    else:
        return html.Label('Error! The State Vectors csv has already been created', className="info-error")


@app.callback(Output('info-download-kep', 'children'),
                [Input("download-kep-elements", "n_clicks")])
def download_kep_vectors(n_clicks):

    if n_clicks == 1:
        propagatorObj.download_kep_vectors()
        return html.Label('The Keplerian Elements csv has been created', className="info-message")
    elif n_clicks == 0:
        pass
    else:
        return html.Label('Error! The Keplerian Elements csv has already been created', className="info-error")


@app.callback(Output('info-download-accel', 'children'),
                [Input("download-accel-vectors", "n_clicks")])
def download_accel_vectors(n_clicks):

    if n_clicks == 1:
        propagatorObj.download_accel_vectors()
        return html.Label('The Acceleration Vectors csv has been created', className="info-message")
    elif n_clicks == 0:
        pass
    else:
        return html.Label('Error! The Acceleration Vectors csv has already been created', className="info-error")


if __name__ == '__main__':
    app.run_server(debug=True)
