print("importing libraries")
import dash
import dash_core_components as dcc
import dash_html_components as html
# from dash_table import DataTable
import plotly.graph_objs as go
import text
import numpy as np
from OptiScan import scan, database, utils
from OptiScan.database import MoleculeDB
import subprocess as sp
from os import path
import os
import flask



def plot_graph(figures, gid, layout=None):
    if not layout:
        return dcc.Graph(id=gid, figure={"data": figures})
    else:
        return dcc.Graph(id=gid, figure={"data": figures, "layout": layout})


def molecule_plot_data(label_data, backbone_data, snrx, snry):
    xs = np.linspace(0, label_data.shape[0]/2, label_data.shape[0])
    return [go.Scatter(x=xs, y=label_data, name="Specific lab."),
            go.Scatter(x=xs, y=backbone_data, name="Backbone lab."),
            go.Scatter(x=snrx/2, y=snry, name="Nicking sites", mode="markers")]


def histogram(backbones):
    backbone_data = np.array([(np.mean(x), np.max(x)) for x in backbones])
    return [go.Histogram(x=[x[0] for x in backbone_data], name="Mean backbone intensity", nbinsx=25,
                         marker=dict(color='rgb(158,202,225)', line=dict(color='rgb(8,48,107)', width=1.5))),
            go.Histogram(x=[x[1] for x in backbone_data], name="Max backbone intensity", nbinsx=25,
                         marker=dict(color='rgb(158,202,225)', line=dict(color='rgb(8,48,107)', width=1.5)))]


def histogram_labels(mols, snr):
    nick_counts = [len(utils.get_bnx_info(x, snr)["nick_distances"]) for x in mols]
    return go.Histogram(x=nick_counts, name="Specific label content", nbinsx=25,
                         marker=dict(color='rgb(158,202,225)', line=dict(color='rgb(8,48,107)', width=1.5)))


def histogram_lens(mols):
    lens = [len(x)/2 for x in mols]
    return go.Histogram(x=lens, name="Molecule length", nbinsx=25,
                         marker=dict(color='rgb(158,202,225)', line=dict(color='rgb(8,48,107)', width=1.5)))


def get_nicking_sites(label_data, snr):
    xs = np.array(utils.get_bnx_info(label_data, snr)["nick_distances"])
    return xs, label_data[xs]


runs = "/Users/akdel/PycharmProjects/test_optiscan/OptiScan/test_data/"

box_style={"box-shadow":"1px 3px 20px -4px rgba(0,0,0,0.75)",
           "border-radius": "5px", "background-color": "#f9f7f7"}

box_style_lg={"top-margin": 25,
              "border-style": "solid",
              "border-color": "rgb(187, 187, 187)",
              "border-width": "1px",
              "border-radius": "5px",
              "background-color": "#edfdff"}

print("imported libraries")

external_stylesheets = ["https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.css"]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets, url_base_pathname="/optiscan/")

print("loaded stylesheet and starting html server")


def get_db_stats(run_info):
    avg_mol_size = 0
    number_of_mols = 0
    max_len = 0
    for run_id in run_info:
        current_run = run_info[run_id]
        for scan in current_run:
            cols = current_run[scan]["columns"]
            for col in cols:
                current_col = cols[col]
                avg_mol_size += current_col["mean_molecule_length"] * current_col["molecule_count"]
                number_of_mols += current_col["molecule_count"]
                if max_len < current_col["longest_shape"]:
                    max_len = current_col["longest_shape"]
    return number_of_mols, (avg_mol_size/number_of_mols) / 2, max_len / 2


def empty_scan_page():
    header = html.H2("Molecule Detection", style={"text-align": "center"})
    text_body = html.P('Adjust parameters below and press "Run OptiScan" to extract and record optical map molecules:',
                       style={"margin": 10})
    styles = {"text-align": "right"}
    p1 = html.Div(children=[html.Div(html.H6("Database Name: ", style=styles), className="two columns"),
                          html.Div(dcc.Input(placeholder='database', type='text',
                                                      value='', id="db-name"),
                                   className="three columns")])
    p2 = html.Div(children=[html.Div(html.H6("Chip Dimensions: ", style=styles), className="three columns"),
                         html.Div(dcc.Dropdown(options=[{'label': '12 x 95', 'value': "12, 95"},
                                                        {'label': '12 x 120', 'value': "12, 120"},
                                                        {"label": "test dimensions", "value": "12, 5"}],
                                                        placeholder="test dimensions",
                                                        value="12, 5", id="chip-dimension"),
                                  className="three columns")],)
    p3 = html.Div(children=[html.Div(html.H6("Platform: ", style=styles), className="two columns"),
                         html.Div(dcc.Dropdown(options=[{'label': 'Irys', 'value': "irys"},
                                               {'label': 'Saphyr', 'value': "saphyr"}],
                                               value="irys",
                                               placeholder="Irys", id="platform"),
                                  className="three columns")],)
    p4 = html.Div(children=[html.Div(html.H6("Run Folders: ", style=styles), className="three columns"),
                          html.Div(dcc.Input(placeholder='test_data/test_run', type='text',
                                                      value='test_data/test_run', id="folders-name"),
                                   className="three columns")])
    p5 = html.Div(children=[html.Div(html.H6("Num. of Threads: ", style=styles), className="two columns"),
                            html.Div(dcc.Input(placeholder="2", value="2",type="text", id="threads"),
                                     className="three columns")], )
    p6 = html.Div(children=[html.Div(html.H6("Organism: ", style=styles), className="three columns"),
                            html.Div(dcc.Input(placeholder='test_organism', type='text',
                                               value='test_organism', id="organism"),
                                     className="three columns")])

    container1 = html.Div(children=[p3, p4], className="container")
    container2 = html.Div(children=[p1, p2], className="container")
    container3 = html.Div(children=[p5, p6], className="container")
    container4 = html.Div([html.Button("Run OptiScan", className="three columns", id="run-optiscan",
                                      style={"background":"lightgreen"}),
                                      html.H6("Molecule detection running..", id="optiscan-running", style={'display': 'none'}, className="three columns"),
                                      html.H6("Molecule detection completed.", id="optiscan-completed", style={'display': 'none'}, className="three columns")],
                         className="container")

    return html.Div([html.Br(), header, html.Br(), container1, container2, container3, container4, html.Br()],
                    style=box_style, id="scan-page")


def empty_database_page():
    header = html.H2("Inspect & Export", style={"text-align": "center"})
    text_body = html.P('Access to molecule database to view scan and molecules statistics and visualize molecules',
                       style={"margin": 15})
    styles = {"text-align": "right", "margin-right": 5}
    sub_h1 = html.H4("1. Connect", style={"margin": 15, "text-align": "center", "padding-top": 25})
    database_names = os.listdir("database")
    p1 = html.Div([html.Div(html.H6("Database Name: ", style=styles), className="two columns"),
                   html.Div(dcc.Dropdown(options=[{"label": x, "value": "database/%s" %  x} for x in database_names if x.endswith(".db")],
                                      value='', id="db-location-access"), className="three columns"),
                   html.Button("Connect", style={"background":"lightblue"}, className="three columns", id="connect-db"),
                   html.Button(html.Img(src="https://cdn2.iconfinder.com/data/icons/dark-action-bar-2/96/refresh-512.png", width=35, height=35),
                               style={"background": "lightblue"}, className="three columns", id="refresh")],  ###################


                  className="container")
    p1_res = html.Div(None, className="container", id="db-summary-field", style=box_style_lg)
    sub_h2 = html.H4("2. Check database content", style={"margin": 15, "text-align": "center", "padding-top": 25})
    p2 = html.Div([html.Div(html.H6("Run id: ", style=styles), className="two columns"),
                   html.Div(dcc.Dropdown(options=[{"label":"", "value":""}], value="", id="run-id"), className="three columns", id="run-field"),
                   html.Div(html.H6("Scan no: ", style=styles), className="two columns"),
                   html.Div(dcc.Dropdown(options=[{"label":"", "value":""}], value="", id="scan-id"), className="three columns", id="scan-field")],
                  className="container")
    p3 = html.Div([html.Div(html.H6("Column id: ", style=styles), className="two columns"),
                   html.Div(dcc.Dropdown(options=[{"label":"", "value":""}], value="", id="column-id"), className="three columns", id="column-field")],
                  className="container")
    p4 = html.Div([html.Br(), html.Br()], id="molecule-stats", className="container", style=box_style_lg)
    p5 = html.Div([html.Div(html.H6("Molecule id: ", style=styles), className="two columns"),
                   html.Div(dcc.Dropdown(options=[{"label": "", "value": ""}], value="", id="molecule-id"),
                            className="three columns"),
                   html.Div(html.H6("SNR: ", style=styles), className="one column"),
                   html.Div(dcc.Slider(min=2.5, max=15, step=0.2, value=3, id="snr"), className="three columns"),
                   html.Div(html.H6(""), id="snr-out", className="one column")],
                  className="container")
    sub_h3 = html.H4("3. Inspect Raw Molecules ", style={"margin": 15, "text-align": "center", "padding-top": 25})
    mol_plot = html.Div(plot_graph(molecule_plot_data(np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)), "mol"),
                        id="molecule-plot")
    container1 = html.Div([p1, p1_res, html.Br(), #dcc.Markdown("------------"),
                           sub_h2, html.Br(), p2, p3, html.Br(), p4, html.Br(), #dcc.Markdown("------------"),
                           sub_h3, html.Br(), p5, mol_plot])
    sub_h3 = html.H4("4. Quick Stats. ", style={"margin": 15, "text-align": "center", "padding-top": 25})
    backbone_hist = html.Div([html.Br(),
                              html.Div([html.Button("Mol. length", className="three columns",
                                                    style={"background":"lightblue"}, id="mol-length-button"),
                                        html.Button("Label content", style={"background": "lightblue"},
                                                    className="three columns", id="label-density-button"),
                                        html.Button("Backbone Avg.", style={"background":"lightblue"},
                                                   className="three columns", id="backbone-avg-button"),
                                        html.Button("Backbone Max.", style={"background":"lightblue"},
                                                   className="three columns", id="backbone-max-button")],
                                       className="row"),
                              html.Br(), html.Div(plot_graph([], "backbone-hist"), className="row", id="histogram-field", hidden=True)],
                             id="backbone-field", className="container")

    header2 = html.H4("5. Export as BNX", style={"text-align": "center"})
    text_body = html.P('BNX conversion tool',
                       style={"margin": 15})
    styles = {"text-align": "right"}
    database_names = os.listdir("database")
    container3 = html.Div([html.Div([html.Div(html.H6("Database Location: ", style=styles), className="two columns"),
                                     html.Div(dcc.Dropdown(
                                         options=[{"label": x, "value": "database/%s" % x} for x in database_names if x.endswith(".db")],
                                         value='', id="db-location-bnx"), className="three columns"),
                                     html.Div(html.H6("Output filename:", style=styles), className="two columns"),
                                     html.Div(dcc.Input(type='text', value='', id="filename"), className="three columns")], className="row"),
                           html.Div([html.Div(html.H6("Signal to noise ratio (SNR): "), className="two columns", style=styles),
                                     html.Div(dcc.Input(id="snr-bnx-input", type="text", value="3.3"), className="three columns"),
                                     html.Div(html.H6("Minimum Length: "), className="two columns", style=styles),
                                     html.Div(dcc.Input(id="min-len-bnx-input", type="text", value="175"), className="three columns")],
                                    className="row"),
                           html.Div([html.Div(html.H6("Max. backbone intensity: "), className="two columns", style=styles),
                                    html.Div(dcc.Input(id="intensity-bnx-input", type="text", value="10000"), className="three columns")],
                                    className="row")],
                          className="container")
    container4 = html.Div(html.Button("Export", className="three columns", id="bnx-convert-button",
                                      style={"background": "lightgreen"}), className="container")
    download_link = html.Div(html.P(""), id="link-field")

    return html.Div([html.Br(), html.Br(),header, sub_h1, html.Br(), container1,
                     sub_h3, backbone_hist, html.Br(),
                     html.Br(), header2, html.Br(),
                     container3, container4, html.Br(), download_link, html.Br()], style=box_style, id="db-page")


def gitlab_link_optiscan():
    return html.A([html.Img(width=30, height=30,
                             src="https://res-1.cloudinary.com/crunchbase-production/image/upload/c_lpad,h_256,w_256,f_auto,q_auto:eco/v1436005432/xnceesad5dbk42jdwlcp.png",
                             ), html.Small("OptiScan library and cli", style={"text-align": "center"})], href="https://gitlab.com/akdel/OptiScan")

app.layout = html.Div([html.Div([
    html.Br(),
    html.H1('OptiScan', style={"text-align": "center", "font-weight": "bold", "font-family": "Georgia, serif",
                               "font-size": "75px"}),
    html.Br(),
    html.P(text.Intro, style={"margin": 10}),
    html.Br(),
    html.Br(),
    empty_scan_page(),
    html.Br(),
    html.Div(id="scan-page-result"),
    html.Br(),
    empty_database_page(),
    html.Br()], className="container"),
    gitlab_link_optiscan()])




@app.callback(dash.dependencies.Output("optiscan-running", "style"),
             [dash.dependencies.Input("run-optiscan", "n_clicks"),
              dash.dependencies.Input("optiscan-completed", "style")],
             [dash.dependencies.State("db-name", "value"),
              dash.dependencies.State("chip-dimension", "value"),
              dash.dependencies.State("platform", "value"),
              dash.dependencies.State("folders-name", "value"),
              dash.dependencies.State("threads", "value"),
              dash.dependencies.State("organism", "value")])
def optiscan_running_response(click, completed, db_name, dim, platform, runs_path, threads, organism_name):
    if not click or not db_name or not dim or not platform or not runs_path:
        return {"display": "none"}
    elif completed["display"] == "block":
        return {"display": "none"}
    else:
        return {"display": "block"}



@app.callback(dash.dependencies.Output("optiscan-completed", "style"),
             [dash.dependencies.Input("run-optiscan", "n_clicks")],
             [dash.dependencies.State("db-name", "value"),
              dash.dependencies.State("chip-dimension", "value"),
              dash.dependencies.State("platform", "value"),
              dash.dependencies.State("folders-name", "value"),
              dash.dependencies.State("threads", "value"),
              dash.dependencies.State("organism", "value")])
def run_optiscan(click, db_name, dim, platform, runs_path, threads, organism_name):
    if not db_name or not dim or not platform or not runs_path:
        return html.H1("")
    else:
        dim = dim.split(",")
        dim = f"{int(dim[0])},{int(dim[1])}"
        print(dim)
        db_name = f"database/{db_name}"
        cmd = f"python3 ./extract_molecules.py {runs_path} {dim} {db_name} {threads} {organism_name} {platform}"
        # scanner(db_name, dim, platform, runs_path)
        sp.check_call(cmd, shell=True)
        return {"display": "block"}



@app.callback(dash.dependencies.Output("run-field", "children"),
              [dash.dependencies.Input("connect-db", "n_clicks")],
              [dash.dependencies.State("db-location-access", "value")])
def update_run_dropdown(clicks, val):
    if clicks and val:
        print(val)
        mc = database.MoleculeConnector(val)
        # print(mc.db_location, mc.load_db_to_class(), mc.db_runs)
        runs = list(mc.db_runs.keys())
        opts = [{"label": runs[x], "value": runs[x]} for x in range(len(runs))]
        mc.db.close()
        print(runs)
        return dcc.Dropdown(options=opts, value="", id="run-id")
    return dcc.Dropdown(options=[{"label": "", "value": ""}], value="", id="run-id")


@app.callback(dash.dependencies.Output("scan-field", "children"),
              [dash.dependencies.Input("run-id", "value")],
              [dash.dependencies.State("db-location-access", "value")])
def update_scan_dropdown(run_id, val):
    if run_id and val:
        mc = database.MoleculeConnector(val)
        run_info = mc.db_runs[run_id]
        mc.db.close()
        opts = [{"label": str(x), "value": str(x)} for x in run_info]
        return dcc.Dropdown(options=opts, value="", id="scan-id")
    return dcc.Dropdown(options=[{"label": "", "value": ""}], value="", id="scan-id")


@app.callback(dash.dependencies.Output("column-field", "children"),
              [dash.dependencies.Input("run-id", "value"),
               dash.dependencies.Input("scan-id", "value")],
              [dash.dependencies.State("db-location-access", "value")])
def update_column_dropdown(run_id, scan_id, val):
    if run_id and scan_id and val:
        mc = database.MoleculeConnector(val)
        column_info = mc.db_runs[run_id][int(scan_id)]["columns"]
        mc.db.close()
        opts = [{"label": str(x), "value": str(x)} for x in column_info]
        return dcc.Dropdown(options=opts, value="", id="column-id")
    return dcc.Dropdown(options=[{"label": "", "value": ""}], value="", id="column-id")

@app.callback(dash.dependencies.Output("molecule-id", "options"),
              [dash.dependencies.Input("connect-db", "n_clicks")],
              [dash.dependencies.State("db-location-access", "value")])
def update_molecule_dropdown(n_clicks, val):
    if n_clicks and val:
        mc = database.MoleculeConnector(val)
        num_mols = mc._get_number_of_molecules_in_db()
        mc.db.close()
        return [{"label": str(x), "value": str(x)} for x in range(1, num_mols-1)]
    return [{"label": "", "value": ""}]


@app.callback(dash.dependencies.Output("molecule-plot", "children"),
              [dash.dependencies.Input("molecule-id", "value"),
               dash.dependencies.Input("snr", "value")],
              [dash.dependencies.State("db-location-access", "value")])
def update_molecule_plot(mol_id, snr, db_name):
    if mol_id and snr and db_name:
        snr = float(snr)
        mc = database.MoleculeConnector(db_name)
        mc.load_db_to_class()
        lab, back = mc.get_single_molecule_from_database(int(mol_id))
        snrx, snry = get_nicking_sites(lab, snr)
        mc.db.close()
        fig = molecule_plot_data(lab, back, snrx, snry)
        layout = dict(xaxis=dict(title="Molecule length (kb)"), yaxis=dict(title="Intensity"), title="Molecule Plot",
                      legend=dict(xanchor="center", yanchor="top"), paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)')
        return plot_graph(fig, "mol", layout=layout)


@app.callback(dash.dependencies.Output("snr-out", "children"),
              [dash.dependencies.Input("snr", "value")])
def show_snr_value(snr):
    if snr:
        snr = str(snr)
        return html.H6(snr, style={"text-align": "right"})
    else:
        return html.H6("")

@app.callback(dash.dependencies.Output("molecule-stats", "children"),
              [dash.dependencies.Input("run-id", "value"),
               dash.dependencies.Input("scan-id", "value"),
               dash.dependencies.Input("column-id", "value")],
              [dash.dependencies.State("db-location-access", "value")])
def update_molecule_stats_table(run_id, scan_id, column_id, db_name):
    if run_id and scan_id and column_id and db_name:
        mc = database.MoleculeConnector(db_name)
        mc.load_db_to_class()
        run_info = mc.db_runs[run_id]
        scan_id, column_id = int(scan_id), int(column_id)
        mc.db.close()
        num_molecules_in_run = 0
        num_molecules_in_scan = 0
        for scan in run_info:
            current_scan = run_info[scan]
            print(current_scan.keys())
            if scan == scan_id:
                for column in current_scan["columns"]:
                    current_column = current_scan["columns"][column]
                    if column == column_id:
                        num_molecules_in_col = current_column["molecule_count"]
                        num_molecules_in_scan += num_molecules_in_col
                        num_molecules_in_run += num_molecules_in_col
                        mean_len = current_column["mean_molecule_length"] //2
                        max_len = current_column["longest_shape"] //2
                    else:
                        num_molecules_in_scan += current_column["molecule_count"]
                        num_molecules_in_run += current_column["molecule_count"]
            else:
                for column in current_scan["columns"]:
                    current_column = current_scan["columns"][column]
                    num_molecules_in_run += current_column["molecule_count"]
        return [html.Br(), html.Div([html.H6("Molecules in Run: %s" % num_molecules_in_run, className="four columns"),
                html.H6("Molecules in Scan: %s" % num_molecules_in_scan, className="four columns"),
                html.H6("Molecules in Column: %s\n" % num_molecules_in_col, className="four columns")], className="container"),
                html.Div([html.H6("Mean molecule length in Column: %skb" % mean_len, className="six columns"),
                html.H6("Max molecule length in Column: %skb" % max_len, className="six columns"), html.Br()], className="container")]
    else:
        return ""

@app.callback(dash.dependencies.Output("db-summary-field", "children"),
              [dash.dependencies.Input("connect-db", "n_clicks")],
              [dash.dependencies.State("db-location-access", "value")])
def update_db_summary(clicked, db_name):
    if clicked and db_name:
        mc = database.MoleculeConnector(db_name)
        mc.load_db_to_class()
        run_info = mc.db_runs
        try:
            _ = mc.get_single_molecule_from_database(1)
        except AttributeError:
            mc.write_molecule_metadata_to_disk()
        mc.db.close()
        all_info = get_db_stats(run_info)
        return [html.Br(), html.Div([html.H6("Number of molecules: %s" % all_info[0], className="four columns", style={"text-align": "center"}),
                html.H6("Average length: %s" % all_info[1], className="four columns", style={"text-align": "center"}),
                html.H6("Longest molecule: %s" % all_info[2], className="four columns", style={"text-align": "center"})],
                        className="container")]
    return ""


@app.callback(dash.dependencies.Output("histogram-field", "hidden"),
              [dash.dependencies.Input("backbone-avg-button", "n_clicks_timestamp"),
               dash.dependencies.Input("backbone-max-button", "n_clicks_timestamp"),
               dash.dependencies.Input("label-density-button", "n_clicks_timestamp"),
               dash.dependencies.Input("mol-length-button", "n_clicks_timestamp")],
              [dash.dependencies.State("db-location-access", "value")])
def show_backbone_histogram(avg_t, max_t, lab_num, mol_len, db_name):
    if (max_t and db_name) or (avg_t and db_name) or (db_name and lab_num) or (db_name and mol_len):
        return False
    else:
        return True


@app.callback(dash.dependencies.Output("histogram-field", "children"),
              [dash.dependencies.Input("backbone-avg-button", "n_clicks_timestamp"),
               dash.dependencies.Input("backbone-max-button", "n_clicks_timestamp"),
               dash.dependencies.Input("label-density-button", "n_clicks_timestamp"),
               dash.dependencies.Input("mol-length-button", "n_clicks_timestamp")],
              [dash.dependencies.State("db-location-access", "value"),
               dash.dependencies.State("snr", "value")])
def show_backbone_histogram(avg_t, max_t, lab_num, mol_len, db_name, snr):
    if not snr:
        snr = 3.
    if not avg_t and max_t and lab_num:
        return ""
    if not avg_t:
        avg_t = 0
    if not max_t:
        max_t = 0
    if not lab_num:
        lab_num = 0
    if not mol_len:
        mol_len = 0
    if db_name:
        mc = database.MoleculeConnector(db_name)
        mc.load_db_to_class()
        if max_t > max(avg_t, lab_num, mol_len):
            backbones = (x[1] for x in mc.yield_molecule_signals_in_all_runs())
            hist_figures = histogram(backbones)
            mc.db.close()
            layout = dict(xaxis=dict(title="Maximum Intensity Bins"), yaxis=dict(title="Num. of Molecules"),
                          title="Maximum Backbone Intensity", bargap=0.2,
                          legend=dict(xanchor="center", yanchor="top"), paper_bgcolor='rgba(0,0,0,0)',
                          plot_bgcolor='rgba(0,0,0,0)')
            return plot_graph([hist_figures[1]], "backbone-hist", layout=layout)
        elif avg_t > max(lab_num, mol_len):
            backbones = (x[1] for x in mc.yield_molecule_signals_in_all_runs())
            hist_figures = histogram(backbones)
            mc.db.close()
            layout = dict(xaxis=dict(title="Average Intensity Bins"), yaxis=dict(title="Num. of Molecules"),
                          title="Average Backbone Intensity", bargap=0.2,
                          legend=dict(xanchor="center", yanchor="top"), paper_bgcolor='rgba(0,0,0,0)',
                          plot_bgcolor='rgba(0,0,0,0)')
            return plot_graph([hist_figures[0]], "backbone-hist", layout=layout)
        elif lab_num > mol_len:
            mols = [x[0] for x in mc.yield_molecule_signals_in_all_runs()]
            hist_figures = histogram_labels(mols, float(snr))
            mc.db.close()
            layout = dict(xaxis=dict(title="Number of labels"), yaxis=dict(title="Num. of Molecules"),
                          title="Number of specific labels per molecule", bargap=0.2,
                          legend=dict(xanchor="center", yanchor="top"), paper_bgcolor='rgba(0,0,0,0)',
                          plot_bgcolor='rgba(0,0,0,0)')
            return plot_graph([hist_figures], "backbone-hist", layout=layout)
        else:
            mols = [x[0] for x in mc.yield_molecule_signals_in_all_runs()]
            hist_figures = histogram_lens(mols)
            mc.db.close()
            layout = dict(xaxis=dict(title="Length (kb)"), yaxis=dict(title="Num. of Molecules"),
                          title="Molecule length distribution", bargap=0.2,
                          legend=dict(xanchor="center", yanchor="top"), paper_bgcolor='rgba(0,0,0,0)',
                          plot_bgcolor='rgba(0,0,0,0)')
            return plot_graph([hist_figures], "backbone-hist", layout=layout)



@app.callback(dash.dependencies.Output("link-field", "children"),
              [dash.dependencies.Input("bnx-convert-button", "n_clicks")],
              [dash.dependencies.State("db-location-bnx", "value"),
               dash.dependencies.State("filename", "value"),
               dash.dependencies.State("snr-bnx-input", "value"),
               dash.dependencies.State("min-len-bnx-input", "value"),
               dash.dependencies.State("intensity-bnx-input", "value")])
def convert_bnx(clicked, db_name, filename, snr, min_len, intensity):
    def filter_func(sig):
        if (sig[0].shape[0] > (float(min_len) * 2)) and (np.mean(sig[1]) < float(intensity)):
            return True
        else:
            return False
    if clicked and db_name and filename and snr and min_len and intensity:
        mc = database.MoleculeConnector(db_name)
        mc.load_db_to_class()
        database.molecules_to_bnxv2(filter(filter_func, mc.yield_molecule_signals_in_all_runs()), 10, 510, "static/%s" % filename,
                                    bnx_template_path="bnx_template.txt", signal_to_noise_ratio=float(snr))
        mc.db.close()
        return html.A("Download %s here" % filename, href="/static/%s" % filename, style={"margin": "10%"})
    else:
        return html.A("")

@app.callback(dash.dependencies.Output("db-location-bnx", "options"),
              [dash.dependencies.Input("connect-db", "n_clicks")],
              [dash.dependencies.State("db-location-access", "options")])
def update_bnx_db_selection(clicked, location_val):
    if clicked and location_val:
        return location_val
    else:
        return ""

@app.callback(dash.dependencies.Output("db-location-access", "options"),
              [dash.dependencies.Input("refresh", "n_clicks")])
def update_db_dropdown(clicked):
    if clicked:
        database_names = os.listdir("database")
        return [{"label": x, "value": "database/%s" % x} for x in database_names if x.endswith(".db")]
    database_names = os.listdir("database")
    return [{"label": x, "value": "database/%s" % x} for x in database_names if x.endswith(".db")]

@app.server.route('/static/<path:path>')
def downlad_file(path):
    root_dir = os.getcwd()
    return flask.send_from_directory(
        os.path.join(root_dir, 'static'), path)




if __name__ == '__main__':
    from sys import argv
    app.run_server(debug=True, host=argv[1], port=argv[2])
