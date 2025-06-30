import json
import os
import plotly.graph_objects as go
import plotly.express as px

from cfit.tree.mutation.Neoantigen import Neoantigen
from cfit.util.Analysis import Analysis
from cfit.util.Log import Log
from cfit.util.Utils import Utils
from cfit.tree.Tree import Tree
from cfit.tree.node.Node import Node
from cfit.tree.SampleTree import SampleTree

#  Run as python3 bin/plot_tree.py -d $HOME/MtSinai/Data -ns 9 -kd_thr 500

def initialize(tree_format):

    parser = Utils.make_cohort_parser()
    args = parser.parse_args()

    if "-" in args.ns:
        ns = [int(x) for x in  args.ns.split("-")[:2]]
        ns = list(range(ns[0], ns[1]+1))
    else:
        ns = [int(x) for x in (args.ns).split(",")]

    Neoantigen.WEPS = 0.0
    npos = "1-9"

    Utils.w = 1
    anl = Analysis()
    anl.set_MHC_version(args.netMHC)

    if args.mapping is None:
        mappingfile = os.path.join(args.dir, "Info_files", "mapping.json")
    else:
        mappingfile = args.mapping
    with open(mappingfile) as f:
        mappingjs = json.load(f)

    if args.config is None:
        if tree_format == 'pairtree':
            configfile = os.path.join(args.dir, "Info_files", "config_pairtree.json")
        else:
            configfile = os.path.join(args.dir, "Info_files", "config_phylowgs.json")
    else:
        configfile = args.config
    with open(configfile) as f:
        configjs = json.load(f)

    alndir = os.path.join(args.dir, configjs["aln_dir"])
    iedbfasta = None
    if "iedb_file" in configjs:
        iedbfasta = os.path.join(args.dir, configjs["iedb_file"])
    else:
        iedbfasta = os.path.join(alndir, "enemy.fasta")

    Log.logger("Initializing...")
    Log.logger("iedbfasta = "+iedbfasta)
    
    if tree_format == 'pairtree':
        anl.initialize_config(configjs, mappingjs, args.dir, kd_thr=args.kd_thr, ns=ns, tree_format='pairtree')
    else:
        anl.initialize_config(configjs, mappingjs, args.dir, kd_thr=args.kd_thr, ns=ns, tree_format='phylowgs')
    return anl


def get_edges(tree):
    edges = []
    for nid in tree.nodes:
        edges.append((tree.nodes[nid].parent.id, nid))
    return edges


def get_node_depths(tree):
    node_depth = dict()
    for nid in tree.nodes:
        if nid == 0:
            node_depth[nid] = 0
            continue
        d = 1
        parent_id = tree.nodes[nid].parent.id
        while parent_id != 0:
            d += 1
            parent_id = tree.nodes[parent_id].parent.id
        node_depth[nid] = d
    return node_depth


def get_edge_coords(X, Y, edges):
    Xe, Ye = [], []
    for edge in edges:
        Xe += [X[edge[0]], X[edge[1]], None]
        Ye += [Y[edge[0]], Y[edge[1]], None]
    return Xe, Ye


def get_general_params(tree):
    n = len(tree.nodes)
    colors = px.colors.sample_colorscale("turbo", [k/(n - 1) for k in range(n)])
    edges = get_edges(tree)
    depths = get_node_depths(tree)
    return n, colors, edges, depths


def plot_general(n, X, Y, Xe, Ye, colors):
    assert(len(colors) == n)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Xe, y=Ye, mode='lines',
                    line=dict(color='#000000', width=3), hoverinfo=None, showlegend=False))
    for i in range(n):
        x = [X[i]]
        y = [Y[i]]
        size = 30
        color = colors[i] 
        if y[0] < 0.03:
            size = 10
            color = 'black'
        fig.add_trace(go.Scatter(name="clone "+str(i), x=x, y=y, mode='markers',
                             marker=dict(symbol='circle', size=size, color=color,
                                        line=dict(color='#000000', width=1))))
    return fig


def plot_mutations(thepatient, edges, depths, n, colors):
    tree = thepatient.trees[0]
    Y = {nid: len(tree.nodes[nid].mutations) for nid in tree.nodes}
    Xe, Ye = get_edge_coords(depths, Y, edges)
    fig = plot_general(n, depths, Y, Xe, Ye, colors)
    fig.update_layout(title='Cluster Depth vs Mutations'+' ('+thepatient.name+')', 
                      xaxis_title="Distance from root", yaxis_title='Number of mutations',
                        showlegend=True)
    fig.show()

   
def plot_frequencies(thepatient, edges, depths, n_nodes, colors, which_freq='inc', n_trees=2):
    samples = {thepatient.samples[i].name: thepatient.samples[i].trees[0] for i in range(n_trees)}
    for sample in samples:
        samptree = samples[sample]
        Y = {samnode.id: samnode.X for samnode in samptree.nodes.values()}
        y_title = 'Inclusive frequency'
        if which_freq == 'exc':
            Y = {samnode.id: samnode.Y for samnode in samptree.nodes.values()}
            y_title = 'Exclusive frequency'
        Xe, Ye = get_edge_coords(depths, Y, edges)
        title = "Cluster Depth vs Frequency ("+thepatient.name+")<br><sup>Sample Name: {}</sup><br><sup>Sample Purity: {:.2f}</sup>".format(sample, samptree.purity)
        fig = plot_general(n_nodes, depths, Y, Xe, Ye, colors)
        fig.update_layout(title=title, xaxis_title="Distance from root", yaxis_title=y_title,
                            showlegend=True)
        fig.show()


def plot_mut_freq(thepatient, edges, n_nodes, colors):
    sample = thepatient.samples[0]
    tree, samptree = thepatient.trees[0], sample.trees[0]
    X = {node.id: len(node.mutations) for node in tree.nodes.values()}
    Y = {samnode.id: samnode.X for samnode in samptree.nodes.values()}
    y_title = "Inclusive frequency"
    x_title = "Number of mutations"
    title = "Mutations vs Frequency ("+thepatient.name+")<br><sup>Sample Name: {}</sup><br><sup>Sample Purity: {:.2f}</sup>".format(sample.name, samptree.purity)
    Xe, Ye = get_edge_coords(X, Y, edges)
    fig = plot_general(n_nodes, X, Y, Xe, Ye, colors)
    fig.update_layout(title=title, xaxis_title=x_title, yaxis_title=y_title,
                            showlegend=True)
    fig.show()


def main():
    tree_format = 'pairtree'
    anl = initialize(tree_format)
    patients = ['11-LTS']
    for patient in patients:
        thepatient = anl.patients[patient]
        n, colors, edges, depths = get_general_params(thepatient.trees[0])
        plot_mutations(thepatient, edges, depths, n, colors)
        plot_frequencies(thepatient, edges, depths, n, colors, 'inc', 2)
        plot_mut_freq(thepatient, edges, n, colors)
    


if __name__ == "__main__":
    main()