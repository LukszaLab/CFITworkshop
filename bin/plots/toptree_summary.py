import json
import os
import argparse
import plotly.graph_objects as go
import plotly.express as px

#from cfit.plot.Clone import Clone
from cfit.plot.CloneTree import CloneTree

from cfit.util.Analysis import Analysis

# How to run:
# hdir=$HOME/MtSinai/Data_met_fs_corrected
# tree_format='pairtree'
# odir=$HOME/MtSinai/Data_met_fs_corrected/Reports
# python3 bin/toptree_summary.py -dir $hdir -tree_format $tree_format -out_dir $odir

################################################################################

AXIS = dict(showline=False, 
            zeroline=False,
            showgrid=False,
            showticklabels=False)

CIRCLE_SIZE = 28
CIRCLE_COLOR = '#818589'

LABEL_SIZE = 12

STAR_SIZE = 20

LINE_WIDTH = 2.5
LINE_COLOR = '#1D2339'

FONT_SIZE = 14

LINE = dict(color=LINE_COLOR, width=LINE_WIDTH)

MARKER = dict(symbol='circle', size=CIRCLE_SIZE, color=CIRCLE_COLOR,
                line=dict(color=LINE_COLOR, width=2))


####################################################################################

def get_parent_vector(cfit_tree):
    n = len(cfit_tree.nodes)
    result = []
    for i in range(1, n):
        result += [cfit_tree.nodes[i].parent.id]
    return result, n


def create_nodes_depth_dict(parent_vector, n):
    nid2depth = dict()
    nid2depth[0] = 0
    for i in range(1, n):
        d = 1
        parent = parent_vector[i-1]
        while parent != 0:
            d += 1
            parent = parent_vector[parent-1]
        nid2depth[i] = d
    return nid2depth

def initialize_clonetree(parent_vector, n, nid2depth):
    nid2node = dict()
    root = Clone(id=0, depth=0)
    root.set_y(0)
    root.set_pid(0)
    root.set_parent(root)
    nid2node[0] = root
    nodes = [root]
    for i in range(n-1):
        new_node = Clone(id=i+1)
        new_node.set_pid(parent_vector[i])
        new_node.set_depth(float(nid2depth[i+1]/4))
        nid2node[i+1] = new_node
        nodes.append(new_node)
    for node in nodes:
        node.set_parent(nid2node[node.pid])
        parent = nid2node[node.pid]
        if parent != node:
            parent.add_child(node)
    for node in nodes:
        node.set_subtree_height()
    leaves = list()
    for node in nodes:
        node.children_ordered = sorted(node.children, key=lambda x: x.subtree_h, reverse=True)
        if len(node.children) == 0:
            leaves.append(node)
    tree = CloneTree(root=root, nodes=nodes, leaves=leaves, h=max(nid2depth.values()))
    tree.nid2node = nid2node
    tree.order_leaves()
    tree.set_ys_of_clones()
    return tree


def find_nodes_with_drivers(tree, driver_genes):
    #print(driver_genes)
    gene_nids = {gene: None for gene in driver_genes}
    for nid in tree.nodes:
        for gene in driver_genes:
            if gene in [mutation.gene for mutation in tree.nodes[nid].exclusiveMutations]:
                gene_nids[gene] = nid
    driver_nids = dict()
    for gene in gene_nids.keys():
        if gene_nids[gene] == None:
            continue
        if gene_nids[gene] not in driver_nids.keys():
            driver_nids[gene_nids[gene]] = [gene]
        else:
            driver_nids[gene_nids[gene]] += [gene]
    return driver_nids


def get_general_params(tree, n):
    colors = px.colors.sample_colorscale("turbo", [k/(n - 1) for k in range(n)])
    edges = []
    for nid in tree.nodes:
        edges.append((tree.nodes[nid].parent.id, nid))
    return colors, edges


def get_edge_coords(X, Y, edges):
    Xe, Ye = [], []
    for edge in edges:
        Xe += [X[edge[0]], X[edge[1]], None]
        Ye += [Y[edge[0]], Y[edge[1]], None]
    return Xe, Ye


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
        if y[0] < 0.03: # clone frequency threshold
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
    fig.update_layout(title='Cluster Depth vs Mutations', 
                      xaxis_title="Distance from root", 
                      yaxis_title='Number of mutations',
                      showlegend=True)
    return fig


def plot_frequencies(thepatient, edges, depths, n_nodes, colors, which_freq='inc', n_trees=2):
    figs = []
    samples = {thepatient.samples[i]: thepatient.samples[i].trees[0] for i in range(n_trees)}
    for sample in samples:
        samptree = samples[sample]
        Y = {samnode.id: samnode.X for samnode in samptree.nodes.values()}
        y_title = 'Inclusive frequency'
        if which_freq == 'exc':
            Y = {samnode.id: samnode.Y for samnode in samptree.nodes.values()}
            y_title = 'Exclusive frequency'
        Xe, Ye = get_edge_coords(depths, Y, edges)
        title = "Cluster Depth vs Frequency (sample "+sample.name+")<br><sup>Sample Purity: {:.2f}</sup><br><sup>Tissue: {}</sup>".format(samptree.purity, sample.tissue)
        fig = plot_general(n_nodes, depths, Y, Xe, Ye, colors)
        fig.update_layout(title=title, xaxis_title="Distance from root", yaxis_title=y_title,
                            showlegend=True)
        figs.append(fig)
    return figs

def plot_mut_freqs(thepatient, edges, n_nodes, colors):
    sample = thepatient.samples[0]
    tree, samptree = thepatient.trees[0], sample.trees[0]
    X = {node.id: len(node.mutations) for node in tree.nodes.values()}
    Y = {samnode.id: samnode.X for samnode in samptree.nodes.values()}
    y_title = "Inclusive frequency"
    x_title = "Number of mutations"
    title = "Mutations vs Frequency (top tree for sample "+sample.name+")<br><sup>Sample Purity: {:.2f}</sup>".format(samptree.purity)
    Xe, Ye = get_edge_coords(X, Y, edges)
    fig = plot_general(n_nodes, X, Y, Xe, Ye, colors)
    fig.update_layout(title=title, xaxis_title=x_title, yaxis_title=y_title,
                            showlegend=True)
    return fig


def plot_clonetree(patient, clonetree, n, driver_nodes):
    pos = dict()
    edges = []
    for node in clonetree.nodes:
        pos[node.id] = (node.depth, node.y)
        if node.id == 0:
            continue
        edges.append((node.pid, node.id))

    Xe, Ye = [], []
    for edge in edges:
        Xe += [pos[edge[0]][0], pos[edge[1]][0], None]
        Ye += [pos[edge[0]][1], pos[edge[1]][1], None]

    labels = [str(i) for i in range(n+1)]

    annotations = []
    for i in range(n):
        annotations.append(
            dict(text=labels[i], x=pos[i][0], y=pos[i][1], xref='x1', yref='y1',
            font=dict(color='#FFFFFF', size=LABEL_SIZE), showarrow=False))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Xe, y=Ye, mode='lines',line=dict(color=LINE_COLOR, width=LINE_WIDTH), 
                             hoverinfo='none', showlegend=False))
    fig.add_trace(go.Scatter(x=[pos[i][0] for i in range(n)], 
                             y=[pos[i][1] for i in range(n)], 
                             mode='markers', 
                             marker=dict(symbol='circle', size=CIRCLE_SIZE, color=CIRCLE_COLOR,
                             line=dict(color=LINE_COLOR, width=1)), 
                             hoverinfo='none', showlegend=False))
    
    m = len(driver_nodes)
    colors = px.colors.sample_colorscale("plasma", m+1, low=0.2, high=0.8)
    star_colors = dict(zip(driver_nodes.keys(), colors))
    
    for node in driver_nodes.keys():
        name = str(driver_nodes[node])
        name = name[1:-1]
        name = name.replace("'", "") 
        (x1, y1) = pos[node]
        (x2, y2)= pos[clonetree.nid2node[node].pid]
        star = ((x1+x2)/2, (y1+y2)/2)
        
        fig.add_trace(go.Scatter(name=name, mode="markers", 
                                 x=[star[0]], y=[star[1]], 
                                 marker_symbol="star", 
                                 marker_color=star_colors[node],
                                 marker_size=STAR_SIZE, hoverinfo='none'))
    
        
    fig.update_layout(title="Top Clone Tree", annotations=annotations, font_size=FONT_SIZE, 
                      legend_title_text='Mutation', showlegend=True, margin=dict(l=50,r=50,b=50,t=60,pad=5),
                        xaxis=AXIS, yaxis=AXIS)
    return fig



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-dir", help="data directory") 
                            # contains:
                            #    mapping.json
                            #    config.json
                            #    VCF data folder
                            #    tree outputs folder
                            #    driver_genes.txt (optional)

    parser.add_argument("-tree_format", default="phylowgs",
                        help="tree format (phylowgs, pairtree, mclust)")
    
    parser.add_argument("-out_dir", 
                        help="path to output directory for html file(s)")
    
    args = parser.parse_args()

    mappingfile = os.path.join(args.dir, "mapping.json") # can change this path to mapping file with specific subset of patients
    configfile = os.path.join(args.dir, "config.json")

    with open(mappingfile, 'r') as f:
        mapping = json.load(f)
    
    with open(configfile, 'r') as f:
        config = json.load(f)
    
    driversfile, drivers = "", ""

    # if driver_genes.txt file present, clone tree will indicate branches along which driver genes acquired mutations
    if os.path.exists(os.path.join(args.dir, "driver_genes.txt")):
        driversfile = os.path.join(args.dir, "driver_genes.txt")
        with open(driversfile, 'r') as f:
            drivers = [line.strip() for line in f.readlines()]
   

    
    anl = Analysis()
    anl.initialize_config(config, mapping, args.dir, kd_thr=500, ns=[9], tree_format=args.tree_format)


    for patient in mapping: # generates graphs report for all patients in mapping file
                            # can edit mapping file to indicate specific subset of patients
        patname = patient["name"]
        thepatient = anl.patients[patname]

        parent_vector, n = get_parent_vector(thepatient.trees[0]) # extract tree topology in order to initialize clonetree for plotting
        nid2depth = create_nodes_depth_dict(parent_vector, n) # --> dict: key = node id, value = distance from root
        clonetree = initialize_clonetree(parent_vector, n, nid2depth)
        
        driver_nodes = find_nodes_with_drivers(thepatient.trees[0], drivers)

        # topology plot:
        clonetree_plot = plot_clonetree(patname, clonetree, n, driver_nodes)
        
        colors, edges = get_general_params(thepatient.trees[0], n)
        m = len(thepatient.samples)

        # other plots
        mut_plot = plot_mutations(thepatient, edges, nid2depth, n, colors)
        inc_freq_plots = plot_frequencies(thepatient, edges, nid2depth, n, colors, 'inc', m)
        exc_freq_plots = plot_frequencies(thepatient, edges, nid2depth, n, colors, 'exc', m)
        mut_freq_plot = plot_mut_freqs(thepatient, edges, n, colors)


        # export interactive plotly graphs to html file
        with open(os.path.join(args.out_dir, patname+'_'+args.tree_format+'_report.html'), 'a') as f:
            f.write("<html> <body> <h1>{} Graphs Report for \
           <font color = #000000>{}</font></h1>\n \
           </body></html>".format(args.tree_format, patname))
            f.write(clonetree_plot.to_html(full_html=False, include_plotlyjs='cdn'))
            f.write(mut_plot.to_html(full_html=False, include_plotlyjs='cdn'))
            for fig in inc_freq_plots:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
            for fig in exc_freq_plots:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
            f.write(mut_freq_plot.to_html(full_html=False, include_plotlyjs='cdn'))


if __name__ == "__main__":
    main()




    


