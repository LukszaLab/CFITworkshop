import json
import os
import plotly.graph_objects as go
import plotly.express as px

from Clone import Clone
from CloneTree import CloneTree


from cfit.tree.mutation.Neoantigen import Neoantigen
from cfit.util.Analysis import Analysis
from cfit.util.Log import Log
from cfit.util.Utils import Utils
from cfit.tree.Tree import Tree
from cfit.tree.node.Node import Node
from cfit.tree.SampleTree import SampleTree

#  Run as python3 bin/tree_report.py -d $HOME/MtSinai/Data_met_fs_corrected -ns 9 -kd_thr 500

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
TITLE = 'Clone Tree'

LINE = dict(color=LINE_COLOR, width=LINE_WIDTH)
MARKER = dict(symbol='circle', size=CIRCLE_SIZE, color=CIRCLE_COLOR,
                line=dict(color=LINE_COLOR, width=2))


####################################################################################

def initialize_cfit(tree_format):

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
        mappingfile = os.path.join(args.dir, "mapping.json")
    else:
        mappingfile = args.mapping
    with open(mappingfile) as f:
        mappingjs = json.load(f)

    if args.config is None:
       configfile = os.path.join(args.dir, "config.json")
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


def get_parent_vector(tree):
    result = []
    for i in range(1, len(tree.nodes)):
        result += [tree.nodes[i].parent.id]
    return result, i

def create_nodes_depth_dict(parent_vector, n):
    nid2depth = dict()
    nid2depth[0] = 0
    for i in range(1,n+1):
        d = 1
        parent = parent_vector[i-1]
        while parent != 0:
            d += 1
            parent = parent_vector[parent-1]
        nid2depth[i] = d
    return nid2depth


def initialize_plottree(parent_vector, n, nid2depth):
    nid2node = dict()
    root = Clone(id=0, depth=0)
    root.set_y(0)
    root.set_pid(0)
    root.set_parent(root)
    nid2node[0] = root
    nodes = [root]
    for i in range(n):
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
    return tree

def order_leaves(tree):
    leafid2ord = dict()
    done = False
    d = -1
    curr_node = tree.root
    while not done:
        d += 1
        while not curr_node in tree.leaves and curr_node.children_ordered:
            curr_node = curr_node.children_ordered.pop(0)
        if curr_node in tree.leaves:
            leafid2ord[curr_node.id] = d
        else:
            d -= 1
        if curr_node.is_root():
            done = True
        else:
            curr_node = curr_node.parent
    tree.leaf2ord = leafid2ord


def get_postorder(tree):
    result = []
    to_visit = [tree.root]
    while len(to_visit) > 0:
        curr_node = to_visit.pop()
        result.append(curr_node)
        for child in curr_node.children:
            to_visit.append(child)
    return result[::-1]


def set_ys_of_clones(tree, order):
    for leaf in tree.leaves:
        leaf.set_y(tree.leaf2ord[leaf.id] * CIRCLE_SIZE)
    for node in order:
        parent = node.parent
        y = 0
        for child in parent.children:
            y += child.y
        y /= len(parent.children)
        parent.set_y(y)


def clone_coordinates(tree, n):
    pos = dict()
    for node in tree.nodes:
        pos[node.id] = (node.depth, node.y)
    assert(len(pos) == n+1)
    X = [pos[i][0] for i in range(n+1)]
    Y = [pos[i][1] for i in range(n+1)]
    return pos, X, Y


def clonetree_edges(tree, pos):
    edges = []
    for node in tree.nodes:
        if node.id == 0:
            continue
        edges.append((node.pid, node.id))
    Xe = []
    Ye = []
    for edge in edges:
        Xe += [pos[edge[0]][0], pos[edge[1]][0], None]
        Ye += [pos[edge[0]][1], pos[edge[1]][1], None]
    return Xe, Ye


def clonetree_mut(tree, mut_node, pos):
    (x1, y1) = pos[mut_node]
    (x2, y2)= pos[tree.nid2node[mut_node].pid]
    return ((x1+x2)/2, (y1+y2)/2)


def plot_clonetree(patient, tree, n, driver_nodes):
    pos, X, Y = clone_coordinates(tree, n) # make sure all coord lists are consistent in ordering!
    Xe, Ye = clonetree_edges(tree, pos)
    labels = [str(i) for i in range(n+1)]
    annotations = make_annotations(X, Y, labels, font_size=LABEL_SIZE, font_color='#FFFFFF')
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Xe, y=Ye, mode='lines',line=dict(color=LINE_COLOR, width=LINE_WIDTH), 
                             hoverinfo='none', showlegend=False))
    fig.add_trace(go.Scatter(x=X, y=Y, mode='markers', 
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
        star = clonetree_mut(tree, node, pos)
        
        fig.add_trace(go.Scatter(name=name, mode="markers", x=[star[0]], y=[star[1]], marker_symbol="star", 
                                 marker_color=star_colors[node],
                                    marker_size=STAR_SIZE, hoverinfo='none'))
    
        
    fig.update_layout(title=TITLE+": "+patient, annotations=annotations, font_size=FONT_SIZE, 
                      legend_title_text='Mutation', showlegend=True, margin=dict(l=50,r=50,b=50,t=60,pad=5),
                        xaxis=AXIS, yaxis=AXIS)
    return fig


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
    #assert(n == len(colors))
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
    fig.update_layout(title='Cluster Depth vs Mutations', 
                      xaxis_title="Distance from root", yaxis_title='Number of mutations',
                        showlegend=True)
    return fig
    # fig.show()

   
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
        title = "Cluster Depth vs Frequency ("+sample.name+")<br><sup>Sample Purity: {:.2f}</sup><br><sup>Tissue: {}</sup>".format(samptree.purity, sample.tissue)
        fig = plot_general(n_nodes, depths, Y, Xe, Ye, colors)
        fig.update_layout(title=title, xaxis_title="Distance from root", yaxis_title=y_title,
                            showlegend=True)
        figs.append(fig)
    return figs


def plot_mut_freq(thepatient, edges, n_nodes, colors):
    sample = thepatient.samples[0]
    tree, samptree = thepatient.trees[0], sample.trees[0]
    X = {node.id: len(node.mutations) for node in tree.nodes.values()}
    Y = {samnode.id: samnode.X for samnode in samptree.nodes.values()}
    y_title = "Inclusive frequency"
    x_title = "Number of mutations"
    title = "Mutations vs Frequency ("+sample.name+")<br><sup>Sample Purity: {:.2f}</sup>".format(samptree.purity)
    Xe, Ye = get_edge_coords(X, Y, edges)
    fig = plot_general(n_nodes, X, Y, Xe, Ye, colors)
    fig.update_layout(title=title, xaxis_title=x_title, yaxis_title=y_title,
                            showlegend=True)
    return fig


def make_annotations(X, Y, text, font_size, font_color):
    L = len(X)
    assert (L == len(Y))
    if len(text) != L:
        raise ValueError('pos and text must have same length')
    annotations = []
    for i in range(L):
        annotations.append(
            dict(text=text[i], x=X[i], y=Y[i], xref='x1', yref='y1',
            font=dict(color=font_color, size=font_size), showarrow=False))
    return annotations


def find_nodes_with_drivers(tree, driver_genes):
    print(driver_genes)
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


def main():
    #TODO: optimize 
    driver_genes = []

    with open("/Users/veramazeeva/MtSinai/Data_met_fs_corrected/driver_genes.txt", "r") as f:
        driver_genes = [line.strip() for line in f.readlines()]
        
    print(driver_genes)
    tree_format = 'pairtree'
    anl = initialize_cfit(tree_format)
    
    patients = ['PAM41']
    for patname in patients:
        thepatient = anl.patients[patname]
        parent_vector, n = get_parent_vector(thepatient.trees[0])
        nid2depth = create_nodes_depth_dict(parent_vector, n)
        clonetree = initialize_plottree(parent_vector, n, nid2depth)
        order_leaves(clonetree)
        postorder = get_postorder(clonetree)
        
        set_ys_of_clones(clonetree, postorder)
        driver_nodes = find_nodes_with_drivers(thepatient.trees[0], driver_genes)
        print(driver_nodes)
        n1, colors, edges, depths = get_general_params(thepatient.trees[0])
        m = len(thepatient.samples)
        
        fig = plot_clonetree(patname, clonetree, n, driver_nodes)
        fig1 = plot_mutations(thepatient, edges, depths, n1, colors)
        figs = plot_frequencies(thepatient, edges, depths, n1, colors, 'inc', m) + plot_frequencies(thepatient, edges, depths, n1, colors, 'exc', m)
        #fig2 = plot_mut_freq(thepatient, edges, n, colors)
        
        with open("reports/"+tree_format+'/'+patname+'_report.html', 'a') as f:
            f.write("<html> <body> <h1>{} Graphs Report for \
           <font color = #000000>{}</font></h1>\n \
           </body></html>".format(tree_format, patname))
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
            f.write(fig1.to_html(full_html=False, include_plotlyjs='cdn'))
            for fig in figs:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
            #f.write(fig2.to_html(full_html=False, include_plotlyjs='cdn'))
            f.close()
    


if __name__ == "__main__":
    main()