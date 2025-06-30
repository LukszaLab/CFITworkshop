import plotly.graph_objects as go
import plotly.express as px

from cfit.plot.CloneTree import CloneTree


class PlotTree():

    """
        PlotTree object stores node coordinates and edges used to plot tree from scatter plot.
    
    Attributes:

        what: str 
            what is plotted

        pos: dict (clone id --> (x-coord, y-coord))
            used to look up position of node by id, which is used to plot node point within scatter plot

        edges: list
            list of tuples (node.id, node parent.id)

        n: int
            number of nodes within the tree

        tree: cfit.Tree object
            contains topology of tree
        
        figs: list
            list of go.Figure() objects

    """


    def __init__(self, patient, what, drivers=None, tid=0, patientid=True):
        
        """
        Params
        ----------

        patient : PatientLine
            patient whose tree will be graphed

        what : str
            indicates what is plotted
            options: "topology", "dist vs incmuts", "dist vs excmuts", "dist vs neos", "muts vs freq", "dist vs incfreq", "dist vs excfreq"

        drivers : list
            list of driver genes
            only applicable for what="topology"
            if specified, will mark branches of clonetree where driver gene mutations are acquired
        
        tid : int
            tree index

        patientid : bool
            indicates whether or not to add patient identifier as part of plot title

        """
        
        if what not in ["topology", "dist vs incmuts", "dist vs excmuts", "muts vs freq", "dist vs incfreq", "dist vs excfreq", "dist vs neos"]:
            raise Exception("please specify valid plot type")

        self.thepatient = patient
        self.what = what
    
        self.tree = patient.trees[tid]
        self.n = len(patient.trees[0].nodes)

        self.edges = [(self.tree.nodes[nid].parent.id, nid) for nid in self.tree.nodes]
        self.pos = dict()

        self.labels = [str(i) for i in range(self.n)]
        
        self.annotations = list()
        self.drivers = drivers

        self.colors = px.colors.sample_colorscale("turbo", self.n)

        self.figs = list()

        self.title = ""
        self.xtitle = ""
        self.ytitle = ""

        name = ""
        if patientid:
            name = "("+self.thepatient.name


        if what == "topology":
            self.title = "Clone Tree" +str(tid+1) +" " +name+")"
            self.__init_from_clonetree()

        elif what == "dist vs incmuts":
            self.title = "Clone Depth vs Mutations " +name+" Tree"+str(tid+1)+")"
            self.__init_from_sampdata(x='dist', y='incmuts')

        elif what == "dist vs excmuts":
            self.title = "Clone Depth vs Mutations " +name+" Tree"+str(tid+1)+")"
            self.__init_from_sampdata(x='dist', y='excmuts')

        elif what == "dist vs neos":
            self.title = "Clone Depth vs Neoantigens " +name+" Tree"+str(tid+1)+")"
            self.__init_from_sampdata(x='dist', y='neos')

        elif what == "muts vs freq":
            self.title = "Clone Mutations vs Frequency " +name+" Tree"+str(tid+1)+")"
            self.__init_from_sampdata(x='incmuts', y='incfreq')

        elif what == "dist vs incfreq":
            self.title = "Clone Depth vs Frequency " +name+" Tree"+str(tid+1)+")"
            self.__init_from_sampdata(x='dist', y='incfreq')

        elif what == "dist vs excfreq":
            self.title = "Clone Depth vs Frequency " +name+" Tree"+str(tid+1)+")"
            self.__init_from_sampdata(x='dist', y='excfreq')

        

    def __init_from_clonetree(self):
        parent_vector = []
        for i in range(1, self.n):
            parent_vector += [self.tree.nodes[i].parent.id]
        clonetree = CloneTree(parent_vector)
        for node in clonetree.nid2node.values():
            self.pos[node.id] = (node.depth, node.y)
        for i in range(self.n):
            self.annotations.append(
                dict(text=self.labels[i], x=self.pos[i][0], y=self.pos[i][1], xref='x1', yref='y1',
                font=dict(color='#FFFFFF', size=12), showarrow=False))
        if self.drivers:
            self.get_driver_nodes()
    


    def get_driver_nodes(self):
        gene_nids = {gene: None for gene in self.drivers}
        for nid in self.tree.nodes:
            for gene in self.drivers:
                if gene in [mutation.gene for mutation in self.tree.nodes[nid].exclusiveMutations]:
                    gene_nids[gene] = nid
        self.drivers = dict()
        for gene in gene_nids:
            if gene_nids[gene] == None:
                continue
            if gene_nids[gene] not in self.drivers.keys():
                self.drivers[gene_nids[gene]] = [gene]
            else:
                self.drivers[gene_nids[gene]] += [gene]



    def __init_from_sampdata(self, x, y):

        if self.what in ["dist vs incmuts", "dist vs excmuts", "dist vs neos"]:
            
            pos_x, pos_y = dict(), dict()

            pos_x = self.get_dists()
            self.xtitle = 'Distance from root'

            if 'incmuts' in y:
                    pos_y = self.get_muts('inc')
                    self.ytitle = 'Cumulative number of mutations'
            
            elif 'excmuts' in y:
                    pos_y = self.get_muts('exc')
                    self.ytitle = 'Number of exclusive mutations'

            elif 'neos' in y:
                    pos_y = self.get_neos()
                    self.ytitle = 'Cumulative number of neoantigens'
            
            for i in range(self.n):
                self.pos[i] = [(pos_x[i], pos_y[i])]
        
        else:
            X, Y = [], []
            for sample in self.thepatient.samples:
                pos_x, pos_y = dict(), dict()

                if x == 'dist':
                    pos_x = self.get_dists()
                    self.xtitle = 'Distance from root'

                elif 'muts' in x:
                    pos_x = self.get_muts('inc')
                    self.xtitle = 'Cumulative number of mutations'
                
                pos_y = self.get_freqs(sample, y[0:3])
                self.ytitle = y[0:3].capitalize() + 'lusive frequency'
     
                X.append(pos_x)
                Y.append(pos_y)

            for i in range(self.n):
                self.pos[i] = [(X[j][i], Y[j][i]) for j in range(len(self.thepatient.samples))]


    def get_dists(self):
        result = dict()
        result[0] = 0
        for node in self.tree.nodes.values():
            if node.id == 0:
                continue
            d = 1
            parent = node.parent
            while parent.id != 0:
                d += 1
                parent = parent.parent
            result[node.id] = d
        return result


    def get_muts(self, mut):
        result = dict()
        if mut == 'exc':
            result = {node.id: len(node.exclusiveMutations) for node in self.tree.nodes.values()}
        else:
            result = {node.id: len(node.mutations) for node in self.tree.nodes.values()}
        return result


    def get_neos(self):
        result = {node.id: len(node.neoantigens) + len(node.fsneoantigens) for node in self.tree.nodes.values()}
        return result
    

    def get_freqs(self, sample, freq):
        result = dict()
        if freq == 'exc':
            result = {sampnode.id: sampnode.Y for sampnode in sample.trees[0].nodes.values()}
        else:
            result = {sampnode.id: sampnode.X for sampnode in sample.trees[0].nodes.values()}
        return result


    def plot(self):
        if self.what == 'topology':
            self.plot_clonetree()
        elif self.what in ['dist vs incmuts', 'dist vs excmuts', 'dist vs neos']:
            self.plot_dist()
        else:
            self.plot_freqs()
            

    def plot_clonetree(self):
        fig = go.Figure()
            
        Xe, Ye = [], []

        for edge in self.edges:
            Xe += [self.pos[edge[0]][0], self.pos[edge[1]][0], None]
            Ye += [self.pos[edge[0]][1], self.pos[edge[1]][1], None]
            

        fig.add_trace(go.Scatter(x=Xe, y=Ye, mode='lines',line=dict(color='#1D2339', width=3), 
                                 hoverinfo='none', showlegend=False))
            
        fig.add_trace(go.Scatter(x=[self.pos[i][0] for i in range(self.n)], 
                                 y=[self.pos[i][1] for i in range(self.n)], 
                                 mode='markers', 
                                 marker=dict(symbol='circle', size=45, color='#818589',
                                 line=dict(color='#1D2339', width=1)), 
                                 hoverinfo='none', showlegend=False))
        
        if self.drivers:
            m = len(self.drivers)
            colors = px.colors.sample_colorscale("plasma", m+1)
            star_colors = dict(zip(self.drivers.keys(), colors))
            for node in self.drivers.keys():
                name = (',').join(self.drivers[node])
                name = name.replace("'", "") 
                (x1, y1) = self.pos[node]
                (x2, y2)= self.pos[self.tree.nodes[node].parent.id]
                star = ((x1+x2)/2, (y1+y2)/2)
            
                fig.add_trace(go.Scatter(name=name, mode="markers", 
                                    x=[star[0]], y=[star[1]], 
                                    marker_symbol="star", 
                                    marker_color=star_colors[node],
                                    marker_size=28, hoverinfo='text', hovertext=name))


        fig.update_layout(title=self.title, annotations=self.annotations, font_size=14,
                          showlegend=True, legend_title_text='Mutation',
                          xaxis=dict(showline=False, zeroline=False, showgrid=False, showticklabels=False), 
                          yaxis=dict(showline=False, zeroline=False, showgrid=False, showticklabels=False))
        
        self.figs.append(fig)

    


    def plot_dist(self):
        fig = go.Figure()

        Xe, Ye = [], []

        for edge in self.edges:
            Xe += [self.pos[edge[0]][0][0], self.pos[edge[1]][0][0], None]
            Ye += [self.pos[edge[0]][0][1], self.pos[edge[1]][0][1], None]

        fig.add_trace(go.Scatter(x=Xe, y=Ye, mode='lines',
                    line=dict(color='#000000', width=3), hoverinfo=None, showlegend=False))
            
        for i in range(self.n):
            (x, y) = self.pos[i][0]
            color = self.colors[i]

            fig.add_trace(go.Scatter(name="clone "+str(i), x=[x], y=[y], mode='markers',
                             marker=dict(symbol='circle', size=30, color=color,
                                        line=dict(color='#000000', width=1))))

        fig.update_layout(title=self.title, 
                      xaxis_title=self.xtitle, 
                      yaxis_title=self.ytitle,
                      showlegend=True)
            
        self.figs.append(fig)



    def plot_freqs(self):
        for i in range(len(self.thepatient.samples)):
            sample = self.thepatient.samples[i]
            sampname = sample.name
            samptissue = sample.tissue
            samptime = sample.timePoint
    

            fig = go.Figure()

            Xe, Ye = [], []

            for edge in self.edges:
                Xe += [self.pos[edge[0]][i][0], self.pos[edge[1]][i][0], None]
                Ye += [self.pos[edge[0]][i][1], self.pos[edge[1]][i][1], None]

            fig.add_trace(go.Scatter(x=Xe, y=Ye, mode='lines',
                    line=dict(color='#000000', width=3), hoverinfo=None, showlegend=False))
                
            for j in range(self.n):
                (x, y) = self.pos[j][i]
                size = 30
                color = self.colors[j]
                if y < 0.03: # clone frequency threshold
                    size = 10
                    color = 'black'

                fig.add_trace(go.Scatter(name="clone "+str(j), x=[x], y=[y], mode='markers',
                             marker=dict(symbol='circle', size=size, color=color,
                                        line=dict(color='#000000', width=1))))
                    
            full_title = self.title+"<br><sup>Sample Name: {}</sup><br><sup>Tissue Type: {}</sup><br><sup>Anatomical Location: {}</sup>".format(sampname, samptime, samptissue)
                
            fig.update_layout(title=full_title, 
                      xaxis_title=self.xtitle, 
                      yaxis_title=self.ytitle,
                      showlegend=True, margin_t=150)
                
            self.figs.append(fig)







