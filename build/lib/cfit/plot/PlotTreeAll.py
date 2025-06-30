import os

from cfit.plot.PlotTree import PlotTree

class PlotTreeAll():

    def __init__(self, patient, drivers=None, outdir=None, tree_format=None, tid=0, save_html=True, show=False):
        """
        Params
        ----------

        patient : PatientLine object
            patient whose trees will be graphed
        
        drivers : list
            list of driver genes
            if specified, will mark branches of clonetree where driver gene mutations are acquired
        
        tree_format : str
            method used to resonstruct phylogeny
            options: "phylowgs", "pairtree", "mclust"

        outdir : str
            path to output directory where html file with graphs will be saved
            only needed if save_html=True

        tid : int
            tree index

        save_html: bool
            indicates whether or not to save generated html

        show : bool
            indicates whether or not to immediately show generated graphs (each opened in new browser window)

        """

        self.figures = []

        self.__init_plottrees(patient, drivers, tid)

        if show:
            self.show_figs()
        
        if save_html:
            if outdir is None:
                raise Exception("with save_html=True, please provide output directory path")
            t = ""
            if tree_format is not None:
                t = tree_format
            self.save_as_html(outdir, patient.name, t)

    
    def __init_plottrees(self, patient, drivers, tid):
        for type in ["topology"]:
#        for type in ["topology", "dist vs incmuts", "dist vs excmuts", "dist vs neos", "muts vs freq", "dist vs incfreq", "dist vs excfreq"]:
            plottree = PlotTree(patient, type, drivers, tid)
            plottree.plot()
            self.figures += plottree.figs

    def show_figs(self):
        for figure in self.figures:
            figure.show()
        
    def save_as_html(self, outdir, patname, tree):
        filename = patname+'_'+tree+'_report.html'
        if tree == "":
            filename = patname+'_report.html'
        with open(os.path.join(outdir, filename), 'w') as f:
            f.write("<html> <body> <h1> {} Tree Graphs Report for \
           <font color = #000000>{}</font></h1>\n \
           </body></html>".format(tree, patname))
            for figure in self.figures:
                f.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))
    
