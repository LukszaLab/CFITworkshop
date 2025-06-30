from cfit.fitness.neo_quality import CARNeoantigenQuality
import numpy as np

class ModelNeoantigenQuality(CARNeoantigenQuality):


    ##############PARAMETERS###################
    kon_diff = 1e-4  ##<-in nanomolar (1e4*1e-9)
    kstep = 1/(1.5*60)
    kload = 1/(30*60)
    rho0 = 1e3
    xm = 1e5/6 #1e5 divided by 6 different alleles
    kout = 1e-5
    kin = 1e-5
    #################END####################


    def __init__(self, iedbfasta="", alndir="", matrix="blosum62",
                 ep_dist_model_name="all_tcr_all_combos_model"):
        '''

        :param iedbfasta: str
            path to the IEDB fasta file

        :param alndir: str
            path to the alignment directory, neoantigens/*

        :param matrix: str

        :param ep_dist_model_name: str

        '''
        CARNeoantigenQuality.__init__(self, iedbfasta=iedbfasta, alndir=alndir, matrix=matrix, ep_dist_model_name=ep_dist_model_name)



    def y_pm(self, kd, kon_diff, xm, kstep, kload, rho0, kin, kout):
        '''
        :param kd: float
            xxx
        :param kon_diff: float
            xxx
        :param xm: float
            xxx
        :param kstep: float
            xxx
        :param kload: float
            xxx
        :param rho0: float
            xxx
        :param kin: float
            xxx
        :param kout: float
            xxx
        :return float
        '''

        koff = kd*kon_diff ##From K_D(nM) to k_off (1/s)
        V_er = 2.5*((4/3)*np.pi*(5e-5)**3)/10 #This compute the volume of the ER (radius in dm)
        p1 = 1/(1 + koff/kstep) #Step 1: proofreading intra ER
        p2 = 1/(kload + koff) #Step 2: proofreading extra ER
        M = (xm/(V_er*6e23))*1e9 #Concentration of MHC inside the ER (now it's in nM)
        kon_eff = self.bivalue_kon_eff(kin) #Effective association rate intra ER
        p3 = (rho0*kin*xm*kon_eff)/(kout + (M*kon_eff)/(1+koff/kstep)) #contribution from intra ER peptide density
        return p1*p2*p3

    #This function gives value of the effective association rate based on intra-cellular dynamics. Only two regimes are considered
    def bivalue_kon_eff(self, kin):
        if kin ==1:
            kon_eff = 1e-9
        if kin == 1e-5:
            kon_eff = 1e-6
        return kon_eff


    def compute_neoantigen_sample_quality(self, neo, sample, a=26, k=1, w=0.5, kdthr=None,
                                          include_R=True, just_R=False):
        '''

        Compute quality of a neoantigen in a sample.

        :param neo: cfit.tree.mutation.Neoantigen

        :param sample: cfit.patient.Sample

        :param a: float
            shift parameter of the R term

        :param k: float
            slope component of the R term

        :param w: float
            relative weight of logC term [ w*log C + (1-w)*log A ]

        :param kdthr: float
            hard threshold on the dissociation constant

        :param include_R: bool
            whether to include R: if false (*)*R model is used

        :param just_R: bool
            whether to include only R: if true R model is used

        :return: float
            neoantigen quality

        '''

        #        if len(neo.mtPeptide) != 9:
        #            return 0
        neo.qattributes.GbSum = -a
        neo_concentration = sample.mutations[neo.mid].expression

        D = self.epitope_dist(neo.wtPeptide, neo.mtPeptide)
        A = self.compute_A(neo, KDnormalize=1, neo_concetration=neo_concentration)

        R = 1.
        if include_R:
            R = self.compute_recognition_probability(neo, a, k)

        neo.qattributes.D = D
        neo.qattributes.A = A
        neo.qattributes.R = R
        neo.qattributes.neotype_prob = self.compute_neotype_prob(neo)
        neo.qattributes.mut_pres_prob = self.compute_mut_presentation(neo)


        if just_R:
            quality = R
        else:
            quality = ((1 - w) * np.log(A) + w * D) * R

        if kdthr is not None:
            if neo.kD > kdthr:
                quality = 0.0
        return quality
