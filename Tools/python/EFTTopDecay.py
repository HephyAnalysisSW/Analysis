''' implement https://arxiv.org/pdf/hep-ph/0605190v2.pdf
    for reweighting the top decay in EFT
'''

from math import pi, sqrt

class EFTTopDecay:

    def __init__( self, MT = 172.5, MW = 80.387, MZ = 91.1876, MB = 4.18, Gf=0.000011663787):

        self.MT = MT
        self.MW = MW
        self.MB = MB
        self.Gf = Gf

        self.aEW = (Gf*MW**2*(1. - MW**2/MZ**2)*sqrt(2.))/pi
        self.ee  = 2.*sqrt(self.aEW)*sqrt(pi)
        self.sth = sqrt(1. - MW**2/MZ**2)
        self.g   = self.ee/self.sth
        self.q  = sqrt(MT**4 + MW**4 + MB**4 - 2*MT**2*MW**2 - 2*MT**2*MB**2 - 2*MW**2*MB**2)/(2.*MT)

        self.xW = MW/MT
        self.xb = MB/MT

    def widths( self, VLRe=1, VLIm=0, VRRe=0, VRIm=0, gLRe=0, gLIm=0, gRRe=0, gRIm=0):
        Gamma_0 = self.g**2*self.q/(32*pi)*( \
             self.MT**2/self.MW**2*( VLRe**2 + VLIm**2 + VRRe**2 + VRIm**2 )*( 1.-self.xW**2-2*self.xb**2-self.xW**2*self.xb**2+self.xb**4) -4.*self.xb*(VLRe*VRRe+VLIm*VRIm)\
           + (gLRe**2 + gLIm**2 + gRRe**2 + gRIm**2 )*(1-self.xW**2+self.xb**2)-4.*self.xb*(gLRe*gRRe+gLIm*gRIm)   \
           - 2*self.MT/self.MW*(VLRe*gRRe+VLIm*gRIm+VRRe*gLRe+VRIm*gLIm)*(1.-self.xW**2-self.xb**2)  \
           + 2*self.MT/self.MW*self.xb*(VLRe*gLRe+VLIm*gLIm+VRRe*gRRe+VRIm*gRIm)*(1.+self.xW**2-self.xb**2))
        Gamma_LR = self.g**2*self.q/(32*pi)*( \
             ( VLRe**2 + VLIm**2 + VRRe**2 + VRIm**2 )*( 1.-self.xW**2+self.xb**2) -4.*self.xb*(VLRe*VRRe+VLIm*VRIm)\
           + self.MT**2/self.MW**2*(gLRe**2 + gLIm**2 + gRRe**2 + gRIm**2 )*(1.-self.xW**2-2*self.xb**2-self.xW**2*self.xb**2+self.xb**4)-4.*self.xb*(gLRe*gRRe+gLIm*gRIm)   \
           - 2*self.MT/self.MW*(VLRe*gRRe+VLIm*gRIm+VRRe*gLRe+VRIm*gLIm)*(1.-self.xW**2-self.xb**2)  \
           + 2*self.MT/self.MW*self.xb*(VLRe*gLRe+VLIm*gLIm+VRRe*gRRe+VRIm*gRIm)*(1.+self.xW**2-self.xb**2))

        Gamma_LR_diff = self.g**2/(64*pi)*self.MT**3/self.MW**2*(-self.xW**2*( VLRe**2+VLIm**2-VRRe**2-VRIm**2) + (gLRe**2+gLIm**2-gRRe**2-gRIm**2)*(1-self.xb**2) \
           + 2*self.xW*(VLRe*gRRe+VLIm*gRIm-VRRe*gLRe-VRIm*gLIm)+2.*self.xW*self.xb*(VLRe*gLRe+VLIm*gLIm-VRRe*gRRe-VRIm*gRIm))*(1.-2*self.xW**2-2*self.xb**2+self.xW**4-2*self.xW**2*self.xb**2+self.xb**4)

        return {'Gamma_0':Gamma_0, 'Gamma_R':Gamma_LR+Gamma_LR_diff, 'Gamma_L':Gamma_LR-Gamma_LR_diff, 'Gamma':Gamma_0+2*Gamma_LR}

if __name__=="__main__":
    import Analysis.Tools.syncer as syncer
    import ROOT

    eftTopWidth = EFTTopDecay()
    widths      = eftTopWidth.widths()
    widths_BSM  = eftTopWidth.widths(gLRe=.2) 
    print("SM",  widths)
    print("BSM", widths_BSM)

    fR = ROOT.TF1("FR", "{FR}*(3/8.)*(1+cos(x))**2".format(FR=widths['Gamma_R']/widths['Gamma']),0,pi)
    fL = ROOT.TF1("FL", "{FL}*(3/8.)*(1-cos(x))**2".format(FL=widths['Gamma_L']/widths['Gamma']),0,pi)
    f0 = ROOT.TF1("F0", "{F0}*(3/4.)*sin(x)**2".format(F0=widths['Gamma_0']/widths['Gamma']),0,pi)
    tot= ROOT.TF1("total", "{F0}*(3/4.)*sin(x)**2+{FL}*(3/8.)*(1-cos(x))**2+{FR}*(3/8.)*(1+cos(x))**2".format(F0=widths['Gamma_0']/widths['Gamma'],FL=widths['Gamma_L']/widths['Gamma'],FR=widths['Gamma_R']/widths['Gamma']),0,pi)


    tot.SetLineColor(ROOT.kBlack)
    fR.SetLineColor(ROOT.kBlue)
    fL.SetLineColor(ROOT.kGreen)
    f0.SetLineColor(ROOT.kRed)
    tot.SetLineStyle(2)
    fR .SetLineStyle(2)
    fL .SetLineStyle(2)
    f0 .SetLineStyle(2)

    fR_BSM = ROOT.TF1("FR", "{FR}*(3/8.)*(1+cos(x))**2".format(FR=widths_BSM['Gamma_R']/widths_BSM['Gamma']),0,pi)
    fL_BSM = ROOT.TF1("FL", "{FL}*(3/8.)*(1-cos(x))**2".format(FL=widths_BSM['Gamma_L']/widths_BSM['Gamma']),0,pi)
    f0_BSM = ROOT.TF1("F0", "{F0}*(3/4.)*sin(x)**2".format(F0=widths_BSM['Gamma_0']/widths_BSM['Gamma']),0,pi)
    tot_BSM= ROOT.TF1("total", "{F0}*(3/4.)*sin(x)**2+{FL}*(3/8.)*(1-cos(x))**2+{FR}*(3/8.)*(1+cos(x))**2".format(F0=widths_BSM['Gamma_0']/widths_BSM['Gamma'],FL=widths_BSM['Gamma_L']/widths_BSM['Gamma'],FR=widths_BSM['Gamma_R']/widths_BSM['Gamma']),0,pi)


    tot_BSM.SetLineColor(ROOT.kBlack)
    fR_BSM.SetLineColor(ROOT.kBlue)
    fL_BSM.SetLineColor(ROOT.kGreen)
    f0_BSM.SetLineColor(ROOT.kRed)

    c1 = ROOT.TCanvas()
    tot.Draw()
    tot.GetYaxis().SetRangeUser(0,1)
    fR.Draw("same")
    fL.Draw("same")
    f0.Draw("same")

    tot_BSM.Draw("same")
    fR_BSM.Draw("same")
    fL_BSM.Draw("same")
    f0_BSM.Draw("same")

    c1.Print('./hel.png')
    
