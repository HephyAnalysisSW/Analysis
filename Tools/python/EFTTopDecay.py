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

        return {'0':Gamma_0, 'R':Gamma_LR+Gamma_LR_diff, 'L':Gamma_LR-Gamma_LR_diff }

if __name__=="__main__":
    eftTopWidth = EFTTopDecay()
    print eftTopWidth.widths(), "total", sum(eftTopWidth.widths().values())
