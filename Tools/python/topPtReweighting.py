''' Functions for top pt reweighting
'''
from math import exp, sqrt

#https://indico.cern.ch/event/463433/contribution/5/attachments/1193635/1733505/MSavitskyi_ToppTRew13TeV_TopModGen.pdf
#updated a and b parameters for 13 TeV: https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#pt_top_Reweighting

def getUnscaledTopPairPtReweightungFunction():
    a = 0.0615
    b = -0.0005

    def f(topPts):
        assert len(topPts)<=2, "Found too many top pts: %i"%len(topPts)
        res=1.
        for pt in topPts:
            res*=exp(a+b*pt)
        return sqrt(res)
    return f

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
def getTopPtDrawString():
    a = 0.0615
    b = -0.0005
    return "sqrt(exp(Sum$(("+str(a)+"+("+str(b)+")*GenPart_pt)*(abs(GenPart_pdgId)==6&&GenPart_status==62&&GenPart_pt<=400))))"

def getTopPtsForReweighting(r):
    res=[]
    for i in range(r.nGenPart):
        try:
            if abs(r.GenPart_pdgId[i])==6 and r.GenPart_status[i]==62 and r.GenPart_pt[i]<400:
                res.append(r.GenPart_pt[i])
        except IndexError: #rare case with too many particles
            pass
    return res
