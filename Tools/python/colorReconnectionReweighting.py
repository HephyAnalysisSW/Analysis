''' Functions for top pt reweighting
'''
from math import exp, sqrt

def getCRWeight(njets):
    if njets < 3:
        return 1.
    elif njets == 3:
        return (1-0.0122834751382)
    elif njets == 4:
        return 1+0.00650110282004
    elif njets >= 5:
        return 1+0.0272257402539
    
def getCRDrawString():
    nJet = "Sum$(Jet_pt>30&&abs(Jet_eta)<2.4&&Jet_id>0)"
    return "(1*({nJet}<3) + (1-0.0122834751382)*({nJet}==3) + (1+0.00650110282004)*({nJet}==4) + (1+0.0272257402539)*({nJet}>=5))".format( nJet = nJet )
