# MetFilter Analysis Recommendations according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
# Flag_BadChargedCandidateFilter is not recommended anymore, under review
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2

def getFilterCut( year, isData=False, ignoreJSON=False, isFastSim=False, skipBadChargedCandidate=True, skipBadPFMuon=False, skipVertexFilter=False, skipWeight=False):

    if year in [2016, "Run2016"]:
        filters             = [ "Flag_goodVertices" ]                         # primary vertex filter
        if not isFastSim:
            filters        += [ "Flag_globalSuperTightHalo2016Filter" ]       # beam halo filter
        filters            += [ "Flag_HBHENoiseFilter" ]                      # HBHE noise filter
        filters            += [ "Flag_HBHENoiseIsoFilter" ]                   # HBHEiso noise filter
        filters            += [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ]   # ECAL TP filter
        if not skipBadPFMuon:
            filters        += [ "Flag_BadPFMuonFilter" ]                      # Bad PF Muon Filter
        if not skipBadChargedCandidate: #recommended to skip for now!!
            filters        += [ "Flag_BadChargedCandidateFilter" ]            # Bad Charged Hadron Filter
        if isData:
            filters        += [ "Flag_eeBadScFilter" ]                        # ee badSC noise filter (data only)

    elif year in [2017, "Run2017"]:
        filters             = [ "Flag_goodVertices" ]                         # primary vertex filter
        if not isFastSim:
            filters        += [ "Flag_globalSuperTightHalo2016Filter" ]       # beam halo filter
        filters            += [ "Flag_HBHENoiseFilter" ]                      # HBHE noise filter
        filters            += [ "Flag_HBHENoiseIsoFilter" ]                   # HBHEiso noise filter
        filters            += [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ]   # ECAL TP filter
#        filters            += [ "Flag_ecalBadCalibReducedMINIAODFilter" ]     # ECAL bad calibration filter update (needs to be re-run on miniAOD)
        filters            += [ "Flag_ecalBadCalibFilter" ]                   # current replacement for Flag_ecalBadCalibReducedMINIAODFilter -> change to Flag_ecalBadCalibFilterv2
        if not skipBadPFMuon:
            filters        += [ "Flag_BadPFMuonFilter" ]                      # Bad PF Muon Filter
        if not skipBadChargedCandidate: #recommended to skip for now!!
            filters        += [ "Flag_BadChargedCandidateFilter" ]            # Bad Charged Hadron Filter
        if isData:
            filters        += [ "Flag_eeBadScFilter" ]                        # ee badSC noise filter (data only)
        if isFastSim:
            filters        += ["Flag_ecalBadCalibFilter"]
        else:
            filters        += ["Flag_ecalBadCalibFilterV2"]

    elif year in [2018, "Run2018"]:
        filters             = [ "Flag_goodVertices" ]                         # primary vertex filter
        if not isFastSim:
            filters        += [ "Flag_globalSuperTightHalo2016Filter" ]       # beam halo filter
        filters            += [ "Flag_HBHENoiseFilter" ]                      # HBHE noise filter
        filters            += [ "Flag_HBHENoiseIsoFilter" ]                   # HBHEiso noise filter
        filters            += [ "Flag_EcalDeadCellTriggerPrimitiveFilter" ]   # ECAL TP filter
#        filters            += [ "Flag_ecalBadCalibReducedMINIAODFilter" ]     # ECAL bad calibration filter update (needs to be re-run on miniAOD)
        filters            += [ "Flag_ecalBadCalibFilter" ]                   # current replacement for Flag_ecalBadCalibReducedMINIAODFilter -> change to Flag_ecalBadCalibFilterv2
        if not skipBadPFMuon:
            filters        += [ "Flag_BadPFMuonFilter" ]                      # Bad PF Muon Filter
        if not skipBadChargedCandidate: #recommended to skip for now!!
            filters        += [ "Flag_BadChargedCandidateFilter" ]            # Bad Charged Hadron Filter
        if isData:
            filters        += [ "Flag_eeBadScFilter" ]                        # ee badSC noise filter (data only)
        if isFastSim:
            filters        += ["Flag_ecalBadCalibFilter"]
        else:
            filters        += ["Flag_ecalBadCalibFilterV2"]

    elif year=="RunII":
        return "((year==2016&&{filter_2016})||(year==2017&&{filter_2017})||(year==2018&&{filter_2018}))".format( 
            filter_2016 = getFilterCut( 2016, isData=isData, ignoreJSON=ignoreJSON, isFastSim=isFastSim, skipBadChargedCandidate=skipBadChargedCandidate, skipBadPFMuon=skipBadPFMuon, skipVertexFilter=skipVertexFilter),
            filter_2017 = getFilterCut( 2017, isData=isData, ignoreJSON=ignoreJSON, isFastSim=isFastSim, skipBadChargedCandidate=skipBadChargedCandidate, skipBadPFMuon=skipBadPFMuon, skipVertexFilter=skipVertexFilter),
            filter_2018 = getFilterCut( 2018, isData=isData, ignoreJSON=ignoreJSON, isFastSim=isFastSim, skipBadChargedCandidate=skipBadChargedCandidate, skipBadPFMuon=skipBadPFMuon, skipVertexFilter=skipVertexFilter),
          )
    else:
        raise NotImplementedError( "No MET filter found for year %i" %year )

    if not skipVertexFilter:
        filters            += [ "PV_ndof>4", "sqrt(PV_x*PV_x+PV_y*PV_y)<=2", "abs(PV_z)<=24" ]

    if isData:
        if not skipWeight:
            filters        += [ "weight>0" ]
        if not ignoreJSON:
            filters        += [ "jsonPassed>0" ]

    return "&&".join(filters)

