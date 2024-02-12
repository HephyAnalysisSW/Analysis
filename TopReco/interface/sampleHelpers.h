#ifndef sampleHelpers_h
#define sampleHelpers_h

#include <vector>
#include <string>

#include <TString.h>

/// Namespace to treat different analysis eras as enum types
namespace Era{

    /// All analysis eras as needed
    enum Era {

      run1_8tev,
      run2_13tev_50ns,
      run2_13tev_25ns,
      run2_13tev_2015_25ns,
      run2_13tev_2016_25ns,
      run2_13tev_25ns_74X,
      run2_13tev_2016,
      run2_13tev_2017,
      run2_13tev_2018,
      run2_13tev_full,

      undefined
    };

    /// Convert an era from string to enum
    Era convert(const TString& era);

    /// Convert an era from enum to string
    TString convert(const Era& era);

    /// Return energy for given era in TeV
    float energyInTev(const Era era);
}



namespace lepSF{

  enum readoutFormat{undefined, pt_vs_eta, eta_vs_pt, pt_vs_absEta, absEta_vs_pt,
                     pt_vs_etaSc, etaSc_vs_pt, pt_vs_absEtaSc, absEtaSc_vs_pt};

  readoutFormat convertReadoutFormat(const std::string& readoutFormat_);

  class SFHisto{
  public:
    SFHisto();
    SFHisto(const std::string& fileName_, const std::string& histName_,
            const std::string& readout_, const std::string& histNameStat_,
            const std::string& histNameSys_);
    std::string getFile()const;
    std::string getHisto()const;
    std::string getStatHisto()const;
    std::string getSysHisto()const;
    readoutFormat getFormat()const;
  private:
    std::string fileName;
    std::string histName;
    std::string histName_Stat;
    std::string histName_Sys;
    readoutFormat readout;
  };
}



namespace triggerSF{

  enum readoutFormat{undefined, el0_pt_vs_el1_pt, el0_abseta_vs_el1_abseta, el0_pt_vs_el0_abseta, el1_pt_vs_el1_abseta,
                     mu0_pt_vs_mu0_abseta, mu1_pt_vs_mu1_abseta, mu0_pt_vs_mu1_pt, mu0_abseta_vs_mu1_abseta,
                     mu0_pt_vs_el0_pt, mu0_abseta_vs_el0_abseta, el0_pt_vs_mu0_pt, el0_abseta_vs_mu0_abseta};

  readoutFormat convertReadoutFormat(const std::string& readoutFormat_);

  class SFHisto{
  public:
    SFHisto();
    SFHisto(const std::string& fileName_, const std::string& histName_,
            const std::string& readout_, const std::string& histNameStat_,
            const std::string& histNameSys_);
    std::string getFile()const;
    std::string getHisto()const;
    std::string getStatHisto()const;
    std::string getSysHisto()const;
    readoutFormat getFormat()const;
  private:
    std::string fileName;
    std::string histName;
    std::string histName_Stat;
    std::string histName_Sys;
    readoutFormat readout;
  };
}



/// Namespace to define enums needed for jet pileup ID
namespace JetPileupId{

    /// Enum for the working points
    enum WorkingPoint{
        none, L, M, T, undefinedWP
    };



    /// Convert a WorkingPoint from string to enum
    WorkingPoint convertWorkingPoint(const std::string& wp);

    /// Convert a WorkingPoint from enum to string
    std::string convertWorkingPoint(const WorkingPoint& wp);

    /// Cut value for corresponding working point, to be used with >= operator
    int cutValue(const WorkingPoint& wp);
}








/// Namespace to define enums needed for b-tag related stuff
/// Also needed for b-tag scale factors, since ROOT dictionary does not want definitions done in ZTopUtils/
namespace Btag{

    // FIXME: csv, csvv2_50ns and perhaps others are obsolete types, should be removed
    /// Enum for the implemented b-tagging algorithms
    enum Algorithm{

        csv,
        cmvav2,
        cmvav2_2016,
        csvv2,
        csvv2_2016,
        csvv2_74X,
        csvv2_50ns,
        csvv2HIP,

        DeepCSV_2016_94X,
        DeepJet_2016_94X,

        CSVv2_2017_94X,
        DeepCSV_2017_94X,
        DeepJet_2017_94X,

        DeepCSV_2018_102X,
        DeepJet_2018_102X,

        undefinedAlgorithm
    };

    /// Enum for the working points
    enum WorkingPoint{
        L, M, T, undefinedWP
    };

    /// Enum for the implemented modes of btag corrections
    enum CorrectionMode{
        noCorrection,               // Do not apply any corrections, i.e. scale factors event SF=1
        greaterEqualOneTagReweight, // Correct selection efficiency for given working point via event SF for >=N or ==N b-tag selection with N any integer >=1
        randomNumberRetag,          // Random-number based tag flipping for b-/c-/l-jets to correct for selection efficiency
        discriminatorReweight,      // Reweight with event-wise SF to describe b-tag discriminator distribution
        undefinedCorrectionMode     // Undefined
    };



    /// Convert an Algorithm from string to enum
    Algorithm convertAlgorithm(const std::string& algo);

    /// Convert an Algorithm from enum to string
    std::string convertAlgorithm(const Algorithm& algo);

    /// Convert a WorkingPoint from string to enum
    WorkingPoint convertWorkingPoint(const std::string& wp);

    /// Convert a WorkingPoint from enum to string
    std::string convertWorkingPoint(const WorkingPoint& wp);

    /// Convert a CorrectionMode from string to enum
    CorrectionMode convertCorrectionMode(const std::string& mode);

    /// Convert a CorrectionMode from enum to string
    std::string convertCorrectionMode(const CorrectionMode& mode);
}







/// Namespace to define enums needed for MET related stuff
namespace Met{

    // FIXME: mva and pfNoHf are obsolete types, should be removed
    /// Enum for the implemented MET algorithms
    enum Algorithm{
        pf,
        pfT1,
        pfT01XY,
        mva,
        pfNoHf,
        puppi,
        undefinedAlgorithm
    };

    /// Convert an Algorithm from string to enum
    Algorithm convertAlgorithm(const std::string& algo);

    /// Convert an Algorithm from enum to string
    std::string convertAlgorithm(const Algorithm& algo);
}







/// Namespace to treat systematics as enum types
namespace Systematic{

    /// All systematic types as needed in any part of the framework
    enum Type{
        nominal,            // nominal, i.e. no systematic variation applied
        mH110,              // Higgs mass of 110 GeV
        mH115,              // Higgs mass of 115 GeV
        mH120,              // Higgs mass of 120 GeV
        mH1225,             // Higgs mass of 122.5 GeV
        mH1275,             // Higgs mass of 127.5 GeV
        mH130,              // Higgs mass of 130 GeV
        mH135,              // Higgs mass of 135 GeV
        mH140,              // Higgs mass of 140 GeV
        mTop166,            // top quark mass of 166 GeV
        mTop169,            // top quark mass of 169 GeV
        mTop171,            // top quark mass of 171 GeV
        mTop173,            // top quark mass of 173 GeV
        mTop175,            // top quark mass of 175 GeV
        mTop178,            // top quark mass of 178 GeV
        lept,               // scale lepton ID/ISO data-to-MC scale factors
        ele,                // scale electron ID/ISO data-to-MC scale factors
        eleID,              // scale electron ID/ISO data-to-MC scale factors
        eleIDStat,          // scale electron ID/ISO data-to-MC scale factors
        eleIDSyst,          // scale electron ID/ISO data-to-MC scale factors
        eleReco,            // scale electron ID/ISO data-to-MC scale factors
        eleRecoStat,        // scale electron ID/ISO data-to-MC scale factors
        eleRecoSyst,        // scale electron ID/ISO data-to-MC scale factors
        muon,               // scale muon ID/ISO data-to-MC scale factors
        muonID,             // scale muon ID/ISO data-to-MC scale factors
        muonIDStat,         // scale muon ID/ISO data-to-MC scale factors
        muonIDSyst,         // scale muon ID/ISO data-to-MC scale factors
        muonIso,            // scale muon ID/ISO data-to-MC scale factors
        muonIsoStat,        // scale muon ID/ISO data-to-MC scale factors
        muonIsoSyst,        // scale muon ID/ISO data-to-MC scale factors
        eleScaleSyst,       // electron scale correction
        eleScaleGain,       // electron scale correction
        eleScaleStat,       // electron scale correction
        eleScaleEt,         // electron scale correction
        eleSmearingPhi,     // electron smearing correction
        eleSmearingRho,     // electron smearing correction
        eleScaleSmearing,   // electron scale+smearing enevlope
        muonScaleStat,      // muon scale correction
        muonScaleZpt,       // muon scale correction
        muonScaleEwk,       // muon scale correction
        muonScaleDeltaM,    // muon scale correction
        muonScaleEwk2,      // muon scale correction
        muonScale,          // muon scale correction envelope
        trig,               // scale trigger data-to-MC scale factors
        trigEta,            // scale trigger data-to-MC scale factors wrt eta (barrel-or-endcap) in antagonistic way
        pu,                 // scale pileup data-to-MC scale factors
        dy,                 // uncertainty on the Drell-Yan same-flavour background
        bg,                 // general background uncertainty
        dynorm,             // uncertainty on the Drell-Yan background estimation normalization
        kin,                // scale kinematic reconstruction scale factors
        jetPileupID,        // jet pileup-ID Data/MC scale factors
        btag,               // scale b-tagging data-to-MC scale factors of the b-/c-jets
        btagPt,             // median method: scale b-tagging data-to-MC scale factors of the b-/c-jets below/above median pt down/up or up/down
        btagEta,            // median method: scale b-tagging data-to-MC scale factors of the b-/c-jets below/above median eta down/up or up/down
        btagLjet,           // scale b-tagging data-to-MC scale factors of the l-jets
        btagLjetPt,         // median method: scale b-tagging data-to-MC scale factors of the l-jets below/above median pt down/up or up/down
        btagLjetEta,        // median method: scale b-tagging data-to-MC scale factors of the l-jets below/above median eta down/up or up/down
        btagBeff,           // scale the b-tagging efficiencies as estimated from MC for b-jets for stat. uncertainty (not applied anywhere, should it be removed?)
        btagCeff,           // scale the b-tagging efficiencies as estimated from MC for c-jets for stat. uncertainty (not applied anywhere, should it be removed?)
        btagLeff,           // scale the b-tagging efficiencies as estimated from MC for l-jets for stat. uncertainty (not applied anywhere, should it be removed?)
        btagDiscrBpurity,   // for b-tag discriminator reweighting: purity of the HF sample used for the LF SF determination
        btagDiscrLpurity,   // for b-tag discriminator reweighting: purity of the LF sample used for the HF SF determination
        btagDiscrBstat1,    // for b-tag discriminator reweighting: scale part 1 of the statistical uncertainty for b-jets
        btagDiscrBstat2,    // for b-tag discriminator reweighting: scale part 2 of the statistical uncertainty for b-jets
        btagDiscrLstat1,    // for b-tag discriminator reweighting: scale part 1 of the statistical uncertainty for l-jets
        btagDiscrLstat2,    // for b-tag discriminator reweighting: scale part 2 of the statistical uncertainty for l-jets
        btagDiscrCerr1,     // for b-tag discriminator reweighting: scale part 1 of the total uncertainty for c-jets
        btagDiscrCerr2,     // for b-tag discriminator reweighting: scale part 2 of the total uncertainty for c-jets
        subjetbtagDiscrBpurity, // subjet b-tag discriminator reweighting: purity of the HF sample used for the LF SF determination
        subjetbtagDiscrLpurity, // subjet b-tag discriminator reweighting: purity of the LF sample used for the HF SF determination
        subjetbtagDiscrBstat1,  // subjet b-tag discriminator reweighting: scale part 1 of the statistical uncertainty for b-jets
        subjetbtagDiscrBstat2,  // subjet b-tag discriminator reweighting: scale part 2 of the statistical uncertainty for b-jets
        subjetbtagDiscrLstat1,  // subjet b-tag discriminator reweighting: scale part 1 of the statistical uncertainty for l-jets
        subjetbtagDiscrLstat2,  // subjet b-tag discriminator reweighting: scale part 2 of the statistical uncertainty for l-jets
        subjetbtagDiscrCerr1,   // subjet b-tag discriminator reweighting: scale part 1 of the total uncertainty for c-jets
        subjetbtagDiscrCerr2,   // subjet b-tag discriminator reweighting: scale part 2 of the total uncertainty for c-jets
        subjetbtagHeavyFlavor,  // subjet b-tagging efficiency
        subjetbtagLightFlavor,  // subjet b-tagging mis-tag rate
        jerEta0,    // scale jet energy resolution scale factors (pt/eta bin0)
        jerEta1,    // scale jet energy resolution scale factors (pt/eta bin1)
        jerEta2Pt0, // scale jet energy resolution scale factors (pt/eta bin2)
        jerEta2Pt1, // scale jet energy resolution scale factors (pt/eta bin3)
        jerEta3Pt0, // scale jet energy resolution scale factors (pt/eta bin4)
        jerEta3Pt1, // scale jet energy resolution scale factors (pt/eta bin5)
        jer,        // scale jet energy resolution scale factors
        jes,        // scale jet energy scale scale factors
        jesAbsoluteStat, //0
        jesAbsoluteScale,//1
        jesAbsoluteFlavMap,//2
        jesAbsoluteMPFBias,//3
        jesFragmentation,   // "HighPtExtra,    //4 -----
        jesSinglePionECAL, //5
        jesSinglePionHCAL, //6
        jesFlavorQCD,        //7
        jesTimePtEta,               // "Time" 8 -----
        jesRelativeJEREC1,    //9
        jesRelativeJEREC2,    //10
        jesRelativeJERHF,   //11
        jesRelativePtBB,       //12
        jesRelativePtEC1,    //13
        jesRelativePtEC2,    //14
        jesRelativePtHF,    //15
        jesRelativeBal, //new
        jesRelativeFSR,       //16
        jesRelativeSample,
        jesRelativeStatFSR,   //new
        jesRelativeStatEC,    //2,    //17 -----
        jesRelativeStatHF,   //18
        jesPileUpDataMC,       //19
        jesPileUpPtRef,
        jesPileUpPtBB,        //20
        jesPileUpPtEC1,        //21-----PileUpPtRef]
        jesPileUpPtEC2,
        jesPileUpPtHF,       //22
            //"PileUpBias,       //23-----
        jesPileUpMuZero,
        jesPileUpEnvelope,
        jesSubTotalPileUp,   //24
        jesSubTotalRelative,   //25
        jesSubTotalPt,       //26
        jesSubTotalScale,       //26
        jesSubTotalMC,       //27
        jesSubTotalAbsolute,       //26
        jesTotalNoFlavor,   //29
            //"TotalNoTime"       // ignoring for the moment these subtotals
            //"TotalNoFlavorNoTime"
            //"Time A-D
        jesFlavorZJet,       //30
        jesFlavorPhotonJet,   //31
        jesFlavorPureGluon,   //32
        jesFlavorPureQuark,   //33
        jesFlavorPureCharm,   //34
        jesFlavorPureBottom,//35
        jesCorrelationGroupMPFInSitu,//36
        jesCorrelationGroupIntercalibration,//37
        jesCorrelationGroupbJES,//38
        jesCorrelationGroupFlavor,//39
        jesCorrelationGroupUncorrelated,//40
        jesHttEta0to5, jesHttEta0to3, jesHttEta3to5, jesHttEC2,
        jesAbsolute, jesAbsoluteYear, jesBBEC1, jesBBEC1Year, jesEC2, jesEC2Year, jesHF, jesHFYear,
        jesUserDefinedHEM1516,
        jesFlavorRealistic,

        frac_tthf,          // correction factor for the fraction of tt+HF events from the template fit
        frac_ttother,       // correction factor for the fraction of tt+Other events from the template fit
        lumi,               // luminosity uncertainty
        xsec_ttbb,          // cross-section uncertainty of ttbb process
        xsec_ttb,           // cross-section uncertainty of ttb process
        xsec_tt2b,          // cross-section uncertainty of tt2b process
        xsec_ttcc,          // cross-section uncertainty of ttcc process
        xsec_ttother,       // cross-section uncertainty of tt+light jets process
        xsec_ttZ,           // cross-section uncertainty of ttZ process
        xsec_ttW,           // cross-section uncertainty of ttW process
        xsec_ttH,           // cross-section uncertainty of ttH process
        xsec_ttG,           // cross-section uncertainty of ttG process
        xsec_ttV,           // cross-section uncertainty of ttV process
        xsec_tt,            // cross-section uncertainty of ttbar process
        xsec_t,             // cross-section uncertainty of single top process
        xsec_ttDM,          // cross-section uncertainty of tt+DM process
        topPtTheory,        // scale top pt as predicted in theoretical ttbar differential cross-section calculations
        topPtFitparam0,     // scale top pt as estimated in ttbar differential cross-section measurements, uncertainty via decorrelated fit parameter p0prime
        topPtFitparam1,     // scale top pt as estimated in ttbar differential cross-section measurements, uncertainty via decorrelated fit parameter p1prime
        topPt,              // scale top pt as estimated in ttbar differential cross-section measurements, uncertainty via switching off and applying twice
        mass,               // variations of masses used in process generation (here top quark mass)
        match,              // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter
        match_ttbb,         // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_tt2b,         // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_ttb,          // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_ttcc,         // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_ttother,      // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        erdon,
        erdonretune,
        gluonmovetune,
        meScale,            // Q2 scale uncertainty in process generation on Matrix Element only
        meScale_ttbb,       // Q2 scale uncertainty in process generation on Matrix Element only (ttbb process)
        meScale_ttb,        // Q2 scale uncertainty in process generation on Matrix Element only (ttb process)
        meScale_tt2b,       // Q2 scale uncertainty in process generation on Matrix Element only(tt2b process)
        meScale_ttcc,       // Q2 scale uncertainty in process generation on Matrix Element only (ttcc process)
        meScale_ttother,    // Q2 scale uncertainty in process generation on Matrix Element only (tt+light jets process)
        meFacScale,            // Q2 factorization scale uncertainty in process generation on Matrix Element only
        meFacScale_ttbb,       // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbb process)
        meFacScale_ttb,        // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttb process)
        meFacScale_tt2b,       // Q2 factorization scale uncertainty in process generation on Matrix Element only(tt2b process)
        meFacScale_ttcc,       // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttcc process)
        meFacScale_ttother,    // Q2 factorization scale uncertainty in process generation on Matrix Element only (tt+light jets process)
        meFacScale_tt,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar process)
        meFacScale_z,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (z+jets process)
        meFacScale_w,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (w+jets process)
        meFacScale_st,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (single top process)
        meFacScale_ttv,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+X process)
        meFacScale_ttz,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+Z process)
        meFacScale_ttw,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+W process)
        meFacScale_ttg,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+Gamma process)
        meFacScale_vv,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (diboson process)
        meFacScale_ttdm,   // Q2 factorization scale uncertainty in process generation on Matrix Element only (tt+DM process)
        meFacScale_htott_res,   // Q2 factorization scale uncertainty in process generation on Matrix Element only (H/A->tt resonant process)
        meFacScale_htott_int,   // Q2 factorization scale uncertainty in process generation on Matrix Element only (H/A->tt interference process)
        meRenScale,            // Q2 renormalization scale uncertainty in process generation on Matrix Element only
        meRenScale_ttbb,       // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbb process)
        meRenScale_ttb,        // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttb process)
        meRenScale_tt2b,       // Q2 renormalization scale uncertainty in process generation on Matrix Element only(tt2b process)
        meRenScale_ttcc,       // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttcc process)
        meRenScale_ttother,    // Q2 renormalization scale uncertainty in process generation on Matrix Element only (tt+light jets process)
        meRenScale_tt,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar process)
        meRenScale_z,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (z+jets process)
        meRenScale_w,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (w+jets process)
        meRenScale_st,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (single top process)
        meRenScale_ttv,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+X process)
        meRenScale_ttz,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+Z process)
        meRenScale_ttw,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+W process)
        meRenScale_ttg,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+Gamma process)
        meRenScale_vv,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (diboson process)
        meRenScale_ttdm,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (tt+DM process)
        meRenScale_htott_res,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (H/A->tt resonant process)
        meRenScale_htott_int,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (H/A->tt interference process)

        psScale,            // Q2 scale uncertainty in process generation on Parton Shower only
        psScale_ttbb,       // Q2 scale uncertainty in process generation on Parton Shower only (ttbb process)
        psScale_ttb,        // Q2 scale uncertainty in process generation on Parton Shower only (ttb process)
        psScale_tt2b,       // Q2 scale uncertainty in process generation on Parton Shower only (tt2b process)
        psScale_ttcc,       // Q2 scale uncertainty in process generation on Parton Shower only (ttcc process)
        psScale_ttother,    // Q2 scale uncertainty in process generation on Parton Shower only (tt+light jets process)
        psISRScale,         // alpha_s^ISR scale uncertainty in process generation on Parton Shower only
        psISRScale_ttbb,    // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_tt2b,    // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_ttb,     // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_ttcc,    // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_ttother, // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale,         // alpha_s^FSR scale uncertainty in process generation on Parton Shower only
        psFSRScale_ttbb,    // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_tt2b,    // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_ttb,     // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_ttcc,    // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_ttother, // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psScaleWeight,      // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        scale,              // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower
        scale_ttbb,         // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (ttbb process)
        scale_ttb,          // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (ttb process)
        scale_tt2b,         // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (tt2b process)
        scale_ttcc,         // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (ttcc process)
        scale_ttother,      // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (tt+light jets process)
        bFrag,              // b quark fragmentation function (up and down)
        bFrag_central,      // b quark fragmentation function (so called central: it is not nominal)
        bFrag_Peterson,     // b quark fragmentation function (alternative Peterson fragmentation function)
        bSemilep,           // b hadron semileptonic decay branching ratios
        ueTune,             // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune
        ueTune_ttbb,        // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_tt2b,        // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_ttb,         // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_ttcc,        // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_ttother,     // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        powhegv2,           // POWHEGV2 event generator matched to PYTHIA8 shower
        powheg,             // POWHEG event generator matched to PYTHIA shower
        powhegv2Herwig,     // POWHEGV2 event generator matched to HERWIG++ shower
        powhegHerwig,       // POWHEG event generator matched to HERWIG shower
        powhegHelac,        // POWHEG ME generator with Helac as one loop provider (usually interfaced with Pythia8)
        powhegOpenloops,     // POWHEG ME generator with Openloops as one loop provider
        amcatnlofxfx,       // aMC@NLO event generator matched via FxFx to PYTHIA8 shower
        mcatnlo,            // MC@NLO event generator
        madgraphmlm,        // Madgraph event generator matched via MLM to PYTHIA8 shower
        cp5,                // cp5 tune for 2016 94X samples, used for Btag reweighting
        perugia11,          // Perugia11 parton shower tune
        perugia11NoCR,      // Perugia11 parton shower tune, no colour-reconnection
        alphasPdf,          // Variation of strong coupling in nominal PDF
        l1prefiring,          // Variation of L1 prefiring correction for 2016/2017
        normPdfGg,          // Normalization due to pdf in gg production
        normPdfGq,          // Normalization due to pdf in gqbar production
        normPdfQq,          // Normalization due to pdf in qqbar production
        normPdfTth,         // Normalization due to pdf gg prodution applied to ttH only (decoupled sys. from tt+XX)
        pdf_pca_1,          // PDF variation from PCA decomposition
        pdf_pca_2,          // PDF variation from PCA decomposition
        pdf,                // PDF variations
        xsec_vv,            // cross-section uncertainty of diboson process
        xsec_v,             // cross-section uncertainty of single boson process (if w+jets and z+jets grouped)
        xsec_w,             // cross-section uncertainty of w+jets process
        xsec_z,             // cross-section uncertainty of z+jets process
        unclustered,        // unclustered energy variation in the MET
        uncorrelatedType,   // Variations that aren't correlated between the FullRunII datasets
        closure,            // Closure test
        allAvailable,       // All systematics which are available
        all,                // All allowed systematics
        undefinedType       // No systematic defined (also not nominal)
    };



    /// Convert a type from string to enum
    Type convertType(const TString& type);

    /// Convert a type from enum to string
    TString convertType(const Type& type);

    /// Convert a vector of types from string to enum
    std::vector<Type> convertType(const std::vector<TString>& types);

    /// Convert a vector of types from string to enum
    std::vector<Type> convertType(const std::vector<std::string>& types);

    /// Convert a vector of types from string to enum
    std::vector<TString> convertType(const std::vector<Type>& types);

    /// All variations as needed in any part of the framework
    enum Variation{up, down, central, undefinedVariation};

    /// Convert a variation from string to enum
    Variation convertVariation(const TString& variation);

    /// Convert a variation from enum to string
    TString convertVariation(const Variation& variation);

    /// Convert a vector of variations from string to enum
    std::vector<Variation> convertVariation(const std::vector<TString>& variations);

    /// Convert a vector of variations from string to enum
    std::vector<Variation> convertVariation(const std::vector<std::string>& variations);

    /// Convert a vector of variations from string to enum
    std::vector<TString> convertVariation(const std::vector<Variation>& variations);







    /// Define for which systematics up/down variations are allowed
    const std::vector<Type> upDownTypes{
        lept, trig, trigEta, pu,
        ele, eleID, eleIDSyst, eleIDStat, eleReco, eleRecoSyst, eleRecoStat,
        muon, muonID, muonIDSyst, muonIDStat, muonIso, muonIsoSyst, muonIsoStat,
        eleScaleSyst, eleScaleGain, eleScaleStat, eleScaleEt, eleSmearingPhi, eleSmearingRho,
        muonScaleEwk, muonScaleStat, muonScaleZpt, muonScaleDeltaM, muonScaleEwk2,
        eleScaleSmearing, muonScale,
        dy, bg, kin,
        dynorm,
        jetPileupID,
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        btagBeff, btagCeff, btagLeff,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
        jer, jerEta0, jerEta1, jerEta2Pt0, jerEta2Pt1, jerEta3Pt0, jerEta3Pt1,
        jes, jesAbsoluteStat, jesAbsoluteScale, jesAbsoluteFlavMap, jesAbsoluteMPFBias, jesFragmentation, jesSinglePionECAL,
        jesSinglePionHCAL, jesFlavorQCD, jesTimePtEta, jesRelativeJEREC1, jesRelativeJEREC2, jesRelativeJERHF, jesRelativePtBB, jesRelativePtEC1,
        jesRelativePtEC2, jesRelativePtHF, jesRelativeBal, jesRelativeSample, jesRelativeFSR, jesRelativeStatFSR, jesRelativeStatEC, jesRelativeStatHF, jesPileUpDataMC,
        jesPileUpPtRef, jesPileUpPtEC1, jesPileUpPtEC2, jesPileUpPtHF, jesPileUpPtBB, jesPileUpMuZero, jesPileUpEnvelope, jesSubTotalPileUp,
        jesSubTotalRelative, jesSubTotalPt, jesSubTotalScale, jesSubTotalMC, jesSubTotalAbsolute, jesTotalNoFlavor,
        //"TotalNoTime"       // ignoring for the moment these subtotals
        //"TotalNoFlavorNoTime"
        //"Time A-D
        jesFlavorZJet, jesFlavorPhotonJet, jesFlavorPureGluon, jesFlavorPureQuark, jesFlavorPureCharm, jesFlavorPureBottom, jesFlavorRealistic,
        jesCorrelationGroupMPFInSitu, jesCorrelationGroupIntercalibration, jesCorrelationGroupbJES, jesCorrelationGroupFlavor, jesCorrelationGroupUncorrelated,
        jesHttEta0to5, jesHttEta0to3, jesHttEta3to5, jesHttEC2,
        jesAbsolute, jesAbsoluteYear, jesBBEC1, jesBBEC1Year, jesEC2, jesEC2Year, jesHF, jesHFYear,
        jesUserDefinedHEM1516,
        frac_tthf, frac_ttother,
        lumi,
        xsec_ttbb, xsec_ttb, xsec_tt2b, xsec_ttcc, xsec_ttother,
        xsec_ttH, xsec_ttZ, xsec_ttW, xsec_ttG, xsec_ttV, xsec_tt, xsec_t,
        xsec_v, xsec_vv,
        xsec_w, xsec_z,
        topPtTheory, topPtFitparam0, topPtFitparam1, topPt,
        mass,
        match,
        match_ttbb, match_ttb, match_tt2b, match_ttcc, match_ttother,
        meScale, meScale_ttbb, meScale_ttb, meScale_tt2b, meScale_ttcc, meScale_ttother,
        meFacScale, meFacScale_ttbb, meFacScale_ttb, meFacScale_tt2b, meFacScale_ttcc, meFacScale_ttother,
        meRenScale, meRenScale_ttbb, meRenScale_ttb, meRenScale_tt2b, meRenScale_ttcc, meRenScale_ttother,
        meFacScale_tt, meFacScale_z, meFacScale_w, meFacScale_st, meFacScale_vv, meFacScale_ttv,
        meFacScale_ttz, meFacScale_ttw, meFacScale_ttg, meFacScale_ttdm, meFacScale_htott_res, meFacScale_htott_int,
        meRenScale_tt, meRenScale_z, meRenScale_w, meRenScale_st, meRenScale_vv, meRenScale_ttv,
        meRenScale_ttz, meRenScale_ttw, meRenScale_ttg, meRenScale_ttdm, meRenScale_htott_res, meRenScale_htott_int,
        psScale,
        psScale_ttbb, psScale_ttb, psScale_tt2b, psScale_ttcc, psScale_ttother,
        psISRScale,
        psISRScale_ttbb, psISRScale_ttb, psISRScale_tt2b, psISRScale_ttcc, psISRScale_ttother,
        psFSRScale,
        psFSRScale_ttbb, psFSRScale_ttb, psFSRScale_tt2b, psFSRScale_ttcc, psFSRScale_ttother,
        scale,
        //scale_ttbb, scale_ttb, scale_tt2b, scale_ttcc, scale_ttother,
        bFrag, bSemilep,
        ueTune,
        ueTune_ttbb, ueTune_ttb, ueTune_tt2b, ueTune_ttbb, ueTune_ttcc, ueTune_ttother,
        alphasPdf,l1prefiring,
        normPdfGg, normPdfGq, normPdfQq, normPdfTth,
        pdf_pca_1, pdf_pca_2, pdf, psScaleWeight,
        unclustered
    };

    /// Define for which systematics central variations are allowed
    /// This is also used to identify for which systematics variation numbers can be assigned
    const std::vector<Type> centralTypes{
        pdf, psScaleWeight
    };



    /// Check the validity of a variation for a given type
    void isValid(const Type& type, const Variation& variation, const int variationNumber =-1);


    const std::vector<Type> jesTypes{
        jes, jesAbsoluteStat, jesAbsoluteScale, jesAbsoluteFlavMap, jesAbsoluteMPFBias, jesFragmentation, jesSinglePionECAL,
        jesSinglePionHCAL, jesFlavorQCD, jesTimePtEta, jesRelativeJEREC1, jesRelativeJEREC2, jesRelativeJERHF, jesRelativePtBB, jesRelativePtEC1,
        jesRelativePtEC2, jesRelativePtHF, jesRelativeBal, jesRelativeSample, jesRelativeFSR, jesRelativeStatFSR, jesRelativeStatEC, jesRelativeStatHF, jesPileUpDataMC,
        jesPileUpPtRef, jesPileUpPtEC1, jesPileUpPtEC2, jesPileUpPtHF, jesPileUpPtBB, jesPileUpMuZero, jesPileUpEnvelope, jesSubTotalPileUp,
        jesSubTotalRelative, jesSubTotalPt, jesSubTotalScale, jesSubTotalMC, jesSubTotalAbsolute, jesTotalNoFlavor,
        //"TotalNoTime"       // ignoring for the moment these subtotals
        //"TotalNoFlavorNoTime"
        //"Time A-D
        jesFlavorZJet, jesFlavorPhotonJet, jesFlavorPureGluon, jesFlavorPureQuark, jesFlavorPureCharm, jesFlavorPureBottom, jesFlavorRealistic,
        jesCorrelationGroupMPFInSitu, jesCorrelationGroupIntercalibration, jesCorrelationGroupbJES, jesCorrelationGroupFlavor, jesCorrelationGroupUncorrelated,
        jesHttEta0to5, jesHttEta0to3, jesHttEta3to5, jesHttEC2,
        jesAbsolute, jesAbsoluteYear, jesBBEC1, jesBBEC1Year, jesEC2, jesEC2Year, jesHF, jesHFYear,
        jesUserDefinedHEM1516,
    };


    /// Define b-tag systematics, valid for all b-tag corrections
    const std::vector<Type> btagTypes{
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        btagBeff, btagCeff, btagLeff,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
    };

    /// Define b-tag systematics, valid for b-tag corrections concerning discriminator reweighting
    const std::vector<Type> btagDiscriminatorReweightTypes{
        btag, btagLjet,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
    };

    /// Define b-tag systematics, valid for b-tag corrections concerning efficiency
    const std::vector<Type> btagEfficiencyCorrectionTypes{
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        btagBeff, btagCeff, btagLeff,
    };

    /// Define ttbar systematics, i.e. variations of the ttbar sample (e.g. mass or scale variations)
    const std::vector<Type> ttbarTypes{
        topPtTheory, topPtFitparam0, topPtFitparam1, topPt,
        mass,
        match,
        //match_ttbb, match_ttb, match_tt2b, match_ttcc, match_ttother,
        erdon, erdonretune, gluonmovetune,
        meScale, meScale_ttbb, meScale_ttb, meScale_tt2b, meScale_ttcc, meScale_ttother,
        meFacScale, meFacScale_ttbb, meFacScale_ttb, meFacScale_tt2b, meFacScale_ttcc, meFacScale_ttother,
        meRenScale, meRenScale_ttbb, meRenScale_ttb, meRenScale_tt2b, meRenScale_ttcc, meRenScale_ttother,
        psScale, psScale_ttbb, psScale_ttb, psScale_tt2b, psScale_ttcc, psScale_ttother,
        psISRScale,
        //psISRScale_ttbb, psISRScale_ttb, psISRScale_tt2b, psISRScale_ttcc, psISRScale_ttother,
        psFSRScale,
        //psFSRScale_ttbb, psFSRScale_ttb, psFSRScale_tt2b, psFSRScale_ttcc, psFSRScale_ttother,
        scale, scale_ttbb, scale_ttb, scale_tt2b, scale_ttcc, scale_ttother,
        bFrag, bFrag_central, bFrag_Peterson, bSemilep,
        ueTune,
        ueTune_ttbb, ueTune_ttb, ueTune_tt2b, ueTune_ttcc, ueTune_ttother,
        powhegv2, powheg, powhegv2Herwig, powhegHerwig, powhegHelac, powhegOpenloops, amcatnlofxfx, mcatnlo, madgraphmlm, cp5, perugia11, perugia11NoCR,
        alphasPdf, pdf, psScaleWeight,
        closure,
    };

    /// Define cross-section uncertainty systematics, which use nominal samples, and change only the scaling
    const std::vector<Type> crossSectionTypes{
        xsec_ttbb, xsec_ttb, xsec_tt2b, xsec_ttcc, xsec_ttother,
        xsec_ttH, xsec_ttZ, xsec_ttW, xsec_ttG, xsec_ttV, xsec_tt, xsec_t,
        xsec_v, xsec_vv, xsec_w, xsec_z
    };

    /// Define uncertainties due to tt+HF fraction scale factor from the fit, which use nominal samples, and change only the scaling
    const std::vector<Type> tthfFractionTypes{
        frac_tthf, frac_ttother
    };

    /// Define systematics that do not require dedicated root files
    const std::vector<Type> fileIndependentTypes{
        xsec_ttbb, xsec_ttb, xsec_tt2b, xsec_ttcc, xsec_ttother,
        xsec_ttH, xsec_ttZ, xsec_ttW, xsec_ttV, xsec_tt, xsec_t,
        xsec_v, xsec_vv,
        xsec_w, xsec_z,
        dynorm,
        frac_tthf, frac_ttother,
        lumi,
        normPdfGg, normPdfGq, normPdfQq, normPdfTth
    };

    const std::vector<Type> uncorrelatedTypes{
        jer, jerEta0, jerEta1, jerEta2Pt0, jerEta2Pt1, jerEta3Pt0, jerEta3Pt1,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
        jesAbsoluteStat, jesRelativeStatEC, jesRelativeStatFSR,
        jesRelativeJEREC1, jesRelativeJEREC2,
        jesRelativePtEC1, jesRelativePtEC2, jesTimePtEta,
        jesAbsoluteYear, jesBBEC1Year, jesEC2Year, jesHFYear, jesRelativeSample

    };

    /// Class for proper handling of systematic
    class Systematic{

     public:

        Systematic();

        Systematic(const Type& type, const Variation& variation, const int variationNumber =-1);

        Systematic(const TString& systematicName);

        ~Systematic(){}

        bool operator<(const Systematic& rhs)const{return this->name() < rhs.name();}

        TString name()const;

        Type type()const{return type_;}

        Variation variation()const{return variation_;}

        int variationNumber()const{return variationNumber_;}

        std::string type_str() const { return std::string(convertType(type_).Data()); }

    private:

        Type type_;

        Variation variation_;

        int variationNumber_;
    };

    /// Set all systematics from a list of allowed types, using the defined variations
    std::vector<Systematic> allowedSystematicsAnalysis(const std::vector<Type>& allowedTypes);

    /// Set up systematics from vector of systematicNames
    std::vector<Systematic> setSystematics(const std::vector<std::string>& systematicNames);

    /// Set up systematic for nominal (i.e. no systematic variation)
    Systematic nominalSystematic();

    /// Set up undefined systematic
    Systematic undefinedSystematic();
}









/// Namespace to treat decay channels as enum types
namespace Channel {

    /// All dileptonic decay channels as needed in any part of the framework
    enum Channel{ee, emu, mumu, sameflavor, combined, tautau, se, smu, met, undefined};


    /// All dileptonic decay channels allowed for analysis step
    /// (allow undefined to select all channels if no option is set, i.e. option is empty)
    const std::vector<Channel> allowedChannelsAnalysis
        {ee, emu, mumu, undefined};

    /// All dileptonic decay channels allowed for plotting step
    const std::vector<Channel> allowedChannelsPlotting
        {ee, emu, mumu, sameflavor, combined};

    /// Real analysis channels, i.e. all channels which describe a real final state
    const std::vector<Channel> realChannels
        {ee, emu, mumu};

    /// Possible Drell-Yan decay channels
    const std::vector<Channel> dyDecayChannels
        {ee, mumu, tautau};


    /// Convert a channel from string to enum
    Channel convert(const TString& channel);

    /// Convert a channel from enum to string
    TString convert(const Channel& channel);

    /// Return the label of a channel as used for drawing
    TString label(const Channel& channel);

    /// Convert a vector of channels from string to enum
    std::vector<Channel> convert(const std::vector<TString>& channels);

    /// Convert a vector of channels from string to enum
    std::vector<Channel> convert(const std::vector<std::string>& channels);

    /// Convert a vector of channels from string to enum
    std::vector<TString> convert(const std::vector<Channel>& channels);
}







namespace common{

    /// Create and assign an output folder depending on the channel and systematic
    TString assignFolder(const char* baseDir, const Channel::Channel& channel, const Systematic::Systematic& systematic, const char* subDir ="");

    /// Access an already existing input folder
    TString accessFolder(const char* baseDir, const Channel::Channel& channel,
                         const Systematic::Systematic& systematic, const bool allowNonexisting =false);

    /// Access the real final state from a filename, ie. only "ee", "emu", "mumu", but not "combined"
    Channel::Channel finalState(const TString& filename);

    /// Find file list for a given channel and systematic, and return its name
    /// In case it does not exist, return empty string
    TString findFilelist(const TString& filelistDirectory,
                         const Channel::Channel& channel,
                         const Systematic::Systematic& systematic);

    /// Find from vector of given systematics those for which a file list exists for all given channels
    std::vector<Systematic::Systematic> findSystematicsFromFilelists(const TString& filelistDirectory,
                                                                     const std::vector<Channel::Channel>& v_channel,
                                                                     const std::vector<Systematic::Systematic>& v_systematic);

    /// Read the file list for given channel and systematic, and return the input file names
    /// In case a vector of patterns is specified, only files containing this pattern in the full path name will be read
    std::vector<TString> readFilelist(const TString& filelistDirectory,
                                      const Channel::Channel& channel,
                                      const Systematic::Systematic& systematic,
                                      const std::vector<TString>& v_pattern =std::vector<TString>());

    /// Read a file for given file name, and return the lines each as element in vector
    /// In case a vector of patterns is specified, only lines containing this pattern will be read
    std::vector<TString> readFile(const TString& filename,
                                  const std::vector<TString>& v_pattern =std::vector<TString>());
}








#endif
