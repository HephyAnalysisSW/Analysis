#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>

#include <TLorentzVector.h>
#include <TMath.h>

#include "Analysis/TopReco/interface/analysisUtils.h"
#include "Analysis/TopReco/src/classes.h"





// --- Several conversion functions -------------------------------------------------------------------------------------

void common::LVtof4(const LV& lv, float* f)
{
    f[0] = lv.E();
    f[1] = lv.Px();
    f[2] = lv.Py();
    f[3] = lv.Pz();
}

std::string common::f2s(const float& f)
{
    char result[100];
    if(std::abs(f) < 5) {
        std::sprintf(result, "%.3f", f);
        std::string s = std::string(result);
        while (s.length() > 0 && s[s.length()-1] == '0') s.erase(s.end()-1);
        if (s.length() > 0 && s[s.length()-1] == '.') s.erase(s.end()-1);
        return s;
    }
    else {
        std::sprintf(result, "%.0f", f);
        return std::string(result);
    }
}

void common::LVtod4(const LV& lv, double* d)
{
    d[0] = lv.E();
    d[1] = lv.Px();
    d[2] = lv.Py();
    d[3] = lv.Pz();
}

std::string common::d2s(const double& d)
{
    char result[100];
    if(std::abs(d) < 5) {
        std::sprintf(result, "%.3f", d);
        std::string s = std::string(result);
        while (s.length() > 0 && s[s.length()-1] == '0') s.erase(s.end()-1);
        if (s.length() > 0 && s[s.length()-1] == '.') s.erase(s.end()-1);
        return s;
    }
    else {
        std::sprintf(result, "%.0f", d);
        return std::string(result);
    }
}



const TLorentzVector common::LVtoTLV(const LV& lv)
{
    return TLorentzVector(lv.X(), lv.Y(), lv.Z(), lv.T());
}



const LV common::TLVtoLV(const TLorentzVector& lv)
{
    LV result;
    result.SetXYZT(lv.X(), lv.Y(), lv.Z(), lv.T());
    return result;
}



const std::vector<TLorentzVector> common::VLVtoVTLV(const VLV& vlv)
{
  std::vector<TLorentzVector> result;

  for(unsigned int i = 0; i < vlv.size(); i++){

    TLorentzVector tlv;
    tlv = TLorentzVector(vlv[i].Px(), vlv[i].Py(), vlv[i].Pz(), vlv[i].E());

    result.push_back(tlv);
  }

  return result;
}




// --- Functions concerning the treatment of indices of vectors (for working with data stored in nTuple branches) -------------

std::vector<int> common::mergeIndices(const std::vector<int>& v_index1, const std::vector<int>& v_index2, const bool allowOverlap)
{
    std::vector<int> result(v_index1);

    for(const int index : v_index2){
        // Check if index is already contained in first vector
        if(std::find(v_index1.begin(), v_index1.end(), index) != v_index1.end()){
            if(allowOverlap) continue;
            else{
                std::cerr<<"ERROR in common::mergeIndices()! Two collections should be merged, but they overlap (not allowed)\n...break\n"<<std::endl;
                exit(92);
            }
        }
        else result.push_back(index);
    }

    return result;
}



std::vector<float> common::parametersLV(const VLV& v_lv, const common::LVParameter& parameter)
{
    std::vector<float> v_variable;
    for(const LV& lv : v_lv){
        if(parameter == LVpt) v_variable.push_back(lv.pt());
        else if (parameter == LVeta) v_variable.push_back(lv.eta());
        else if (parameter == LVet) v_variable.push_back(lv.Et());
        else{
            std::cerr<<"Error in common::parametersLV()! Lorentz vector parameter is not implemented\n...break\n";
            exit(638);
        }
    }
    return v_variable;
}

void common::orderIndices(int& index1, int& index2, const VLV& v_lv, const common::LVParameter& parameter, const bool absoluteValue)
{
    const std::vector<float> v_variable = parametersLV(v_lv, parameter);
    orderIndices(index1, index2, v_variable, absoluteValue);
}

void common::orderIndices(std::vector<int>& v_index, const VLV& v_lv, const common::LVParameter& parameter, const bool absoluteValue)
{
    const std::vector<float> v_variable = parametersLV(v_lv, parameter);
    orderIndices(v_index, v_variable, absoluteValue);
}

void common::selectIndices(std::vector<int>& v_index, const VLV& v_lv, const common::LVParameter& parameter, const float cutValue, const bool lowerThreshold)
{
    const std::vector<float> v_variable = parametersLV(v_lv, parameter);
    selectIndices(v_index, v_variable, cutValue, lowerThreshold);
}

int common::extremumIndex(const VLV& v_lv, const common::LVParameter& parameter, const bool maximumValue)
{
    const std::vector<float> v_variable = parametersLV(v_lv, parameter);
    return extremumIndex(v_variable, maximumValue);
}

void common::orderLV(LV& lv1, LV& lv2, const LV& inputLv1, const LV& inputLv2, const common::LVParameter& parameter, const bool absoluteValue)
{
    float variable1;
    float variable2;
    if(parameter == LVpt){
        variable1 = inputLv1.pt();
        variable2 = inputLv2.pt();
    }
    else if(parameter == LVeta){
        variable1 = inputLv1.eta();
        variable2 = inputLv2.eta();
    }
    else if(parameter == LVet){
        variable1 = inputLv1.Et();
        variable2 = inputLv2.Et();
    }
    else{
        std::cerr<<"Error in common::orderLV()! Lorentz vector parameter is not implemented\n...break\n";
        exit(639);
    }

    if(absoluteValue){
        variable1 = std::abs(variable1);
        variable2 = std::abs(variable2);
    }

    if (variable1 > variable2) {
        lv1 = inputLv1;
        lv2 = inputLv2;
    } else {
        lv1 = inputLv2;
        lv2 = inputLv1;
    }
}

void common::selectLeptonIndicesID(std::vector<int>& v_index, const std::vector<int>* const v_lepID, const int min_lepID, const int only_absPdgID, const std::vector<int>* const v_pdgID)
{
  std::vector<int> result;

  if(v_lepID == nullptr)
  {
    const std::string log("null pointer to vector of \"lepton-ID\" values");
    throw std::runtime_error("common::selectLeptonIndicesID -- "+log);
  }

  if(only_absPdgID > 0)
  {
    if(v_pdgID == nullptr)
    {
      const std::string log("null pointer to vector of \"lepPdgIds\" values");
      throw std::runtime_error("common::selectLeptonIndicesID -- "+log);
    }
  }

  for(const int index : v_index)
  {
    if(index >= int(v_lepID->size()))
    {
      const std::string log("lepton index ("+std::to_string(index)+") out-of-range in vector of \"lepton-ID\" values");
      throw std::runtime_error("common::selectLeptonIndicesID -- "+log);
    }

    // if only applying selection on one lepton flavor, all other leptons are retained (pass=true)
    if(only_absPdgID > 0)
    {
      if(index >= int(v_pdgID->size()))
      {
        const std::string log("lepton index ("+std::to_string(index)+") out-of-range in vector of \"lepPdgIds\" values");
        throw std::runtime_error("common::selectLeptonIndicesID -- "+log);
      }

      if(std::abs(v_pdgID->at(index)) != only_absPdgID)
      {
        result.emplace_back(index);
        continue;
      }
    }

    if(v_lepID->at(index) >= min_lepID)
    {
      result.emplace_back(index);
    }
  }

  v_index.clear();
  v_index = result;
}

void common::selectMuonIndicesID(std::vector<int>& v_index, const std::vector<int>* const v_lepID, const int min_lepID, const std::vector<int>* const v_pdgID)
{
  return common::selectLeptonIndicesID(v_index, v_lepID, min_lepID, 13, v_pdgID);
}

void common::selectElectronIndicesID(std::vector<int>& v_index, const std::vector<int>* const v_lepID, const int min_lepID, const std::vector<int>* const v_pdgID)
{
  return common::selectLeptonIndicesID(v_index, v_lepID, min_lepID, 11, v_pdgID);
}

// Function to correct PF-MET phi
LV common::metWithXYCorrection(const LV& uncormet, const int runnb, const Era::Era& era, const bool isMC, const int nPV)
{
  enum TheRunEra{y2016B,y2016C,y2016D,y2016E,y2016F,y2016G,y2016H,y2017B,y2017C,y2017D,y2017E,y2017F,y2018A,y2018B,y2018C,y2018D,y2016MC,y2017MC,y2018MC};

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#xy_Shift_Correction_MET_phi_modu
  // https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection.h

  int npv = nPV;
  if(npv>100) npv = 100;
  int runera = -1;
  bool usemetv2 = false;

  int year = -1;

  if(era == Era::run2_13tev_25ns || era == Era::run2_13tev_2015_25ns || era == Era::run2_13tev_2016_25ns ||era== Era::run2_13tev_2016){
      year = 2016;
  }
  else if(era == Era::run2_13tev_2017)
  {
      year = 2017;
  }
  else if(era == Era::run2_13tev_2018)
  {
    year = 2018;
  }else{
      const std::string log("could not convert era to year");
      throw std::runtime_error("common::metXYCorrection -- "+log);
  }

  if(isMC && year == 2016) runera = y2016MC;
  else if(isMC && year == 2017) {runera = y2017MC; usemetv2 =true;}
  else if(isMC && year == 2018) runera = y2018MC;

  else if(!isMC && runnb >=272007 &&runnb<=275376  ) runera = y2016B;
  else if(!isMC && runnb >=275657 &&runnb<=276283  ) runera = y2016C;
  else if(!isMC && runnb >=276315 &&runnb<=276811  ) runera = y2016D;
  else if(!isMC && runnb >=276831 &&runnb<=277420  ) runera = y2016E;
  else if(!isMC && runnb >=277772 &&runnb<=278808  ) runera = y2016F;
  else if(!isMC && runnb >=278820 &&runnb<=280385  ) runera = y2016G;
  else if(!isMC && runnb >=280919 &&runnb<=284044  ) runera = y2016H;

  else if(!isMC && runnb >=297020 &&runnb<=299329 ){ runera = y2017B; usemetv2 =true;}
  else if(!isMC && runnb >=299337 &&runnb<=302029 ){ runera = y2017C; usemetv2 =true;}
  else if(!isMC && runnb >=302030 &&runnb<=303434 ){ runera = y2017D; usemetv2 =true;}
  else if(!isMC && runnb >=303435 &&runnb<=304826 ){ runera = y2017E; usemetv2 =true;}
  else if(!isMC && runnb >=304911 &&runnb<=306462 ){ runera = y2017F; usemetv2 =true;}

  else if(!isMC && runnb >=315252 &&runnb<=316995 ) runera = y2018A;
  else if(!isMC && runnb >=316998 &&runnb<=319312 ) runera = y2018B;
  else if(!isMC && runnb >=319313 &&runnb<=320393 ) runera = y2018C;
  else if(!isMC && runnb >=320394 &&runnb<=325273 ) runera = y2018D;

  else {
    //Couldn't find data/MC era => no correction applied
    const std::string log("could not determine data/MC era");
    throw std::runtime_error("common::metXYCorrection -- "+log);
  }

  float METxcorr(0.),METycorr(0.);

  if(!usemetv2){//Current recommendation for 2016 and 2018
    if(runera==y2016B) METxcorr = -(-0.0478335*npv -0.108032);
    if(runera==y2016B) METycorr = -(0.125148*npv +0.355672);
    if(runera==y2016C) METxcorr = -(-0.0916985*npv +0.393247);
    if(runera==y2016C) METycorr = -(0.151445*npv +0.114491);
    if(runera==y2016D) METxcorr = -(-0.0581169*npv +0.567316);
    if(runera==y2016D) METycorr = -(0.147549*npv +0.403088);
    if(runera==y2016E) METxcorr = -(-0.065622*npv +0.536856);
    if(runera==y2016E) METycorr = -(0.188532*npv +0.495346);
    if(runera==y2016F) METxcorr = -(-0.0313322*npv +0.39866);
    if(runera==y2016F) METycorr = -(0.16081*npv +0.960177);
    if(runera==y2016G) METxcorr = -(0.040803*npv -0.290384);
    if(runera==y2016G) METycorr = -(0.0961935*npv +0.666096);
    if(runera==y2016H) METxcorr = -(0.0330868*npv -0.209534);
    if(runera==y2016H) METycorr = -(0.141513*npv +0.816732);
    if(runera==y2017B) METxcorr = -(-0.259456*npv +1.95372);
    if(runera==y2017B) METycorr = -(0.353928*npv -2.46685);
    if(runera==y2017C) METxcorr = -(-0.232763*npv +1.08318);
    if(runera==y2017C) METycorr = -(0.257719*npv -1.1745);
    if(runera==y2017D) METxcorr = -(-0.238067*npv +1.80541);
    if(runera==y2017D) METycorr = -(0.235989*npv -1.44354);
    if(runera==y2017E) METxcorr = -(-0.212352*npv +1.851);
    if(runera==y2017E) METycorr = -(0.157759*npv -0.478139);
    if(runera==y2017F) METxcorr = -(-0.232733*npv +2.24134);
    if(runera==y2017F) METycorr = -(0.213341*npv +0.684588);
    if(runera==y2018A) METxcorr = -(0.362865*npv -1.94505);
    if(runera==y2018A) METycorr = -(0.0709085*npv -0.307365);
    if(runera==y2018B) METxcorr = -(0.492083*npv -2.93552);
    if(runera==y2018B) METycorr = -(0.17874*npv -0.786844);
    if(runera==y2018C) METxcorr = -(0.521349*npv -1.44544);
    if(runera==y2018C) METycorr = -(0.118956*npv -1.96434);
    if(runera==y2018D) METxcorr = -(0.531151*npv -1.37568);
    if(runera==y2018D) METycorr = -(0.0884639*npv -1.57089);
    if(runera==y2016MC) METxcorr = -(-0.195191*npv -0.170948);
    if(runera==y2016MC) METycorr = -(-0.0311891*npv +0.787627);
    if(runera==y2017MC) METxcorr = -(-0.217714*npv +0.493361);
    if(runera==y2017MC) METycorr = -(0.177058*npv -0.336648);
    if(runera==y2018MC) METxcorr = -(0.296713*npv -0.141506);
    if(runera==y2018MC) METycorr = -(0.115685*npv +0.0128193);
  }
  else {//these are the corrections for v2 MET recipe (currently recommended for 2017)
    if(runera==y2016B) METxcorr = -(-0.0374977*npv +0.00488262);
    if(runera==y2016B) METycorr = -(0.107373*npv +-0.00732239);
    if(runera==y2016C) METxcorr = -(-0.0832562*npv +0.550742);
    if(runera==y2016C) METycorr = -(0.142469*npv +-0.153718);
    if(runera==y2016D) METxcorr = -(-0.0400931*npv +0.753734);
    if(runera==y2016D) METycorr = -(0.127154*npv +0.0175228);
    if(runera==y2016E) METxcorr = -(-0.0409231*npv +0.755128);
    if(runera==y2016E) METycorr = -(0.168407*npv +0.126755);
    if(runera==y2016F) METxcorr = -(-0.0161259*npv +0.516919);
    if(runera==y2016F) METycorr = -(0.141176*npv +0.544062);
    if(runera==y2016G) METxcorr = -(0.0583851*npv +-0.0987447);
    if(runera==y2016G) METycorr = -(0.0641427*npv +0.319112);
    if(runera==y2016H) METxcorr = -(0.0706267*npv +-0.13118);
    if(runera==y2016H) METycorr = -(0.127481*npv +0.370786);
    if(runera==y2017B) METxcorr = -(-0.19563*npv +1.51859);
    if(runera==y2017B) METycorr = -(0.306987*npv +-1.84713);
    if(runera==y2017C) METxcorr = -(-0.161661*npv +0.589933);
    if(runera==y2017C) METycorr = -(0.233569*npv +-0.995546);
    if(runera==y2017D) METxcorr = -(-0.180911*npv +1.23553);
    if(runera==y2017D) METycorr = -(0.240155*npv +-1.27449);
    if(runera==y2017E) METxcorr = -(-0.149494*npv +0.901305);
    if(runera==y2017E) METycorr = -(0.178212*npv +-0.535537);
    if(runera==y2017F) METxcorr = -(-0.165154*npv +1.02018);
    if(runera==y2017F) METycorr = -(0.253794*npv +0.75776);
    if(runera==y2018A) METxcorr = -(0.362642*npv +-1.55094);
    if(runera==y2018A) METycorr = -(0.0737842*npv +-0.677209);
    if(runera==y2018B) METxcorr = -(0.485614*npv +-2.45706);
    if(runera==y2018B) METycorr = -(0.181619*npv +-1.00636);
    if(runera==y2018C) METxcorr = -(0.503638*npv +-1.01281);
    if(runera==y2018C) METycorr = -(0.147811*npv +-1.48941);
    if(runera==y2018D) METxcorr = -(0.520265*npv +-1.20322);
    if(runera==y2018D) METycorr = -(0.143919*npv +-0.979328);
    if(runera==y2016MC) METxcorr = -(-0.159469*npv +-0.407022);
    if(runera==y2016MC) METycorr = -(-0.0405812*npv +0.570415);
    if(runera==y2017MC) METxcorr = -(-0.182569*npv +0.276542);
    if(runera==y2017MC) METycorr = -(0.155652*npv +-0.417633);
    if(runera==y2018MC) METxcorr = -(0.299448*npv +-0.13866);
    if(runera==y2018MC) METycorr = -(0.118785*npv +0.0889588);
  }

  float CorrectedMET_x = uncormet.E() *cos( uncormet.Phi())+METxcorr;
  float CorrectedMET_y = uncormet.E() *sin( uncormet.Phi())+METycorr;

  float CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  float CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;

  LV cormet;
  cormet.SetCoordinates(CorrectedMET, uncormet.Eta() ,CorrectedMETPhi, uncormet.M());
  return cormet;
}
