#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <string>
#include "headers.hh"

using namespace ROOT;
using namespace ROOT::VecOps;

RVecI sortedIndices(RVecF variable){
  // return sortedIndices based on variable
  return Reverse(Argsort(variable));
}

RVecI CreateTrigIndex(const RVecF Muon_eta,
		      const RVecF Muon_phi,
		      const RVecF TrigObj_eta,
		      const RVecF TrigObj_phi,
		      const float minDR)
{
  RVecI Muon_TrigIdx(Muon_eta.size(), 999);
  float dR = 999.9;
  for (int iMu=0; iMu < Muon_eta.size(); iMu++){
    float tmpDR = minDR;
    for (int iTr=0; iTr < TrigObj_eta.size(); iTr++){
      dR = DeltaR(Muon_eta[iMu], TrigObj_eta[iTr], Muon_phi[iMu], TrigObj_phi[iTr]);
      if (dR<tmpDR){
	tmpDR = dR;
	Muon_TrigIdx[iMu] = iTr;
      }
    }
  }
  return Muon_TrigIdx;
}

RVec<std::pair<int,int>> CreateTPPair(const RVec<Int_t> &Tag_muons,
                                      const RVec<Int_t> &Probe_Candidates,
                                      const int doOppositeCharge,
                                      const RVec<Int_t> &Tag_Charge,
                                      const RVec<Int_t> &Probe_charge,
                                      const int isSameCollection = true)
{
  RVec<std::pair<int,int>> TP_pairs;
  for (int iLep1 = 0; iLep1 < Tag_muons.size(); iLep1++) {
    if (!Tag_muons[iLep1]) continue;
    for(int iLep2 = 0; iLep2 < Probe_Candidates.size(); iLep2++){
      if (isSameCollection && (iLep1 == iLep2)) continue;
      if (!Probe_Candidates[iLep2]) continue;
      if (doOppositeCharge && (Tag_Charge[iLep1] == Probe_charge[iLep2])) continue;
      std::pair<int,int> TP_pair = std::make_pair(iLep1, iLep2);
      TP_pairs.push_back(TP_pair);
    }          
  }              
  return TP_pairs;
}

RVec<Float_t> getTPVariables(RVec<std::pair<int,int>> TPPairs,
                        RVec<Float_t> &Muon_pt, RVec<Float_t> &Muon_eta, RVec<Float_t> &Muon_phi,
                        RVec<Float_t> &Cand_pt, RVec<Float_t> &Cand_eta, RVec<Float_t> &Cand_phi,
                        int option = 4 /* 1 = pt, 2 = eta, 3 = phi, 4 = mass*/ )
{
  RVec<Float_t> TPVariables;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int tag_index = TPPair.first;
    int probe_index = TPPair.second;
    ROOT::Math::PtEtaPhiMVector tag(  Muon_pt[tag_index],   Muon_eta[tag_index],   Muon_phi[tag_index],   0.106); /* Muon mass is used*/
    ROOT::Math::PtEtaPhiMVector probe(Cand_pt[probe_index], Cand_eta[probe_index], Cand_phi[probe_index], 0.106); /* Muon mass is used*/
    if (option == 1) TPVariables.push_back( (tag + probe).pt() );
    else if (option == 2 ) TPVariables.push_back( (tag + probe).eta() );
    else if (option == 3) TPVariables.push_back( (tag + probe).phi() );
    else if (option == 4) TPVariables.push_back( (tag + probe).mass() );
  }
  return TPVariables;
  
}

RVecB getDuplicatedProbes(RVec<std::pair<int,int>> TPPairs, RVec<Float_t> &Cand_pt)
{
  RVecB Cand_duplicated(false, Cand_pt.size());
  RVecI Cand_index;
  for (int i=0;i<TPPairs.size();i++){
    Cand_index.push_back(TPPairs.at(i).second);
  }
  for (int i=0;i<Cand_index.size();i++){
    if (Sum(Cand_index==Cand_index[i])>1) Cand_duplicated[Cand_index[i]] = true;
  }
  return Cand_duplicated;  
}

template <typename T>
RVec<T> getVariables(RVec<std::pair<int,int>> TPPairs,
                     RVec<T>  &Cand_variable,
                     int option /*1 for tag, 2 for probe*/)
{
    RVec<T>  Variables(TPPairs.size(), 0);    
    for (int i = 0; i < TPPairs.size(); i++){
        std::pair<int, int> TPPair = TPPairs.at(i);
        T variable;
        if (option==1)      variable = Cand_variable.at(TPPair.first);
        else if (option==2) variable = Cand_variable.at(TPPair.second);
        Variables[i] = variable;
    }
    return Variables;
}


RVecI doMuonMatch(RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Track_pt, RVecF Track_eta, RVecF Track_phi, float minDR=0.1, bool doRelPt=true, float minRelPt=0.1){
  float tmp_dR = 999;
  float dR = 999;
  int candIdx = -1;
  float tmp_RelPt = 999;
  
  RVecI Track_muonIdx(Track_pt.size(), -1);
  
  if (Track_pt.size()==0 || Muon_pt.size()==0){
    return Track_muonIdx;
  }
  
  for (int i=0; i<Track_pt.size(); i++){
    tmp_dR = 999;
    tmp_RelPt = 999;
    dR = 999;
    candIdx = -1;
    
    for (int j=0; j<Muon_pt.size(); j++){
      
      tmp_dR = sqrt((Muon_eta[j]-Track_eta[i])*(Muon_eta[j]-Track_eta[i])+(Muon_phi[j]-Track_phi[i])*(Muon_phi[j]-Track_phi[i]));
      tmp_RelPt = (Muon_pt[j]-Track_pt[i])/Muon_pt[j];
      
      if (doRelPt){
	if (tmp_dR<dR && tmp_dR<minDR && tmp_RelPt<minRelPt){
	  dR = tmp_dR;
	  candIdx = j;
	}
      }else{
	if (tmp_dR<dR && tmp_dR<minDR){
	  dR = tmp_dR;
	  candIdx = j;
	}
      }
    }
    Track_muonIdx[i] = candIdx;
  } 
  return Track_muonIdx;
}

template <typename T>
RVec<T> getMuonVars(RVecI GeneralTrack_muonIdx, RVec<T> Muon_var, int isBool){
  RVec<T> GeneralTrack_var(GeneralTrack_muonIdx.size(), 0);
  for (int i=0; i<GeneralTrack_muonIdx.size(); i++){
    if (GeneralTrack_muonIdx[i]<0){
      if (!isBool) GeneralTrack_var[i] = -99;
      continue;
    }
    GeneralTrack_var[i] = Muon_var[GeneralTrack_muonIdx[i]];    
  }
  //cout << "--------------" << endl;
  //cout << Muon_var << endl;
  //cout << GeneralTrack_muonIdx << endl;
  //cout << GeneralTrack_var << endl;
  return GeneralTrack_var;
}
