### Very first implementation of Muon TnP configurtion file
variables = {
    'probe' : ["charge", "pt", "eta", "phi", "dxy", "dz", "isTracker", "isStandalone", "isGlobal", "HLTIsoMu24", "tightId", "mediumId", "looseId", "highPtId", "pfIsoId", "pfAbsIso04_all"],
    "tag"   : ["charge", "pt", "eta", "phi", "dxy", "dz", "isTracker", "isStandalone", "isGlobal", "HLTIsoMu24", "tightId", "mediumId", "looseId", "highPtId", "pfIsoId", "pfAbsIso04_all"],
    "save"  : ["pair_mass","pair_pt","pair_eta","pair_phi","eventIdx",
               "probe_charge","probe_pt","probe_eta","probe_isGlobal","probe_HLTIsoMu24","probe_dxy","probe_dz","probe_isTracker","probe_isStandalone","probe_tightId","probe_mediumId","probe_looseId","probe_highPtId","probe_pfIsoId","probe_pfAbsIso04_all","probe_phi",
               "tag_charge","tag_pt","tag_eta","tag_isGlobal","tag_HLTIsoMu24","tag_dxy","tag_dz","tag_isTracker","tag_isStandalone","tag_tightId","tag_mediumId","tag_looseId","tag_highPtId","tag_pfIsoId","tag_pfAbsIso04_all","tag_phi",
               ],
}

selection = {
    "probe" : "GeneralTrack_pt > 2 && abs(GeneralTrack_eta) < 2.4", # Muon_isStandalone
    "tag"   : "Muon_pt > 15 && abs(Muon_eta) < 2.4 &&  Muon_tightId", # Muon_isTrig
    "pair"  : "All_TPmass > 60" 
}

samples = {
    "Run2022" : {
        'input'  : "/eos/user/r/rbhattac/Muon_POG_NanoAOD_v2/Muon/NanoMuonPOGData2022F/231103_232820/0000/",
        "output" : "/eos/cms/store/group/phys_muon//sblancof/nanoAOD/tnp.parquet",
        "lumiMask": "https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json",
        "doSplit": True,
        "root_file": "tnp.root"
    }
}
