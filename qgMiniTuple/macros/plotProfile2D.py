#!/usr/bin/env python
import os,time
for xAndyVar in ["eta_pt", 
                 "rho_pt", "rho_eta", 
                 "deltaRmin_eta", "deltaRmin_pt", 
                 "closebyJetsInCone_eta", "closebyJetsInCone_pt",
                 "closebyJetsInCone10GeV_eta", "closebyJetsInCone10GeV_pt",
                 "ratioDoubleCone_eta", "ratioDoubleCone_pt", 
                 "ptDoubleCone_eta", "ptDoubleCone_pt"]:
  os.system("./plotProfile2D " + xAndyVar + " &")
  time.sleep(15)
