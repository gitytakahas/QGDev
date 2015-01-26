#!/usr/bin/env python
import os,time
for xAndyVar in ["additionalJets_eta", "additionalJets_pt"]:
  os.system("./plotProfile2D " + xAndyVar + " &")
  time.sleep(15)
