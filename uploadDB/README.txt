1.  Get subscribed to cms-cond-dropbox (https://e-groups.cern.ch/e-groups/)
2.  Follow the next steps on lxplus in a CMSSW environment ("cmsenv")
3.  Check payloads with "./checkPayloads.sh <db-file>"
4.  Upload to oracle://cms_orcoff_prep/CMS_COND_TEMP using "./uploadQG.py -f <db-file>"
5.  Check if succeeded using "cmscond_list_iov -c oracle://cms_orcoff_prep/CMS_COND_TEMP -t offline" ???
6.  Upload to oracle://cms_orcon_prod/CMS_COND_31X_PHYSICSTOOLS using "./uploadQG.py -f <db-file> -o"
7.  Check if succeeded using "cmscond_list_iov -c frontier://FrontierProd/CMS_COND_PAT_000 -a  | grep QGLikelihood"
