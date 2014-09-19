# Script which executes crab -create and add modifies the crab scripts to use QG.db files
crab -create
DIR=$(ls -dt */ | head -1)
for f in $CMSSW_BASE/src/QGDev/qgMiniTuple/data/*.db
do
  dbfile=$(basename "$f")
  echo 'Modifying crab script to use '$dbfile
  sed -i '/cat $RUNTIME_AREA\/inputsReport.txt/a\ \necho "Copy QG .db file to working area"\ncp $SOFTWARE_DIR/src/QGDev/qgMiniTuple/data/'$dbfile' .\nls .\n' $DIR/job/CMSSW.sh
done
