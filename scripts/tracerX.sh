#! /bin/bash
FOLDER="$2"
SIMCHA="$1"
GENELIST="$3"
SEGMENTEDFOLDER="$2/segmented"
MEDICCFOLDER="$2/medicc2"
FITNESSFOLDER="$2/fitness"

eval "$(conda shell.bash hook)"

#Run segmentation
FILESEG=$FOLDER/simcha_ready_for_segmentation.tsv
if [ -f "$FILESEG" ] && [ "$4" == "True" ]; 
then
conda activate simcha
cd $SIMCHA
dotnet run --project SimChA -P $FOLDER/simcha_ready_for_segmentation.tsv -s -O $SEGMENTEDFOLDER
else
echo "$FILESEG does not exist"
fi

#Run medicc
FILEMED=$SEGMENTEDFOLDER/consistent_CNs.tsv
if [ -f "$FILEMED" ] && [ "$4" == "True" ]; 
then
sed -i 's/NA/0/g' "$SEGMENTEDFOLDER/consistent_CNs.tsv"
conda activate medicc_env
medicc2 "$SEGMENTEDFOLDER/consistent_CNs.tsv" $MEDICCFOLDER --plot none
else
echo "$FILEMED does not exist"
fi

#Run fitness calculation
FILEFIT=$MEDICCFOLDER/consistent_CNs_final_cn_profiles.tsv
if [ -f "$FILEFIT" ] && [ "$4" == "True" ]; 
then
conda activate simcha
dotnet run --project SimChA -P $MEDICCFOLDER/consistent_CNs_final_cn_profiles.tsv -O $FITNESSFOLDER -D $GENELIST
else
echo "$FILEFIT does not exist."
fi