#MIGEC="java -Xmx8G -jar /raid/users/pk/program_files/migec/migec-1.2.9.jar"
#$MIGEC CheckoutBatch -cute barcodes.txt checkout/
#$MIGEC Histogram checkout/ histogram/

#in=checkout/*_R*.fastq.gz
#for x in $in
#do
#$MIGEC Assemble -c --mask 0:0 $x . assembly/
#done

home=`pwd`

in=./assembly/*_R1.t5.fastq.gz

echo $in
for i in $in
do
echo $i
s=${i#./assembly/}
x=${s%_R1.t5.fastq.gz}
echo $x
#rm $x -r
mkdir $x
cd $x


#echo ""
#echo "/////---align---/////"
#/raid/users/pk/program_files/mixcr/mixcr-2.1.11/mixcr align -OvParameters.geneFeatureToAlign=VTranscript -OreadsLayout=Collinear --report alignmentReport.log "/raid/users/pk/working_dir/assembly/"$x"_R1.t5.fastq.gz" "/raid/users/pk/working_dir/assembly/"$x"_R2.t5.fastq.gz" ""$x"_output_file.vdjca"

#echo ""
#echo "/////---assemble---/////"
#/raid/users/pk/program_files/mixcr/mixcr-2.1.11/mixcr assemble -OseparateByC=true --report assembleReport.log ""$x"_output_file.vdjca" ""$x"_clones.clns"


#echo ""
#echo "/////---exportClones---/////"
#/raid/users/pk/program_files/mixcr/mixcr-2.1.11/mixcr exportClones -s -c TRB -o -t -count -fraction -nFeature CDR3 -aaFeature CDR3 -vGene -dGene -jGene -positionOf VEndTrimmed -positionOf DBeginTrimmed -positionOf DEndTrimmed -positionOf JBeginTrimmed ""$x"_clones.clns" ""$x"_clones.txt"

echo ""
echo "/////---R script run---/////"
Rscript "$home/R_ggplots_pipe.R" --no-save --no-restore --args `pwd` $x

echo "/End/"
cd ..
done



