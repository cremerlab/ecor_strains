EXPDATE="2022-12-17"
INDEX=../Ecoli_K12_MG1655_kallisto.idx;

for dir in ../../../data/transcriptomics/${EXPDATE:0}*/*/; do
    cd $dir;
    kallisto quant --plaintext -i $INDEX -o ./ ./*.gz 
    echo $dir
    cd -
done
