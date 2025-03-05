./getPCfromVCF.sh ../testVCFs/hlaa.vcf.gz 6:29944513-29945558 outputFromVCF
./getPCfromPlinkDirectory.sh ../testPlink/  chr6:29944513-29945558 outputFromPlink
./getSingleVariantFromPlinkDirectory.sh ../testPlink/ chr6:32529369:C:A SNVoutputFromPlink.txt

