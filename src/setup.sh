cd /project/home/lawless/data/db/dbnsfp/


zcat  dbNSFP4.4a_variant.chr6.gz | grep TNFAIP3 > tnfaip3
lawless@sphn-momic-login:~/data/db/dbnsfp$ wc -l tnfaip3
6691 tnfaip3
zcat  dbNSFP4.4a_variant.chr6.gz | head -n1 > tnfaip3_head


zcat  dbNSFP4.4a_variant.chr4.gz | grep NFKB1 > nfkb1
wc -l nfkb1
7123 nfkb1
zcat  dbNSFP4.4a_variant.chr4.gz | head -n1 > nfkb1_head


