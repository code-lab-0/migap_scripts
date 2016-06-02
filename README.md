# MiGAP scripts

## merge.pl
ORF同定、tRNA同定、16SrRNA同定、rRNA同定の結果から中間結果を生成するスクリプト。BioPerlが必要。

### Usage
```usage
perl merge.pl -f <input fasta file> -m <MetaGeneAnnotator Result> -t <tRNAscan-SE Result> -r <RNAmmer Result> -b <16S rRNA Search Result> -p <output file>
```

## fa4tr.pl
同定されたORFに対してCOG DB, RefSeq DBでBLAST検索を行った結果から、ヒットしなかったORFを抽出するスクリプト。BioPerlが必要。

### Usage
```usage
perl fa4tr.pl --fa <CDS AA Sequence> --cog <COG Search Result> --ref <RefSeq Search Result> --out <Output fasta file>
```

## merge-a.pl
中間結果とBLAST検索（COG）・BLAST検索（RefSeq）・BLAST検索（TrEMBL）の結果から最終結果を生成するスクリプト。BioPerlが必要。

### Usage
```usage
perl merge-a.pl --cog <COG Search Result> --trembl <Trembl Search Result> --ref <RefSeq Search Result> --whog <whog file in COG DB dir> --myva <myva file in COG DB dir> --gpff <gpff files in RefSeq DB dir> --prefix <Output file>
```
