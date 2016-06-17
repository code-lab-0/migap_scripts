#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::DB::Fasta;

our $LINE_CHAR_FA = 80;
our $LINE_CHAR_GE = 79;

our $RUN_DATE = "dd-mm-yyyy";
our $COG_DB = "COG{20030417}";
our $TRE_DB = "TrEMBL{20140417}";
our $REF_DB = "RefSeq{20140911}";

my @cFiles = (); # Blast to cog result
my @tFiles = (); # Blast to trembl result
my @rFiles = (); # Blast to refseq result
#my $oFile  = ""; # taxonomy file
my $wFile  = ""; # COG whog file
my $mFile  = ""; # COG myva=gb file
my @gFiles = (); # RefSeq gpff.gz files
my $prefix = ""; # input/output prefix

GetOptions("cog=s{,}"     => \@cFiles,
           "trembl=s{,}"  => \@tFiles,
           "ref=s{,}"     => \@rFiles,
#          "org=s"        => \$oFile,
           "whog=s"       => \$wFile,
           "myva=s"       => \$mFile,
           "gpff=s{,}"    => \@gFiles,
           "prefix=s"     => \$prefix);

#our $ORG2TID = &readOFile($oFile);
our $ACC2COG = &readCOGFiles($wFile, $mFile);

my $data;
my $useRef;
foreach my $file(@cFiles){
  ($data, $useRef) = &readBlastP($data, $useRef, $file, $COG_DB, 60, 60);
}
foreach my $file(@rFiles){
  ($data, $useRef) = &readBlastP($data, $useRef, $file, $REF_DB, 60, 60);
}
foreach my $file(@tFiles){
  ($data, $useRef) = &readBlastP($data, $useRef, $file, $TRE_DB, 30, 30);
}

our $ACC2REFSEQ = &readGFiles($useRef, \@gFiles);

my $inGbkFile = $prefix . ".gbk";
my $inEmblFile = $prefix . ".embl";
my $inAnntFile = $prefix. ".annt";
my $inCsvFile = $prefix. ".csv";
my $inNtFaFile = $prefix . "-na.fasta";
my $inAaFaFile = $prefix . "-aa.fasta";

my $outGbkFile = $prefix . "-a.gbk";
my $outEmblFile = $prefix . "-a.embl";
my $outAnntFile = $prefix. "-a.annt";
my $outCsvFile = $prefix. "-a.csv";
my $outDdbjFile = $prefix. "-a.ddbj";

&addBlastGbk($data, $inGbkFile, $outGbkFile, 0);
&addBlastGbk($data, $inGbkFile, $outDdbjFile, 1);
&addBlastEmbl($data, $inEmblFile, $outEmblFile);
&addBlastAnnt($data, $inAnntFile, $outAnntFile);
&addBlastCsv($data, $inCsvFile, $outCsvFile);


#---------------------------------------------------
# functions

# read cog files
sub readCOGFiles {
  my ($wFile, $mFile) = @_;
  my %tmp;
  my %acc2cog;

  open my $fh, "<", $wFile or die;
  while(<$fh>){
    chomp;
    next if /^$/;
    if(/^___/){
      %tmp = ();
    }
    elsif(/^\[(\w+)\] \s*(COG\d+)\s+(\S.*\S)\s*$/){
      $tmp{"code"} = $1;
      $tmp{"cid"} = $2;
      $tmp{"func"} = $3;
    }
    elsif(/\S+:\s+(\S.+)$/){
      my @arr = split(/\s+/, $1);
      foreach my $acc(@arr){
        $acc2cog{$acc}{"code"} = $tmp{"code"};
        $acc2cog{$acc}{"cid"} = $tmp{"cid"};
        $acc2cog{$acc}{"func"} = $tmp{"func"};
      }
    }
  }
  close $fh;

  open $fh, "<", $mFile or die;
  while(<$fh>){
    chomp;
    next if /^$/;
    my ($acc, $gi) = split(/\s+/, $_);
    $acc2cog{$acc}{"gi"} = $gi;
  }
  close $fh;

  return \%acc2cog;
}

# org to taxonomy id
sub readOFile {
  my $file = shift @_;
  my %org2tid;

  open my $fh, "<", $file or die;
  while(<$fh>){
    chomp;
    next if /^$/ || /^locus/;
    my @arr = split(/\t/, $_);
    $org2tid{$arr[2]} = $arr[1];
  }
  close $fh;

  return \%org2tid;
}

sub readGFiles {
  my ($useRef, $files) = @_;
  my $acc2refseq;
  foreach my $file(@$files){
    $acc2refseq = &readGFile($acc2refseq, $useRef, $file);
  }
  return $acc2refseq;
}

sub readGFile {
  my ($acc2refseq, $useRef, $file) = @_;
  my %col = ("source"=>1, "Protein"=>1, "Region"=>1,
             "Site"=>1, "CDS"=>1);
  my %nodup;
  my %tmp = ();
  open my $fin, "zcat $file 2>/dev/null |" or die;
  while(<$fin>){
    chomp;
    if(/^VERSION\s+(\S+)/){
      if(exists $$useRef{$1}){
        $tmp{"acc"} = $1;
      }
    }
    elsif(exists $tmp{"acc"} && /^\s{5}(\S+)\s{10,}(.+)$/){
      my ($c, $v) = ($1, $2);
      if($c eq "Site" && $v !~ /\)/){
        $tmp{"write"} = -1;
      }
      elsif(exists $col{$c}){
        $tmp{"write"} = 1;
      }
    }
#    elsif(exists $tmp{"write"} && /^CONTIG/){
    elsif(exists $tmp{"write"} && /^ORIGIN/){    # CONTIGの存在しないエントリーがあるのでORIGINに変更
      $$acc2refseq{ $tmp{"acc"} } = $tmp{"arr"};
      %tmp = ();
      %nodup = ();
    }
    elsif(exists $tmp{"write"} && /^\s{21}(\S.+)/){
      my $v = $1;
      if($tmp{"write"} == -1 && $v =~ /\)/){
        $tmp{"write"} = 1;
      }
      elsif($tmp{"write"} == 1){
        my $tmparr = exists $tmp{"arr"} ? $tmp{"arr"} : [];
        if($v =~ /^\/transl_table/ || $v =~ /^\/locus_tag/){
          next;
        }
        elsif($v =~ /^\/db_xref/ || $v =~ /^\/site_type/
          || $v =~ /^\/note/){
          next if exists $nodup{$v};
          $nodup{$v} = 1;
        }

        if($v =~ /^\//){
          push @$tmparr, $v;
        }
        else{
          my $prev = pop @$tmparr;
          $prev .= " $v";
          push @$tmparr, $prev;
        }
        $tmp{"arr"} = $tmparr;
      }
    }
  }
  close $fin;
  return $acc2refseq;
}

# parse blastp result
sub readBlastP {
  my ($data, $useRef, $file, $dbname, $idenTh, $ovTh) = @_;

  my %seqinfo = ();
  my %info = ();
  my $ind = 0;

  open my $fin, '<', $file or die;
  while(<$fin>){
    chomp;

    if(/^Query= (\S+)/){
      $seqinfo{"que"} = $1;
    }
    elsif(/^>\s*(\S.+)$/){
      $seqinfo{"tar"} = $1;
      $seqinfo{"inTar"} = 1;
    }
    elsif(/^Length=(\d+)/ && exists $seqinfo{"inTar"}){
      $seqinfo{"len"} = $1;
      delete $seqinfo{"inTar"}
    }
    elsif(exists $seqinfo{"inTar"}){
      $seqinfo{"tar"} .= $_;
    }
    elsif(/^ Identities = (\d+)\/(\d+).+ Gaps = (\d+)\/\d+/){
      $info{"iden"} = $1;
      $info{"alen"} = $2;
      $info{"gap"} = $3;
    }
    elsif(/^(Query\s+(\d+)\s+)(\S+)\s+(\d+)/){
      if(! exists $info{"qsta"}){
        $info{"qsta"} = $2;
        $info{"qseq"} = "";
        $info{"mid"} = "";
        $info{"sseq"} = "";
        $ind = length($1);
      }
      $info{"qseq"} .= $3;
      $info{"qend"} = $4;
      $info{"inAln"} = 1;
    }
    elsif(/^\s{$ind}(.+)$/ && exists $info{"inAln"} && $info{"inAln"} == 1){
      $info{"mid"} .= $1;
    }
    elsif(/^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/){
      if(! exists $info{"ssta"}){
        $info{"ssta"} = $1;
      }
      $info{"sseq"} .= $2;
      $info{"send"} = $3;
      $info{"inAln"} = 0;
    }
    elsif(/^ Score =\s*(\d+\.?\d*) bits.+Expect = (\S+),/){
      if(exists $info{"sco"}){
        my $ov = &calcOv($info{"alen"}, length($info{"qseq"}),
                         length($info{"sseq"}));
        if($info{"iden"} >= $idenTh && $ov >= $ovTh){
          $info{"str"} = "+";
          ($data, $useRef) = &putBlast($data, $useRef,
                             \%seqinfo, \%info, $dbname);
        }
      }
      %info = ();
      $info{"sco"} = sprintf("%.1f", $1);
      $info{"eval"} = $2;
    }
    elsif(/^Lambda/){
      if(exists $info{"sco"}){
        my $ov = &calcOv($info{"alen"}, length($info{"qseq"}),
                         length($info{"sseq"}));
        if($info{"iden"} >= $idenTh && $ov >= $ovTh){
          $info{"str"} = "+";
          ($data, $useRef) = &putBlast($data, $useRef,
                             \%seqinfo, \%info, $dbname);
        }
      }
      %info = ();
    }

  }
  close $fin;

  return $data, $useRef;
}

sub calcOv {
  my ($alen, $qlen, $slen) = @_;
  my $longer = $qlen <= $slen ? $slen : $qlen;
  my $ov = int($alen * 100 / $longer);
  return $ov; 
}

# put blast data to hash
sub putBlast {
  my ($data, $useRef, $seqinfo, $info, $dbname) = @_;

  my ($que) = ($$seqinfo{"que"});
  my ($tar, $len) = ($$seqinfo{"tar"}, $$seqinfo{"len"});
  my ($sta, $end, $strand) = ($$info{"qsta"}, $$info{"qend"}, $$info{"str"});

  # when it is annotated by cog or/and refseq, skip trembl
  if(exists $$data{$que} && $dbname =~ /^TrEMBL/i){
    return $data, $useRef;
  }

  my $tmp = exists $$data{$que} ? $$data{$que} : [];

  my $blast_type = "AA";

  # blast:
    # AA|NA, db, target, length, identity, alignment length,
    # gap, gapopen, qsta, qend, ssta, send
    # eval, bitsco, qseq, hseq, midline

  if($dbname =~ /^COG/i && $tar =~ /^\S+\|(\S+) /){
    $tar = $1;
  }

  my $test = $$info{"sseq"} . $$info{"qseq"};
  $$info{"gop"} = 0;
  $$info{"gop"} = ($test =~ s/\-+//g) if $test =~ /\-/;

  my @blast = ($blast_type, $dbname, $tar, $len,
               $$info{"iden"}, $$info{"alen"},
               $$info{"gap"}, $$info{"gop"},
               $$info{"qsta"}, $$info{"qend"},
               $$info{"ssta"}, $$info{"send"},
               $$info{"eval"}, $$info{"sco"},
               $$info{"qseq"}, $$info{"sseq"}, $$info{"mid"}
              );

  my @arr = ( "blast=" . join("<>", @blast) );

  if($dbname =~ /^COG/i){
    if(exists $$ACC2COG{$tar}{"code"}){
      push @arr, "classification=\"" . $$ACC2COG{$tar}{"code"} . "\"";
      push @arr, "classification=\"" . $$ACC2COG{$tar}{"cid"} . "\"";
      push @arr, "function=\"" . $$ACC2COG{$tar}{"func"} . "\"";
    }
    if(exists $$ACC2COG{$tar}{"gi"}){
      push @arr, "db_xref=\"GI:" . $$ACC2COG{$tar}{"gi"} . "\"";
    }
  }
  elsif($dbname =~ /^TrEMBL/i){
    if($tar =~ /^tr/){
      my @tmparr = split(/\|/, $tar);
      push @arr, "db_xref=\"" . $tmparr[1] . "\"";
      if($tmparr[2] =~ /^(.+) OS=(.+) GN=(\S+) PE=(\d+) SV=(\d+)/){
        push @arr, "product=\" . $1 . \"";
        push @arr, "note=\"organism:" . $2 . "\"";
        push @arr, "gene=\"" . $3 . "\"";
        push @arr, "note=\"protein existence:" . $4 . "\"";
        push @arr, "note=\"sequence version:" . $5 . "\"";
      }
    }
  }
  elsif($dbname =~ /^RefSeq/i){
    if($tar =~ /^gi\|(\S+)\|ref\|(\S+)\|\s*((\S.*)\s*\[(.+)\])\|/){
      my ($gi, $ref, $nt, $pr, $org) = ($1, $2, $3, $4, $5);
      push @arr, "db_xref=\"GI:$gi\"";
      push @arr, "db_xref=\"ref:$ref\"";
      $$useRef{$ref} = 1;
      push @arr, "note=\"$nt\"";

#     push @arr, "organism=\"$org\"";
#     if($org =~ / sp\.\s+(\S.+)/){
#       push @arr, "strain=\"" . $1 . "\"";
#     }
#     if(exists $$ORG2TID{$org}){
#       push @arr, "db_xref=\"taxon:" . $$ORG2TID{$org} . "\"";
#     }
#     push @arr, "product=\"" . $pr . "\"";
    }
  }

  push @$tmp, \@arr;
  $$data{$que} = $tmp;

  return $data, $useRef;
}

# insert break into str
sub breakStr {
  my ($str, $ind, $line_char) = @_;
  my $parsed = "";
  my $limit = $line_char - length($ind);

  while(length($str) > $limit){
    $parsed .= sprintf("\n%s%s", $ind, substr($str, 0, $limit));
    $str = substr($str, $limit);
  }
  $parsed .= "\n$ind$str" if $str ne "";
  $parsed =~ s/^\n//;

  return $parsed;
}

# add blast info to Genbank
sub addBlastGbk {
  my ($data, $inFile, $outFile, $isDdbj) = @_;

  my $ind = " " x 21;

  my %tmp = ();
  my $not_basecount = 1;   # ORIGIN中にannotationを埋め込ませない。
  open my $fout, ">", $outFile or die;
  open my $fin, "<", $inFile or die;
  while(<$fin>){
    chomp;
    if(/^\s{5}(\S+)\s+(\S+)/ && $not_basecount){
      if(exists $tmp{"isCDS"} && $tmp{"isCDS"} == 1 &&
         exists $$data{$tmp{"pos"}}){
        my $alnarr = $$data{$tmp{"pos"}};
        for my $aln(@$alnarr){
          my $acc = "";
          for my $a(@$aln){
            print $fout &breakStr(sprintf("/%s", $a),
              $ind, $LINE_CHAR_GE), "\n";
            $acc = $1 if $a =~ /db_xref="ref:(\S+)"/;
          }
          if(exists $$ACC2REFSEQ{$acc}){
            my $tmpref = $$ACC2REFSEQ{$acc};
            for my $a(@$tmpref){
              print $fout &breakStr(sprintf("%s", $a),
                $ind, $LINE_CHAR_GE), "\n";
            }
          }
        }
      }
      %tmp = ("name"=>$1, "pos"=>$2, "isCDS"=>0);
      $tmp{"isCDS"} = 1 if $tmp{"name"} eq "CDS";
    }
    elsif(/^BASE COUNT/){
      $not_basecount = 0;
      if(exists $tmp{"isCDS"} && $tmp{"isCDS"} == 1 &&
         exists $$data{$tmp{"pos"}}){
        my $alnarr = $$data{$tmp{"pos"}};
        for my $aln(@$alnarr){
          my $acc = "";
          for my $a(@$aln){
            print $fout &breakStr(sprintf("/%s", $a),
              $ind, $LINE_CHAR_GE), "\n";
            $acc = $1 if $a =~ /db_xref="ref:(\S+)"/;
          }
          if(exists $$ACC2REFSEQ{$acc}){
            my $tmpref = $$ACC2REFSEQ{$acc};
            for my $a(@$tmpref){
              print $fout &breakStr(sprintf("%s", $a),
                $ind, $LINE_CHAR_GE), "\n";
            }
          }
        }
      }
    }
    print $fout $_, "\n";
  }
  close $fin;
  close $fout;
}

# add blast info to Embl
sub addBlastEmbl {
  my ($data, $inFile, $outFile) = @_;

  my $ind = "FT" . " " x 19;

  my %tmp = ();
  open my $fout, ">", $outFile or die;
  open my $fin, "<", $inFile or die;
  while(<$fin>){
    chomp;
    if(/^FT\s{3}(\S+)\s+(\S+)/){
      if(exists $tmp{"isCDS"} && $tmp{"isCDS"} == 1 &&
         exists $$data{$tmp{"pos"}}){
        my $alnarr = $$data{$tmp{"pos"}};
        for my $aln(@$alnarr){
          my $acc = "";
          for my $a(@$aln){
            print $fout &breakStr(sprintf("/%s", $a),
              $ind, $LINE_CHAR_GE), "\n";
            $acc = $1 if $a =~ /db_xref="ref:(\S+)"/;
          }
          if(exists $$ACC2REFSEQ{$acc}){
            my $tmpref = $$ACC2REFSEQ{$acc};
            for my $a(@$tmpref){
              print $fout &breakStr(sprintf("%s", $a),
                $ind, $LINE_CHAR_GE), "\n";
            }
          }
        }
      }
      %tmp = ("name"=>$1, "pos"=>$2, "isCDS"=>0);
      $tmp{"isCDS"} = 1 if $tmp{"name"} eq "CDS";
    }
    elsif(/^XX/){
      if(exists $tmp{"isCDS"} && $tmp{"isCDS"} == 1 &&
         exists $$data{$tmp{"pos"}}){
        my $alnarr = $$data{$tmp{"pos"}};
        for my $aln(@$alnarr){
          for my $a(@$aln){
            print $fout &breakStr(sprintf("/%s", $a),
              $ind, $LINE_CHAR_GE), "\n";
          }
        }
      }
    }
    print $fout $_, "\n";
  }
  close $fin;
  close $fout;
}

# add blast info to Annt
sub addBlastAnnt {
  my ($data, $inFile, $outFile) = @_;

  my %tmp = ();
  open my $fout, ">", $outFile or die;
  open my $fin, "<", $inFile or die;
  while(<$fin>){
    chomp;
    if(/^\t(\S+)\t(\S+)/){
      if(exists $tmp{"isCDS"} && $tmp{"isCDS"} == 1 &&
         exists $$data{$tmp{"pos"}}){
        my $alnarr = $$data{$tmp{"pos"}};
        for my $aln(@$alnarr){
          my $newaln = &gbkToAnnt($aln);
          for my $a(@$newaln){
            print $fout "\t\t\t$a\n";
          }
        }
      }
      %tmp = ("name"=>$1, "pos"=>$2, "isCDS"=>0);
      $tmp{"isCDS"} = 1 if $tmp{"name"} eq "CDS";
    }
    print $fout $_, "\n";
  }
  close $fin;
  close $fout;
}

sub gbkToAnnt {
  my $aln = shift @_;

  my $acc = "";
  my @merge = ();

  for my $a(@$aln){
    push @merge, $a;
    $acc = $1 if $a =~ /db_xref="ref:(\S+)"/;
  }
  if(exists $$ACC2REFSEQ{$acc}){
    my $tmpref = $$ACC2REFSEQ{$acc};
    for my $a(@$tmpref){
      push @merge, $a;
    }
  }

  my $j = -1;
  my @arr;
  for(my $i=0; $i<@merge; $i++){
    if($merge[$i] =~ /^\S+=/){
      $j ++;
      $merge[$i] =~ s/^\///;
      $merge[$i] =~ s/=/\t/;
      $arr[$j] = $merge[$i];
    }
    else{
      $arr[$j] .= $merge[$i];
    }
  }

  my @res;
  for $a(@arr){
    if($a =~ /^function/ || $a =~ /^db_xref/ ||
       $a =~ /^note/ || $a =~ /^product/){
      push @res, $a;
    }
  }

  return \@res;
}

sub addBlastCsv {
  my ($data, $inFile, $outFile) = @_;

  my @header = qw(Feature Location /anticodon= /artificial_location=); # 1-4
  push @header, ("/blast=") x 2; # 5-6
  push @header, ("/calculated_mol_wt=") x 2; # 7-8
  push @header, "/chromosome="; # 9
  push @header, ("/classification=") x 2; # 10-11
  push @header, qw(/coded_by= /codon_start=); # 12-13
  push @header, qw(/collected_by= /collection_date= /country=); # 14-16
  push @header, ("/culture_collection=") x 2; # 17-18
  push @header, ("/db_xref=") x 21; # 19-39
  push @header, ("/EC_number=") x 2; # 40-41
  push @header, ("/function=") x 2; # 42-43
  push @header, qw(/gene= /isolation_source= /lat_lon=); # 44-46
  push @header, qw(/migap_code= /name= /new=); # 47-49
  push @header, ("/note=") x 31; # 50-80
  push @header, qw(/organism= /plasmid=); # 81-82
  push @header, ("/product=") x 2; # 83-84
  push @header, ("/region_name=") x 15; # 85-99
  push @header, "/serotype="; # 100
  push @header, ("/site_type=") x 4; # 101-104
  push @header, qw(/strain= /sub_strain= /synonym=); # 105-107
  push @header, qw(/transl_table= /translation= Sequence); # 108-110

  open my $fout, ">", $outFile or die;
  print $fout join(",", @header), "\n";

  my $header_edited = &editHeader(\@header);

  open my $fin, "<", $inFile or die;
  while(<$fin>){
    chomp;
    next if /^Feature/;

    # escape ',' in tRNA data
    if(/^tRNA/){
      s/,aa:/___aa:/;
      s/,seq:/___seq:/;
    }
    my @in = split(/,/, $_, -1);
    $in[2] =~ s/___/,/g;

    my @out = ("") x 110;
    @out[0, 1, 2] = @in[0, 1, 2];
    @out[12, 49, 50, 82, 107, 108, 109] = @in[4, 7, 8, 9, 10, 11, 12];

    my $pos = $in[1];
    if(exists $$data{$pos}){
      my $alnarr = $$data{$pos};
      for my $aln(@$alnarr){
        my $acc = "";
        for my $a(@$aln){
          if($a =~ /db_xref="ref:(\S+)"/){
            $acc = $1;
          }
        }
        if(exists $$ACC2REFSEQ{$acc}){
          my $tmpref = $$ACC2REFSEQ{$acc};
          push @$aln, @$tmpref;
        }

        for my $a(@$aln){
          if($a =~ /^\/?(.+)=(.+)$/){
            my ($k, $v) = ($1, $2);
            my $ind = &findInd($header_edited, \@out, $k);
            $out[$ind] = $v if $ind > -1;
          }
        }
      }
    }
    print $fout join(",", @out), "\n";
  }
  close $fin;
  close $fout;
}

sub editHeader {
  my $header = shift @_;
  my $res = [];
  for my $h(@$header){
    $h =~ s/^\///;
    $h =~ s/=\s*$//;
    push @$res, $h;
  }
  return $res;
}

sub findInd {
  my ($h, $o, $k) = @_;

  my $i;
  for($i=0; $i<@$h; $i++){
    if($$h[$i] eq $k){
      if($$o[$i] eq "" || $$h[$i+1] ne $k){
        return $i;
      }
    }
  }
  return -1;
}

