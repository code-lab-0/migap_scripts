#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::DB::Fasta;

our $LINE_CHAR_FA = 80;
our $LINE_CHAR_GE = 79;

our $RUN_DATE = "dd-mm-yyyy";
our $RRNA_DB = "16S rRNA{20090220}";

our $CSV_HEADER = "Feature,Location,/anticodon=,/blast=,/codon_start=,"
. "/migap_code=,/new=,/note=,/note=,/product=,/transl_table=,"
. "/translation=,Sequence\n";

our $GFF_HEADER =<< "EOF";
##gff-version 3
##date $RUN_DATE
EOF

our $TRANSL_TABLE = "11";
our $CODON_START = "1";

my $fFile = "";  # fasta
my @mFiles = (); # mga result
my @tFiles = (); # tRNAscan-SE result
my @rFiles = (); # RNAmmer result
my @bFiles = (); # Blast to 16srna result
my $prefix = ""; # output prefix

GetOptions("fa=s"          => \$fFile,
           "mga=s{,}"      => \@mFiles,
           "trnascan=s{,}" => \@tFiles,
           "rnammer=s{,}"  => \@rFiles,
           "blast16s=s{,}" => \@bFiles,
           "prefix=s"      => \$prefix);

my $data;
$data = &readMFiles($data, \@mFiles);
$data = &readTFiles($data, \@tFiles);
$data = &readRFiles($data, \@rFiles, \@bFiles);
$data = &sortByPos($data);

my $SEQ_TYPE = "DNA circular";
my $SEQ_NOTE = "molecular_form:circular";
my $info = &readFFile($fFile);

my $outGbkFile = $prefix . ".gbk";
my $outEmblFile = $prefix . ".embl";
my $outNtFaFile = $prefix . "-na.fasta";
my $outAaFaFile = $prefix . "-aa.fasta";
my $outAnntFile = $prefix. ".annt";
my $outCsvFile = $prefix. ".csv";
my $outGffFile = $prefix. ".gff";

&writeGbkEmbl($fFile, $outGbkFile, $outEmblFile, $data, $info);
&writeFa($fFile, $outNtFaFile, $outAaFaFile, $data, $info);
&writeAnnt($fFile, $outAnntFile, $data, $info);
&writeCsvGff($fFile, $outCsvFile, $outGffFile, $data, $info);


#---------------------------------------------------
# functions

# return nt or aa sequence from fasta
sub getSeq {
  my ($fa, $seq, $sta, $end, $strand, $type) = @_;

  my $db = Bio::DB::Fasta->new($fa);

  my $obj = Bio::Seq->new( -alphabet => "dna",
            -seq => $db->seq($seq, $sta, $end) );

  my $str = "";

  if($strand eq "+"){
    if($type eq "a"){
      $str = $obj->translate(-codontable_id=>$TRANSL_TABLE)->seq;
    }
    else{
      $str = $obj->seq;
    }
  }
  else{
    if($type eq "a"){
      $str = $obj->revcom->translate(-codontable_id=>$TRANSL_TABLE)->seq;
    }
    else{
      $str = $obj->revcom->seq;
    }
  }

  $str =~ s/\*$//;
  return $str;
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

# read original fasta file
sub readFFile {
  my ($file) = @_;
  my $info;
  my $seq = "";
  my %stat = ("bp"=>0, "a"=>0, "c"=>0, "g"=>0, "t"=>0);

  open my $fin, '<', $file or die;
  while(<$fin>){
    s/\r$//;
    chomp;
    next if /^$/;

    if(/^>(\S+)/){
      $$info{$seq} = &addInfo(\%stat) if $seq ne "";
      $seq = $1;
      %stat = ("bp"=>0, "a"=>0, "c"=>0, "g"=>0, "t"=>0);
    }
    else{
      $stat{"bp"} += length($_);
      map { $stat{$_}++ } split(//, lc($_));
    }
  }
  close $fin;

  $$info{$seq} = &addInfo(\%stat) if $seq ne "";
  return $info;
}

# push data into hash
sub addInfo {
  my ($stat) = @_;
  my $hash = {
    "bp" => $$stat{"bp"},
    "a" => $$stat{"a"},
    "c" => $$stat{"c"},
    "g" => $$stat{"g"},
    "t" => $$stat{"t"}
  };
  return $hash;
}

# read mga files
sub readMFiles {
  my ($data, $files) = @_;
  foreach my $file(@$files){
    $data = &readMFile($data, $file);
  }
  return $data;
}

# read mga file
sub readMFile {
  my ($data, $file) = @_;
  my $status = 0;  # 1: header, 2: data
  my $seq = "";

  open my $fin, '<', $file or die;
  while(<$fin>){
    chomp;

    if(/^# (\S+)/ && $status == 0){
      $status = 1;
      $seq = $1;
    }
    elsif(/^\w/ && $status == 1){
      $status = 2;
    }

    if($status == 2){
      next if /^$/;
      my @arr = split(/\t/, $_);
      my $tmp = exists $$data{$seq} ? $$data{$seq} : [];
      my ($id, $strand, $frame, $sco) = @arr[0,3,4,6];
      my ($sta, $end) = @arr[1..2];

      push @$tmp, {"name"=>"CDS", "sta"=>$sta, "end"=>$end, 
        "tool"=>"MetaGeneAnnotator", "id"=>$id, "frame"=>$frame,
        "sco"=>$sco, "strand"=>$strand};

      my ($rsta, $rend, $rsco) = @arr[8..10];
      if($rsta ne "-"){
        push @$tmp, {"name"=>"RBS", "sta"=>$rsta, "end"=>$rend, 
          "sco"=>$rsco, "tool"=>"MetaGeneAnnotator", "strand"=>$strand};
      }

      $$data{$seq} = $tmp;
    }
  }
  close $fin;
  return $data;
}

# read trnascan-seq files
sub readTFiles {
  my ($data, $files) = @_;
  foreach my $file(@$files){
    $data = &readTFile($data, $file);
  }
  return $data;
}

# read trnascan-seq file
sub readTFile {
  my ($data, $file) = @_;

  my ($seq, $sta, $end, $type, $codon, $strand, $csta, $cend) = ();

  open my $fin, '<', $file or die;
  while(<$fin>){
    chomp;

    if(/^Seq:/ || /^Str:/){
      next;
    }
    elsif(/^(\S+)\.trna\w+\s*\((\d+)-(\d+)\)/){
      ($seq, $sta, $end, $strand) = ($1, $2, $3, "+");
      ($sta, $end, $strand) = ($3, $2, "-") if $sta > $end;
    }
    elsif(/^Type: (\S+)\s+Anticodon: (\w{3}|\?\?\?) at .+\((\d+)\-(\d+)\)/){
      ($type, $codon, $csta, $cend) = ($1, $2, $3, $4);
      ($csta, $cend) = ($4, $3) if $csta > $cend;
    }
    elsif(/^$/){
      my $tmp = exists $$data{$seq} ? $$data{$seq} : [];

      push @$tmp, {"name"=>"tRNA", "sta"=>$sta, "end"=>$end,
        "type"=>$type, "codon"=>$codon, "tool"=>"tRNAscan-SE",
        "strand"=>$strand, "csta"=>$csta, "cend"=>$cend};

      $$data{$seq} = $tmp;

      ($seq, $sta, $end, $type, $codon, $strand, $csta, $cend) = ();
    }
  }
  close $fin;
  return $data;
}

# read rnammer files
sub readRFiles {
  my ($data, $rfiles, $bfiles) = @_;

  my $rdata;
  foreach my $rfile(@$rfiles){
    $rdata = &readRFile($rdata, $rfile);
  }

  my $bdata;
  foreach my $bfile(@$bfiles){
    $bdata = &readBFile($bdata, $bfile);
  }

  while(my ($seq, $rtmp) = each %$rdata){
    my $tmp = exists $$data{$seq} ? $$data{$seq} : [];
    foreach my $rf( @$rtmp ){
      my $ovlp = 0;
      if(exists $$bdata{$seq}){
        my $btmp = $$bdata{$seq};
        foreach my $bf( @$btmp ){
          if($$bf{"sta"} <= $$rf{"sta"} && $$rf{"end"} <= $$bf{"end"}){
            $ovlp = 1;
            next;
          }
        }
      }
      push @$tmp, $rf if $ovlp == 0;
    }
    $$data{$seq} = $tmp;
  }

  while(my ($seq, $btmp) = each %$bdata){
    my $tmp = exists $$data{$seq} ? $$data{$seq} : [];
    foreach my $bf( @$btmp ){
      push @$tmp, $bf;
    }
    $$data{$seq} = $tmp;
  }

  return $data;
}

# read rnammer file
sub readRFile {
  my ($data, $file) = @_;
  open my $fin, '<', $file or die;
  while(<$fin>){
    chomp;
    next if /^#/ || /^$/;
    my @arr = split(/\t/, $_);

    my ($seq, $sta, $end, $strand, $type) = @arr[0,3,4,6,8];
    $type =~ s/s_rRNA/S ribosomal RNA/;

    my $tmp = exists $$data{$seq} ? $$data{$seq} : [];

    push @$tmp, {"name"=>"rRNA", "sta"=>$sta, "end"=>$end,
                 "type"=>$type, "strand"=>$strand, "tool"=>"RNAmmer"};

    $$data{$seq} = $tmp;
  }
  close $fin;
  return $data;
}

# read blast result
sub readBFile {
  my ($data, $file) = @_;
  $data = &readBlastN($data, $file, "rRNA", "ribosomal RNA-16S", $RRNA_DB);
  return $data;
}

# parse blastn result
sub readBlastN {
  my ($data, $file, $name, $type, $dbname) = @_;

  my %seqinfo = ();
  my %info = ();
  my $ind = 0;

  open my $fin, '<', $file or die;
  while(<$fin>){
    chomp;

    if(/^Query= (\S+)/){
      $seqinfo{"seq"} = $1;
    }
    elsif(/^>\s*(\S+)/){
      my @tmparr = split(/\|/, $1);
      $seqinfo{"tar"} = pop @tmparr;
    }
    elsif(/^Length=(\d+)/ && exists $seqinfo{"tar"}){
      $seqinfo{"len"} = $1;
    }
    elsif(/^ Identities = (\d+)\/(\d+).+ Gaps = (\d+)\/\d+/){
      $info{"iden"} = $1;
      $info{"alen"} = $2;
      $info{"gap"} = $3;
    }
    elsif(/^ Strand=(\w+)\/(\w+)/){
      $info{"str"} = $1 eq $2 ? "+" : "-";
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
    elsif(/^ Score =\s*(\d+\.?\d*) bits.+Expect = (\S+),?/){
      if(exists $info{"sco"}){
        $data = &putBlast($data, \%seqinfo, \%info, $name, $type, $dbname);
      }
      %info = ();
      $info{"sco"} = sprintf("%.1f", $1);
      $info{"eval"} = $2;
    }
    elsif(/^Lambda/){
      if(exists $info{"sco"}){
        $data = &putBlast($data, \%seqinfo, \%info, $name, $type, $dbname);
      }
      %info = ();
    }

  }
  close $fin;

  return $data;
}

# put blast data to hash
sub putBlast {
  my ($data, $seqinfo, $info, $name, $type, $dbname) = @_;

  my ($seq, $tar, $len) = ($$seqinfo{"seq"}, $$seqinfo{"tar"}, $$seqinfo{"len"});
  my ($sta, $end, $strand) = ($$info{"qsta"}, $$info{"qend"}, $$info{"str"});

  my $tmp = exists $$data{$seq} ? $$data{$seq} : [];

  my $blast_type = "NA";
  $$info{"qseq"} = lc($$info{"qseq"});
  $$info{"sseq"} = lc($$info{"sseq"});

  $$info{"qend"} -= ($$info{"qsta"} - 1);
  $$info{"qsta"} = 1;

  # blast:
    # AA|NA, db, target, length, identity, alignment length,
    # gap, gapopen, qsta, qend, ssta, send
    # eval, bitsco, qseq, hseq, midline

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

  push @$tmp, {"name"=>$name, "sta"=>$sta, "end"=>$end,
               "type"=>$type, "strand"=>$strand, "tool"=>"",
               "blast"=>join("<>", @blast)};

  $$data{$seq} = $tmp;
  return $data;
}

# sort hash by position
sub sortByPos {
  my $data = shift @_;
  my $sorted;

  while(my ($seq, $tmp) = each %$data){
    my @tmp_sorted = sort { $a->{"sta"} <=> $b->{"sta"} } @$tmp;
    $$sorted{$seq} = \@tmp_sorted;
    delete $$data{$seq};
  }
  return $sorted;
}

# write genome sequence to genbank and embl
sub writeGenomeSeq {
  my ($fa, $seq, $foutGbk, $foutEmbl, $type) = @_;

  my $indnum = 1;
  my $status = 0;
  my @nt = ();

  open my $ffa, '<', $fa or die;
  while(<$ffa>){
    s/\r$//;
    chomp;
    next if /^$/;

    if(/^>(\S+)/){
      $status = 0;
      @nt = ();
      if($1 eq $seq){
        $status = 1;
        next;
      }
    }
    next if $status == 0;

    push @nt, split(//, lc($_));

    while(@nt > 59){
      my @sub = ();
      for(my $i=0; $i<6; $i++){
        next if @nt < 1;
        my $num = @nt > 9 ? 10 : @nt;
        my @tmp = splice(@nt, 0, $num);
        push @sub, join("", @tmp);
      }

      print $foutGbk sprintf("%9s %s\n", $indnum, join(" ", @sub));
      print $foutEmbl " " x 6 . join(" ", @sub) . "\n";

      $indnum += 60;
    }
  }
  close $ffa;

  if(@nt > 0){
    my @sub = ();
    for(my $i=0; $i<6; $i++){
      next if @nt < 1;
      my $num = @nt > 9 ? 10 : @nt;
      my @tmp = splice(@nt, 0, $num);
      push @sub, join("", @tmp);
    }

    print $foutGbk sprintf("%9s %s\n", $indnum, join(" ", @sub));
    print $foutEmbl " " x 6 . join(" ", @sub) . "\n";

  }
}

# write genbank and embl
sub writeGbkEmbl {
  my ($fa, $outGbk, $outEmbl, $data, $info) = @_;

  my $indGbk = " " x 21;
  my $indEmbl = "FT" . " " x 19;

  open my $foutGbk, '>', $outGbk or die;
  open my $foutEmbl, '>', $outEmbl or die;

  while(my ($seq, $tmp) = each %$data){
    my $s = $$info{$seq};

    print $foutGbk "LOCUS" . " " x 7 . $seq . " " . $$s{"bp"} . " bp";
    print $foutGbk " " x 4 . $SEQ_TYPE . " " x 20 . $RUN_DATE . "\n";
    print $foutGbk "DEFINITION  $seq\n";
    print $foutGbk "FEATURES" . " " x 13 . "Location/Qualifiers\n";
    print $foutGbk " " x 5 . "source" . " " x 10;
    print $foutGbk $CODON_START . ".." . $$s{"bp"} . "\n";
    print $foutGbk $indGbk . "/note=\"" . $SEQ_NOTE . "\"\n";

    print $foutEmbl sprintf("ID   %-15s standard;  %s; ; %s BP.\n", $seq, $SEQ_TYPE, $$s{"bp"});
    print $foutEmbl "XX\n";
    print $foutEmbl "DE   $seq\n";
    print $foutEmbl sprintf("FH   %-15s Location/Qualifiers\n", "Key");
    print $foutEmbl sprintf("FT   %-15s %s..%s\n", "source", $CODON_START, $$s{"bp"});
    print $foutEmbl $indEmbl . "/note=\"$SEQ_NOTE\"\n";

    foreach my $f( @$tmp ){
      print $foutGbk sprintf("%4s %-15s ", "", $$f{"name"});
      print $foutEmbl sprintf("FT   %-15s ", $$f{"name"});

      if($$f{"strand"} eq "+"){
        print $foutGbk sprintf("%s..%s\n", $$f{"sta"}, $$f{"end"});
        print $foutEmbl sprintf("%s..%s\n", $$f{"sta"}, $$f{"end"});
      }
      else{
        print $foutGbk sprintf("complement(%s..%s)\n", $$f{"sta"}, $$f{"end"});
        print $foutEmbl sprintf("complement(%s..%s)\n", $$f{"sta"}, $$f{"end"});
      }

      if($$f{"name"} eq "CDS"){
        my $aa = &getSeq($fa, $seq, $$f{"sta"}, $$f{"end"}, $$f{"strand"}, "a");

        print $foutGbk $indGbk . "/note=\"" . $$f{"id"} . "\"\n";
        print $foutGbk $indGbk . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";
        print $foutGbk $indGbk . "/transl_table=" . $TRANSL_TABLE . "\n";
        print $foutGbk $indGbk . "/codon_start=" . $CODON_START . "\n";
        my $aaGbk = &breakStr(sprintf("/translation=\"%s\"", $aa), $indGbk, $LINE_CHAR_GE);
        print $foutGbk $aaGbk . "\n";

        print $foutEmbl $indEmbl . "/note=\"" . $$f{"id"} . "\"\n";
        print $foutEmbl $indEmbl . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";
        print $foutEmbl $indEmbl . "/transl_table=" . $TRANSL_TABLE . "\n";
        print $foutEmbl $indEmbl . "/codon_start=" . $CODON_START . "\n";
        my $aaEmbl = &breakStr(sprintf("/translation=\"%s\"", $aa), $indEmbl, $LINE_CHAR_GE);
        print $foutEmbl $aaEmbl . "\n";
      }
      elsif($$f{"name"} eq "RBS"){
        print $foutGbk $indGbk . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";
        print $foutEmbl $indEmbl . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";
      }
      elsif($$f{"name"} eq "rRNA"){
        if($$f{"tool"} ne ""){
          print $foutGbk $indGbk . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";
          print $foutEmbl $indEmbl . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";
        }
        print $foutGbk $indGbk . "/product=\"" . $$f{"type"} . "\"\n";
        print $foutEmbl $indEmbl . "/product=\"" . $$f{"type"} . "\"\n";
        if(exists $$f{"blast"}){
          my $blastGbk = &breakStr(sprintf("/blast=%s", $$f{"blast"}), $indGbk, $LINE_CHAR_GE);
          print $foutGbk $blastGbk . "\n";
          my $blastEmbl = &breakStr(sprintf("/blast=%s", $$f{"blast"}), $indEmbl, $LINE_CHAR_GE);
          print $foutEmbl $blastEmbl . "\n";
        }
      }
      elsif($$f{"name"} eq "tRNA"){
        my $codonpos;
        if($$f{"strand"} eq "+"){
          $codonpos .= sprintf("%d..%d", $$f{"csta"}, $$f{"cend"});
        }
        else{
          $codonpos .= sprintf("complement(%d..%d)", $$f{"csta"}, $$f{"cend"});
        }

        my $anticodon = "/anticodon=\"(pos:" . $codonpos . ","
                      . "aa:". $$f{"type"} .",seq:" . $$f{"codon"} . ")\"";

        my $anticodonGbk = &breakStr($anticodon, $indGbk, $LINE_CHAR_GE);
        print $foutGbk $anticodonGbk . "\n";
        print $foutGbk $indGbk . "/product=\"tRNA-" . $$f{"type"} . "\"\n";
        print $foutGbk $indGbk . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";

        my $anticodonEmbl = &breakStr($anticodon, $indEmbl, $LINE_CHAR_GE);
        print $foutEmbl $anticodonEmbl . "\n";
        print $foutEmbl $indEmbl . "/product=\"tRNA-" . $$f{"type"} . "\"\n";
        print $foutEmbl $indEmbl . "/note=\"identified by " . $$f{"tool"} . "; putative\"\n";
      }
    }

    print $foutGbk "BASE COUNT  ";
    print $foutGbk sprintf("%d a %d c %d g %d t\n", $$s{"a"}, $$s{"c"}, $$s{"g"}, $$s{"t"});

    print $foutEmbl "XX\n";
    print $foutEmbl sprintf("SQ   Sequence %d BP; %d A; %d C; %d G; %d T; 0 other;\n",
      $$s{"bp"}, $$s{"a"}, $$s{"c"}, $$s{"g"}, $$s{"t"});

    print $foutGbk "ORIGIN\n";
    &writeGenomeSeq($fa, $seq, $foutGbk, $foutEmbl);

    print $foutGbk "//\n";
    print $foutEmbl "//\n";
  }
  close $foutGbk;
  close $foutEmbl;
}

# write nt and aa fasta
sub writeFa {
  my ($fa, $outNtFa, $outAaFa, $data, $info) = @_;

  open my $foutNt, '>', $outNtFa or die;
  open my $foutAa, '>', $outAaFa or die;

  while(my ($seq, $tmp) = each %$data){
    my $s = $$info{$seq};

    foreach my $f( @$tmp ){
      next if $$f{"name"} ne "CDS";

      if($$f{"strand"} eq "+"){
        print $foutNt sprintf(">%s..%s\n", $$f{"sta"}, $$f{"end"});
        print $foutAa sprintf(">%s..%s\n", $$f{"sta"}, $$f{"end"});
      }
      else{
        print $foutNt sprintf(">complement(%s..%s)\n", $$f{"sta"}, $$f{"end"});
        print $foutAa sprintf(">complement(%s..%s)\n", $$f{"sta"}, $$f{"end"});
      }

      my $nt = &getSeq($fa, $seq, $$f{"sta"}, $$f{"end"}, $$f{"strand"}, "n");
      print $foutNt &breakStr($nt, "", $LINE_CHAR_FA), "\n";

      my $aa = &getSeq($fa, $seq, $$f{"sta"}, $$f{"end"}, $$f{"strand"}, "a");
      print $foutAa &breakStr($aa, "", $LINE_CHAR_FA), "\n";
    }
  }

  close $foutNt;
  close $foutAa;
}

# write annt
sub writeAnnt {
  my ($fa, $outAnnt, $data, $info) = @_;
  my $entryid = 0;

  open my $foutAnnt, '>', $outAnnt or die;

  while(my ($seq, $tmp) = each %$data){
    my $s = $$info{$seq};

    $entryid ++;
    print $foutAnnt join("\t", "entry" . $entryid, "source",
      $CODON_START . ".." . $$s{"bp"}), "\n";

    print $foutAnnt sprintf("\t\t\tnote\t\"%s\"\n", $SEQ_NOTE);

    foreach my $f( @$tmp ){

      if($$f{"strand"} eq "+"){
        print $foutAnnt sprintf("\t%s\t%s..%s",
          $$f{"name"}, $$f{"sta"}, $$f{"end"});
      }
      else{
        print $foutAnnt sprintf("\t%s\tcomplement(%s..%s)",
          $$f{"name"}, $$f{"sta"}, $$f{"end"});
      }

      if($$f{"name"} eq "CDS"){
        print $foutAnnt "\tnote\t\"", $$f{"id"}, "\"\n";
        print $foutAnnt "\t\t\tnote\t\"identified by ", $$f{"tool"}, "; putative\"\n";
        print $foutAnnt "\t\t\ttransl_table\t", $TRANSL_TABLE, "\n";
        print $foutAnnt "\t\t\tcodon_start\t", $CODON_START, "\n";

        my $aa = &getSeq($fa, $seq, $$f{"sta"}, $$f{"end"}, $$f{"strand"}, "a");
        print $foutAnnt "\t\t\ttranslation\t\"", $aa, "\"\n";

        if(0){  # break each 80 chars
          my $aa_prefix = "";
          $aa = "\"" . $aa . "\"";
          print $foutAnnt "\t\t\ttranslation\t";
          while(length($aa) > 0){
            my $subLen = length($aa) >= 80 ? 80 : length($aa);
            print $foutAnnt $aa_prefix, substr($aa, 0, $subLen);
            $aa = substr($aa, $subLen);
            $aa_prefix = "\n";
          }
          print $foutAnnt "\n";
        }

      }
      elsif($$f{"name"} eq "RBS"){
        print $foutAnnt "\tnote\t\"identified by ", $$f{"tool"}, "; putative\"\n";
      }
      elsif($$f{"name"} eq "rRNA"){
        if($$f{"tool"} ne ""){
          print $foutAnnt "\tnote\t\"identified by ", $$f{"tool"}, "; putative\"\n\t\t";
        }
        print $foutAnnt "\tproduct\t\"", $$f{"type"}, "\"\n";
      }
      elsif($$f{"name"} eq "tRNA"){
        my $codonpos;
        if($$f{"strand"} eq "+"){
          $codonpos .= sprintf("%d..%d", $$f{"csta"}, $$f{"cend"});
        }
        else{
          $codonpos .= sprintf("complement(%d..%d)", $$f{"csta"}, $$f{"cend"});
        }

        print $foutAnnt "\tanticodon\t\"(pos:", $codonpos, ",";
        print $foutAnnt "aa:", $$f{"type"}, ",seq:", $$f{"codon"}, ")\"\n";
        print $foutAnnt "\t\t\tproduct\t\"tRNA-", $$f{"type"}, "\"\n";
        print $foutAnnt "\t\t\tnote\t\"identified by ", $$f{"tool"}, "; putative\"\n";
      }
    }
  }
  close $foutAnnt;
}

# write csv
sub writeCsvGff {
  my ($fa, $outCsv, $outGff, $data, $info) = @_;

  open my $foutCsv, '>', $outCsv or die;
  print $foutCsv $CSV_HEADER;

  open my $foutGff, '>', $outGff or die;
  print $foutGff $GFF_HEADER;

  while(my ($seq, $tmp) = each %$data){
    my $s = $$info{$seq};

    my @arr = ("") x 13;
    @arr[0, 1, 7] = ("source", join("..", $CODON_START, $$s{"bp"}), '"' . $SEQ_NOTE . '"');
    print $foutCsv join(",", @arr), "\n";

    my @gff = (".") x 9;
    @gff[0, 2, 3, 4, 6, 8] = ($seq, "source", $CODON_START, $$s{"bp"},
                              "+", "note=" . $SEQ_NOTE);
    print $foutGff join("\t", @gff), "\n";

    foreach my $f( @$tmp ){

      @arr = ("") x 13;
      $arr[0] = $$f{"name"};

      @gff = (".") x 9;
      @gff[0, 2, 3, 4, 6] = ($seq, $$f{"name"},
                             $$f{"sta"}, $$f{"end"}, $$f{"strand"});
      $gff[1] = $$f{"tool"} if $$f{"tool"} ne "";

      if($$f{"strand"} eq "+"){
        $arr[1] = sprintf("%s..%s", $$f{"sta"}, $$f{"end"});
      }
      else{
        $arr[1] = sprintf("complement(%s..%s)", $$f{"sta"}, $$f{"end"});
      }

      if($$f{"name"} eq "CDS"){
        $arr[4] = $CODON_START;
        $arr[7] = "\"" . $$f{"id"} . "\"";
        $arr[8] = "\"identified by " . $$f{"tool"} . "; putative\"";
        $arr[10] = $TRANSL_TABLE;
        my $aa = &getSeq($fa, $seq, $$f{"sta"}, $$f{"end"}, $$f{"strand"}, "a");
        $arr[11] = '"' . $aa . '"';
        $gff[8] = join(";", "note=" . $arr[7], "note=" . $arr[8],
                       "transl_table=" . $arr[10], "codon_start=" . $arr[4],
                       "translation=" . $arr[11]);
      }
      elsif($$f{"name"} eq "RBS"){
        $arr[7] = "\"identified by " . $$f{"tool"} . "; putative\"";
        $gff[8] = "note=" . $arr[7];
      }
      elsif($$f{"name"} eq "rRNA"){
        $gff[8] = "";
        if($$f{"tool"} ne ""){
          $arr[7] = "\"identified by " . $$f{"tool"} . "; putative\"";
          $gff[8] .= "note=" . $arr[7];
        }
        $gff[8] .= "product=" . $arr[9];
        if(exists $$f{"blast"}){
          $arr[3] = $$f{"blast"};
          $gff[8] .= "blast=" . $arr[7];
        }
        $arr[9] = "\"" . $$f{"type"} . "\"";
      }
      elsif($$f{"name"} eq "tRNA"){
        my $codonpos;
        if($$f{"strand"} eq "+"){
          $codonpos .= sprintf("%d..%d", $$f{"csta"}, $$f{"cend"});
        }
        else{
          $codonpos .= sprintf("complement(%d..%d)", $$f{"csta"}, $$f{"cend"});
        }

        $arr[2] = sprintf("\"(pos:%s,aa:%s,seq:%s)\"",
                          $codonpos, $$f{"type"}, $$f{"codon"});
        $arr[7] = "\"identified by " . $$f{"tool"} . "; putative\"";
        $arr[9] = "\"tRNA-" . $$f{"type"} . "\"";
        $gff[8] = join(";", "anticodon=" . $arr[2], "product=" . $arr[9],
                       "note=" . $arr[7]);
      }

      my $nt = &getSeq($fa, $seq, $$f{"sta"}, $$f{"end"}, $$f{"strand"}, "n");
      $arr[12] = $nt;

      print $foutCsv join(",", @arr), "\n";
      print $foutGff join("\t", @gff), "\n";
    }
  }
  close $foutCsv;
  close $foutGff;
}

