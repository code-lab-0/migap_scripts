#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $fFile  = ""; # fasta
my @cFiles = (); # Blast to cog result
my @rFiles = (); # Blast to refseq result
my $out    = ""; # output

GetOptions("fa=s"         => \$fFile,
           "cog=s{,}"     => \@cFiles,
           "ref=s{,}"     => \@rFiles,
           "out=s"        => \$out);

my @inFiles = (@cFiles, @rFiles);
my $mapped = &seekBlastPfiles(\@inFiles, 60, 60);

&generateFa($fFile, $out, $mapped);

#---------------------------------------------------
# functions

# fasta with not mapped to cog or refseq
sub generateFa {
  my ($in, $out, $mapped) = @_;
  my $flag = 0;
  open my $fout, ">", $out or die;
  open my $fin, "<", $in or die;
  while(<$fin>){
    chomp;
    my $line = $_;
    if(/^>(\S+)/){
      $flag = 0;
      $flag = 1 if exists $$mapped{$1};
    }
    print $fout $line, "\n" if $flag == 0;
  }
  close $fin;
  close $fout;
}

# extract query id from blast result
sub seekBlastPfiles {
  my ($files, $idenTh, $ovTh) = @_;
  my $mapped;
  foreach my $file(@$files){
    $mapped = &seekBlastP($mapped, $file, $idenTh, $ovTh);
  }
  return $mapped;
}

# parse blastp result
sub seekBlastP {
  my ($mapped, $file, $idenTh, $ovTh) = @_;

  my %seqinfo = ();
  my %info = ();

  open my $fin, '<', $file or die;
  while(<$fin>){
    chomp;

    if(/^Query= (\S+)/){
      $seqinfo{"que"} = $1;
    }
    elsif(/^>\s*(\S.+)$/){
      $seqinfo{"inTar"} = 1;
    }
    elsif(/^Length=(\d+)/ && exists $seqinfo{"inTar"}){
      $seqinfo{"len"} = $1;
      delete $seqinfo{"inTar"}
    }
    elsif(/^ Identities = (\d+)\/(\d+).+ Gaps = (\d+)\/\d+/){
      $info{"iden"} = $1;
      $info{"alen"} = $2;
    }
    elsif(/^(Query\s+(\d+)\s+)(\S+)\s+(\d+)/){
      if(! exists $info{"qsta"}){
        $info{"qsta"} = $2;
        $info{"qseq"} = "";
        $info{"sseq"} = "";
      }
      $info{"qseq"} .= $3;
    }
    elsif(/^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/){
      if(! exists $info{"ssta"}){
        $info{"ssta"} = $1;
      }
      $info{"sseq"} .= $2;
    }
    elsif(/^ Score =\s*(\d+\.?\d*) bits.+Expect = (\S+),/){
      if(exists $info{"sco"}){
        my $ov = &calcOv($info{"alen"}, length($info{"qseq"}),
                         length($info{"sseq"}));
#print STDERR join("\t", $seqinfo{"que"}, $info{"iden"}, $ov), "\n";
        if($info{"iden"} >= $idenTh && $ov >= $ovTh){
          $$mapped{$seqinfo{"que"}} = 1;
        }
      }
      %info = ();
      $info{"sco"} = sprintf("%.1f", $1);
    }
    elsif(/^Lambda/){
      if(exists $info{"sco"}){
        my $ov = &calcOv($info{"alen"}, length($info{"qseq"}),
                         length($info{"sseq"}));
#print STDERR join("\t", $seqinfo{"que"}, $info{"iden"}, $ov), "\n";
        if($info{"iden"} >= $idenTh && $ov >= $ovTh){
          $$mapped{$seqinfo{"que"}} = 1;
        }
      }
      %info = ();
    }

  }
  close $fin;

  return $mapped;
}

sub calcOv {
  my ($alen, $qlen, $slen) = @_;
  my $longer = $qlen <= $slen ? $slen : $qlen;
  my $ov = int($alen * 100 / $longer);
  return $ov; 
}

