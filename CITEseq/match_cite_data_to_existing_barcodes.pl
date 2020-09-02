#Barcode structure: 
#R1 ^[16] = cell barcode
#R2 ^10[15] = antibody barcode
#I1 = umi.

if (@ARGV !=3 ) {die "Usage: grab_antibody_codes_efficient.pl Read1.fastq Read2.fastq cellbarcodes.txt\n";}


# Reference antibody barcodes
%refcodes =();
#open($ifhr, "/home/wizard/Tallulah/HumanLiverDeepDive/100819QG_TotalSeq-C_human_Panel_list.csv") or die $!;
open($ifhr, "/cluster/projects/macparland/TA/DeepDive/100819QG_TotalSeq-C_human_Panel_list.csv") or die $!;
<$ifh1>; # header
while(<$ifhr>) {
  chomp;
  @stuff = split(/\t/);
  $code = $stuff[6];
  $code =~ s/"//g;
#  print $stuff[2]."\t".$code."\n";
  $refcodes{$code}=0;
}
close($ifhr);
#exit();

%refcells=();
open($ifhb, $ARGV[2]) or die $!;
while(<$ifhb>) {
   chomp;
   $code =$_;
   $code =~ s/"//g;
 #  print $stuff[2]."\t".$code."\n";
   $refcells{$code}=0;
}


# write in chunks;
my @lines = ();

# Read in reads & extract barcodes;
my $nocell = 0;#number of reads not matched to a cell;
my $noantibody = 0;#number of reads not matched to an antibody;
my $nreads =0; #number of reads;

%out =();
open($ifh2, $ARGV[1]) or die $!;
open($ifh1, $ARGV[0]) or die $!;
my $cellcode = "";
my $umi = "";
my $antibody = "";
while (<$ifh2>) {
  my $write_it = 0;
  my $r1_line = <$ifh1>;
    chomp;
    if ($_ =~ /^\@/){
      $nreads++;
      @stuff = split(/\s+/);
#umi from readname
      $umi = $stuff[1];
      $umi =~/.+:([ATCG]+)/;
      $umi = $1;
# cell barcode from R1
      my $r1_seq = <$ifh1>;
      $cellcode = substr($r1_seq, 0, 16);
      if (!exists($refcells{$cellcode})) {
        foreach my $cell (keys(%refcells)) {
          my $count = ( $cell ^ $cellcode ) =~ tr/\0//;
          if ($count >= length($cellcode)-1) { #allow upto 1 mismatch
            $cellcode = $cell;
            last;
          }
        }
      }
      if (!exists($refcells{$cellcode})) {
        $nocell++;
        next;
      }

# antibody barcode from R2
      my $seq = <$ifh2>;
      $antibody = substr($seq, 10, 15);
      if (exists($refcodes{$antibody})) {
        $refcodes{$antibody}++;
        $write_it =1;
      } else {
        foreach my $code (keys(%refcodes)) {
        if ($seq =~ $code) {
          $write_it =1;
          $refcodes{$code}++;
          $antibody=$code;
        }
        }
      }
      if (!$write_it) {$noantibody++;}
    }
    if ($write_it) {
      if (@lines > 30000) {
        foreach $l (@lines) {
          print $l;
        }
        @lines = ();
      }
      push(@lines, "$cellcode\t$antibody\t$umi\n");
    }
}


close($ifh1);
close($ifh2);


foreach $l (@lines) {
   print $l;
}

print STDERR ("Reads: $nreads\n Nocell: $nocell\n Noantibody: $noantibody\n");

#foreach my $code (keys(%refcodes)) {
#  print STDERR ($code." ".$refcodes{$code}."\n");
#}

