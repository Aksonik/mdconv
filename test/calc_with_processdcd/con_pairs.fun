#!/bin/perl

my $cnt=();
my $ncnt=();

sub start {
$f=0;
$pairs=0; 	### total number of pairs

$pn_11=0;
$pn_22=0;
$pn_33=0;
$pn_12=0;
$pn_13=0;
$pn_23=0;

open(my $file,'<',"con_p1.list");
while(<$file>)
{chomp;
 push @list_p1, $_;
}
print join ",",@list_p1,"\n";
close $file;

open(my $file,'<',"con_p2.list");
while(<$file>)
{chomp;
 push @list_p2, $_;
}
print join ",",@list_p2,"\n";
close $file;

open(my $file,'<',"con_p3.list");
while(<$file>)
{chomp;
push @list_p3, $_;
}
print join ",",@list_p3,"\n";
close $file;

}

sub end {

open(my $file,'>',"pairs.res");
for($i=1;$i<=$pairs;$i++)
 {print $file "$pair_t[$i] $pair_p[$i] $pair_pp[$i]\n";
 }
close $file;

open(my $file,'>',"pairs11.res");
for($i=1;$i<=$pn_11;$i++)
 {print $file "$pair_t_11[$i] $pair_p_11[$i] $pair_pp_11[$i]\n";
 }
close $file;

open(my $file,'>',"pairs22.res");
for($i=1;$i<=$pn_22;$i++)
 {print $file "$pair_t_22[$i] $pair_p_22[$i] $pair_pp_22[$i]\n";
 }
close $file;

open(my $file,'>',"pairs33.res");
for($i=1;$i<=$pn_33;$i++)
 {print $file "$pair_t_33[$i] $pair_p_33[$i] $pair_pp_33[$i]\n";
 }
close $file;

open(my $file,'>',"pairs12.res");
for($i=1;$i<=$pn_12;$i++)
 {print $file "$pair_t_12[$i] $pair_p_12[$i] $pair_pp_12[$i]\n";
 }
close $file;

open(my $file,'>',"pairs13.res");
for($i=1;$i<=$pn_13;$i++)
 {print $file "$pair_t_13[$i] $pair_p_13[$i] $pair_pp_13[$i]\n";
 }
close $file;

open(my $file,'>',"pairs23.res");
for($i=1;$i<=$pn_23;$i++)
 {print $file "$pair_t_23[$i] $pair_p_23[$i] $pair_pp_23[$i]\n";
 }
close $file;

open(my $file,'>',"con.res");
for($i=1;$i<=$f;$i++)
 {print $file "$con_t[$i] $con[$i]\n";
 }
close $file;

open(my $file,'>',"con11.res");
for($i=1;$i<=$f;$i++)
 {print $file "$con_t[$i] $con11[$i]\n";
 }
close $file;

open(my $file,'>',"con22.res");
for($i=1;$i<=$f;$i++)
 {print $file "$con_t[$i] $con22[$i]\n";
 }
close $file;

open(my $file,'>',"con33.res");
for($i=1;$i<=$f;$i++)
 {print $file "$con_t[$i] $con33[$i]\n";
 }
close $file;

open(my $file,'>',"con12.res");
for($i=1;$i<=$f;$i++)
 {print $file "$con_t[$i] $con12[$i]\n";
 }
close $file;

open(my $file,'>',"con13.res");
for($i=1;$i<=$f;$i++)
 {print $file "$con_t[$i] $con13[$i]\n";
 }
close $file;

open(my $file,'>',"con23.res");
for($i=1;$i<=$f;$i++)
 {print $file "$con_t[$i] $con23[$i]\n";
 }
close $file;

}

sub analyze {

$f=$f+1;		### frame
$dt=1.00;
$tt=$f*$dt;

print "time: $tt\n";

my $mol=shift;
my $boxa=shift;
my $boxb=shift;
my $boxc=shift;

#print "box: $boxa $boxb $boxc\n";

my $totchains=$#{$mol->{chain}};		### total number of chains (from zero)
my $totat=$#{$mol->{chain}->[0]->{atom}};	### total number of atoms (from zero)

print "totat: $totat\n";

$cc=0;		### tot. number of contacts
$cc11=0;		### only between 1vii
$cc22=0;		### only between 1gbq
$cc33=0;		### only between 1ubq

$cc12=0;
$cc13=0;
$cc23=0;

for (my $p=1; $p<=$totat; $p++)				### protein A
{for (my $pp=$p+1; $pp<=$totat+1; $pp++)			### protein B
 {

  my $ax=$mol->{chain}->[0]->{atom}->[$p-1]->{xcoor};
  my $ay=$mol->{chain}->[0]->{atom}->[$p-1]->{ycoor};
  my $az=$mol->{chain}->[0]->{atom}->[$p-1]->{zcoor};

  my $aax=$mol->{chain}->[0]->{atom}->[$pp-1]->{xcoor};
  my $aay=$mol->{chain}->[0]->{atom}->[$pp-1]->{ycoor};
  my $aaz=$mol->{chain}->[0]->{atom}->[$pp-1]->{zcoor};

  my $dx=$ax-$aax;
  my $dy=$ay-$aay;
  my $dz=$az-$aaz;

  $dx-=$boxa*&nint($dx/$boxa);
  $dy-=$boxb*&nint($dy/$boxb);
  $dz-=$boxc*&nint($dz/$boxc);

  my $d=$dx*$dx+$dy*$dy+$dz*$dz;

# rg:
# 1vii: 1.063
# 1gb1: 1.188
# 1ubq: 1.294

  if( $p ~~ @list_p1 )
  {$rgA=9.06365}
  if( $pp ~~ @list_p1 )
  {$rgB=9.06365}
  if( $p ~~ @list_p2 )
  {$rgA=10.3346}
  if( $pp ~~ @list_p2 )
  {$rgB=10.3346}
  if( $p ~~ @list_p3 )
  {$rgA=11.67230}
  if( $pp ~~ @list_p3 )
  {$rgB=11.67230}

  $A=1.00;
  $D=7;

  $rc=($A*($rgA+$rgB)+$D)**2;

  print "f: $f pA: $p pB: $pp d: $d c: $rc rgyr: $rgA $rgB\n";

  if($d<=$rc)					### cutoff: sqrt(49)=7A
  {$cc=$cc+1;

   if(( $p ~~ @list_p1 )&&( $pp ~~ @list_p1 ))
   {$cc11=$cc11+1;
   }

   if(( $p ~~ @list_p2 )&&( $pp ~~ @list_p2 ))
   {$cc22=$cc22+1;
   }

   if(( $p ~~ @list_p3 )&&( $pp ~~ @list_p3 ))
   {$cc33=$cc33+1;
   }

   if((( $p ~~ @list_p1 )&&( $pp ~~ @list_p2 ))||(( $p ~~ @list_p2 )&&( $pp ~~ @list_p1 )))
   {$cc12=$cc12+1;
   }

   if((( $p ~~ @list_p1 )&&( $pp ~~ @list_p3 ))||(( $p ~~ @list_p3 )&&( $pp ~~ @list_p1 )))
   {$cc13=$cc13+1;
   }

   if((( $p ~~ @list_p2 )&&( $pp ~~ @list_p3 ))||(( $p ~~ @list_p3 )&&( $pp ~~ @list_p2 )))
   {$cc23=$cc23+1;
   }


#    print "f: $f pA: $p pB: $pp d: $d c: $rc\n";
#    print "con: $tt $p $pp\n";

    $pairs=$pairs+1;

    $pair_t[$pairs]=$tt;	### save time
    $pair_p[$pairs]=$p;		### save 1st protein from the contact pair
    $pair_pp[$pairs]=$pp;	### save 2nd protein from the contact pair

   if(( $p ~~ @list_p1 )&&( $pp ~~ @list_p1 ))
   {$pn_11=$pn_11+1;
    $pair_p_11[$pn_11]=$p;
    $pair_pp_11[$pn_11]=$pp;
    $pair_t_11[$pn_11]=$tt;
   }

   if(( $p ~~ @list_p2 )&&( $pp ~~ @list_p2 ))
   {$pn_22=$pn_22+1;
    $pair_p_22[$pn_22]=$p;
    $pair_pp_22[$pn_22]=$pp;
    $pair_t_22[$pn_22]=$tt;
   }

   if(( $p ~~ @list_p3 )&&( $pp ~~ @list_p3 ))
   {$pn_33=$pn_33+1;
    $pair_p_33[$pn_33]=$p;
    $pair_pp_33[$pn_33]=$pp;
    $pair_t_33[$pn_33]=$tt;
   }

   if((( $p ~~ @list_p1 )&&( $pp ~~ @list_p2 ))||(( $p ~~ @list_p2 )&&( $pp ~~ @list_p1 )))
   {$pn_12=$pn_12+1;
    $pair_p_12[$pn_12]=$p;
    $pair_pp_12[$pn_12]=$pp;
    $pair_t_12[$pn_12]=$tt;
   }

   if((( $p ~~ @list_p1 )&&( $pp ~~ @list_p3 ))||(( $p ~~ @list_p3 )&&( $pp ~~ @list_p1 )))
   {$pn_13=$pn_13+1;
    $pair_p_13[$pn_13]=$p;
    $pair_pp_13[$pn_13]=$pp;
    $pair_t_13[$pn_13]=$tt;
   }

   if((( $p ~~ @list_p2 )&&( $pp ~~ @list_p3 ))||(( $p ~~ @list_p3 )&&( $pp ~~ @list_p2 )))
   {$pn_23=$pn_23+1;
    $pair_p_23[$pn_23]=$p;
    $pair_pp_23[$pn_23]=$pp;
    $pair_t_23[$pn_23]=$tt;
   }


  }


 }
}

$con_t[$f]=$tt;		### save time
$con[$f]=$cc;		### save total number of contacts

$con11[$f]=$cc11;
$con22[$f]=$cc22;
$con33[$f]=$cc33;

$con12[$f]=$cc12;
$con13[$f]=$cc13;
$con23[$f]=$cc23;

# print "tot con: $con_t[$f] $con[$f]\n";

}

sub nint {
  my $x = $_[0]; 
  my $n = int($x);
  if ( $x > 0 ) {
    if ( $x-$n > 0.5) {
      return $n+1;
    }
    else {
      return $n;
    }
  }
  else {
    if ( $n-$x > 0.5) {
      return $n-1;
    }
    else {
      return $n;
    }
  }
}
1;
