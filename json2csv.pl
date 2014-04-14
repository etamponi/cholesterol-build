#!/usr/bin/env perl

use JSONSAAP;
use strict;

my $jsonFile = $ARGV[0];
my $dataset = $ARGV[1];

my ($result, $Binding, $SProtFT, $Interface, $Relaccess, $Impact, $HBonds,
    $SPhobic, $CPhilic, $BCharge, $SSGeom, $Voids, $MLargest, $NLargest,
    $Clash, $Glycine, $Proline, $CisPro);

open (IN, $jsonFile) || die "Cannot open $jsonFile!\n";

my $content = "";

while(<IN>) {
    $content .= $_;
}

my $jsonText = JSONSAAP::Decode($content);

my ($type, $error) = JSONSAAP::Check($jsonText);

if($error ne "")
{
    print "$error\n";
    exit 1;
}

my ($uniprotac, $res, $nat, $mut) = JSONSAAP::IdentifyUniprotSaap($jsonText);
my @jsonSaaps = JSONSAAP::GetSaapArray($jsonText);

foreach my $jsonSaaps (@jsonSaaps) {
    my ($file, $pdbcode, $chain, $resnum, $mutation) 
        = JSONSAAP::IdentifyPDBSaap($jsonSaaps);
    my ($structuretype,$resolution,$rfactor) 
        = JSONSAAP::GetPDBExperiment($jsonSaaps);
    
    my @analyses = JSONSAAP::ListAnalyses($jsonSaaps);
    foreach my $analysis (@analyses) {
        my $pResults = JSONSAAP::GetAnalysis($jsonSaaps, $analysis);
        
        if ($analysis eq "Binding") {
            $Binding = Binding($pResults);
        } elsif ($analysis eq "SProtFT") {
            $SProtFT = SProtFT($pResults);
        } elsif ($analysis eq "Interface") {
            ($Interface,$Relaccess) = Interface($pResults);
        } elsif ($analysis eq "Impact") {
            $Impact = Impact($pResults);
        } elsif ($analysis eq "HBonds") {
            $HBonds = HBonds($pResults);
        } elsif ($analysis eq "SurfacePhobic") {
            $SPhobic = SPhobic($pResults);
        } elsif ($analysis eq "CorePhilic") {
            $CPhilic = CPhilic($pResults);
        } elsif ($analysis eq "BuriedCharge") {
            $BCharge = BCharge($pResults);
        } elsif ($analysis eq "SSGeom") {
            $SSGeom = SSGeom($pResults);
        } elsif ($analysis eq "Voids") {
            ($Voids,$NLargest,$MLargest) = Voids($pResults);
        } elsif ($analysis eq "Clash") {
            $Clash = Clash($pResults);
        } elsif ($analysis eq "Glycine") {
            $Glycine = Glycine($pResults);
        } elsif ($analysis eq "Proline") {
            $Proline = Proline($pResults);
        } elsif ($analysis eq "CisPro") {
            $CisPro = CisPro($pResults);
        }    
    }    
    
    print "$Binding,$SProtFT,$Interface,$Relaccess,$Impact,$HBonds,$SPhobic,";
    print "$CPhilic,$BCharge,$SSGeom,$Voids,$MLargest,$NLargest,$Clash,";
    print "$Glycine,$Proline,$CisPro,$dataset\n";
}

#-------------------------------------------------------------------------
sub Binding
{
    my($pResults) = @_;
    my $result;
    
    if ($$pResults{'Binding-BOOL'} ne "")
    {
        if    ($$pResults{'Binding-BOOL'} eq 'OK' ) {$result='0';}
        elsif ($$pResults{'Binding-BOOL'} eq 'BAD') {$result='1';} 
    }
    else 
    {
        $result = '?';
    }   
    return($result);
}        
#-------------------------------------------------------------------------
sub SProtFT
{
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'SProtFT-FEATURES'} ne "")
    {
        $result = $$pResults{'SProtFT-FEATURES'} ;
        $result = substr($result, 0, 1).",".substr($result, 1, 1).",".substr($result, 2, 1).",".substr($result, 3, 1).",".substr($result, 4, 1).",".substr($result, 5, 1).",".substr($result, 6, 1).",".substr($result, 7, 1).",".substr($result, 8, 1).",".substr($result, 9, 1).",".substr($result, 10, 1).",".substr($result, 11, 1).",".substr($result, 12, 1);
    }
    else 
    {
        $result = '?,?,?,?,?,?,?,?,?,?,?,?,?';
    }       
    return($result); 
}
#-------------------------------------------------------------------------
sub Interface
{
    my($pResults) = @_;
    my $result1 = ""; 
    my $result2 = "";
    
    if (($$pResults{'Interface-RELACCESS'} ne "") && ($$pResults{'Interface-RELACCESS-MOL'} ne ""))
    {
        $result1 = $$pResults{'Interface-RELACCESS'} - $$pResults{'Interface-RELACCESS-MOL'};
        # Relative accessibility 
        $result2= $$pResults{'Interface-RELACCESS'};
    }
    else 
    {
        $result1 = '?'; 
        $result2 = '?';
    } 
    return($result1, $result2);
}  
#-------------------------------------------------------------------------
sub Impact
{ 
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'Impact-CONSSCORE'} ne "")
    {
        $result = $$pResults{'Impact-CONSSCORE'} ;
    }
    else 
    {
        $result = '?';
    }
    return($result); 
}
#-------------------------------------------------------------------------
#When HBonds-ENERGY is a number I can't remember if positive or negative is good. Consequently you need to look at the range of numbers and make sure that BAD and OK are assigned appropriate numbers. For example, if a BAD HBond has an energy of (for example) +100 and a good HBond has an energy of -100 then when HBonds-ENERGY is NULL you should use -100 for OK and +100 for BAD.

sub HBonds
{
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'HBonds-ENERGY'} ne "")
    {
        if ($$pResults{'HBonds-ENERGY'} eq "NULL")
        {
            if ($$pResults{'HBonds-BOOL'} eq 'OK')
            {
                $result = '20';
            }
            if ($$pResults{'HBonds-BOOL'} eq 'BAD')
            {
                $result = '-20';
            }
        }
        else 
        {
            $result = $$pResults{'HBonds-ENERGY'};
        }
    }
    else
    {
        $result = '?';
    }        
    
    return($result);
}
#-------------------------------------------------------------------------
sub SPhobic
{
    my($pResults) = @_;           
    my $result = "";
    
    if (($$pResults{'SurfacePhobic-NATIVE-HPHOB'} ne "") && ($$pResults{'SurfacePhobic-MUTANT-HPHOB'} ne ""))
    {
        if (($$pResults{'SurfacePhobic-MUTANT-HPHOB'}) > ($$pResults{'SurfacePhobic-NATIVE-HPHOB'}))
        {
            $result = $$pResults{'SurfacePhobic-NATIVE-HPHOB'} - $$pResults{'SurfacePhobic-MUTANT-HPHOB'}; 
        }
        else
        {
            $result = '0';
        }
    }
    else 
    {
        $result = '?'; 
    }
    return($result);           
}           
#-------------------------------------------------------------------------
sub CPhilic
{
    my($pResults) = @_;
    my $result = "";
    
    if (($$pResults{'CorePhilic-NATIVE-HPHOB'} ne "") && ($$pResults{'CorePhilic-MUTANT-HPHOB'} ne ""))
    {
        if ($$pResults{'CorePhilic-MUTANT-HPHOB'} < $$pResults{'CorePhilic-NATIVE-HPHOB'})
        {
            $result = $$pResults{'CorePhilic-NATIVE-HPHOB'} - $$pResults{'CorePhilic-MUTANT-HPHOB'}; 
        }
        else
        {
            $result = '0'; 
        }    
    }
    else 
    {
        $result = '?';
    }
    return($result);         
}
#-------------------------------------------------------------------------
sub BCharge
{
    my($pResults) = @_;
    my $result="";
    
    if (($$pResults{'BuriedCharge-NATIVE-CHARGE'} ne "") && ($$pResults{'BuriedCharge-MUTANT-CHARGE'} ne ""))
    {
        $result = $$pResults{'BuriedCharge-NATIVE-CHARGE'} - $$pResults{'BuriedCharge-MUTANT-CHARGE'}; 
    }
    else 
    {
        $result = '?';
    }
    return($result);        
}
#-------------------------------------------------------------------------
sub SSGeom
{
    my($pResults) = @_; ;
    my $result="";
    
    if ($$pResults{'SSGeom-BOOL'} ne "")
    {
        if    ($$pResults{'SSGeom-BOOL'} eq 'OK')  {$result='0';}
        elsif ($$pResults{'SSGeom-BOOL'} eq 'BAD') {$result='1';}    
    }
    else
    {  
        $result = '?';
    }
    
    return($result);
}
#-------------------------------------------------------------------------
sub Voids
{
    my($pResults) = @_;
    my $result1 = ""; 
    my $result2 = ""; my $result2Ref; my $Ntotal=0;
    my $result3 = ""; my $result3Ref; my $Mtotal=0;
    
    if (($$pResults{'Voids-NATIVE-LARGEST'} ne "") && ($$pResults{'Voids-MUTANT-LARGEST'} ne ""))
    {
        $result1 = $$pResults{'Voids-NATIVE-LARGEST'} - $$pResults{'Voids-MUTANT-LARGEST'}; 
		
    }
    else
    {
        $result1 = '0';  # shall i change this to zero? or set it to ?
    } 
    
    if ($$pResults{'Voids-NATIVE'} ne "")
    {
        $result2Ref = $$pResults{'Voids-NATIVE'}; # this is an array 
        $Ntotal = scalar (@$result2Ref);
        
        $result2 = readArrayRef($result2Ref,$Ntotal);

    }
    else 
    {
        $result2 = '0,0,0,0,0,0,0,0,0,0';
    } 
    
    if ($$pResults{'Voids-MUTANT'} ne "")
    {
        $result3Ref = $$pResults{'Voids-MUTANT'}; # this is an array
        $Mtotal = scalar (@$result3Ref);
        
        $result3 = readArrayRef($result3Ref,$Mtotal);
		
	    }
    else
    {
        $result3 =  '0,0,0,0,0,0,0,0,0,0';
    }
    
    return($result1, $result2, $result3);        
}
#-------------------------------------------------------------------------
sub readArrayRef
{
my($result2Ref, $Ntotal) = @_;
my $result;
if    ($Ntotal == 10){$result ="$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],$result2Ref->[7],$result2Ref->[8],$result2Ref->[9]";}
        elsif ($Ntotal == 9) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],$result2Ref->[7],$result2Ref->[8],0";}
        elsif ($Ntotal == 8) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],$result2Ref->[7],0,0"; }
		elsif ($Ntotal == 7) {$result ="$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],$result2Ref->[6],0,0,0"; }
        elsif ($Ntotal == 6) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],$result2Ref->[5],0,0,0,0";}
        elsif ($Ntotal == 5) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],$result2Ref->[4],0,0,0,0,0";}
		elsif ($Ntotal == 4) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],$result2Ref->[3],0,0,0,0,0,0";}
        elsif ($Ntotal == 3) {$result = "$result2Ref->[0],$result2Ref->[1],$result2Ref->[2],0,0,0,0,0,0,0"; }
        elsif ($Ntotal == 2) {$result = "$result2Ref->[0],$result2Ref->[1],0,0,0,0,0,0,0,0";  }
		else                 {$result = "$result2Ref->[0],0,0,0,0,0,0,0,0,0"; }
      return ($result);
}		
#-------------------------------------------------------------------------
sub Clash
{
    my($pResults) = @_;
    my $result="";
    
    if ($$pResults{'Clash-ENERGY'} ne "") 
    {
        $result = $$pResults{'Clash-ENERGY'}
    }
    else
    {
        $result = '?';
    }
    return($result);
}
#-------------------------------------------------------------------------
sub Glycine
{
    my($pResults) = @_;
    my $result="";
    
    if (($$pResults{'Glycine-NATIVE'} eq "GLY") && ($$pResults{'Glycine-MUTANT'} ne "GLY"))
    {
        if(($$pResults{'Glycine-NATIVE-ENERGY'} ne "") && ($$pResults{'Glycine-MUTANT-ENERGY'} ne ""))
        {
            $result = $$pResults{'Glycine-NATIVE-ENERGY'} - $$pResults{'Glycine-MUTANT-ENERGY'}; 
        }
        else 
        {
            $result = '?';
        }
    }
    else
    {
        #if not mutant to Gly we set the diffrent in energy to -10 or -100
        $result = '-100';
    }
    return($result);        
}
#-------------------------------------------------------------------------
sub Proline
{
    my($pResults) = @_;
    my $result = "";
    
    if (($$pResults{'Proline-MUTANT'} eq "PRO") && ($$pResults{'Proline-NATIVE'} ne "PRO"))
    {
        if (($$pResults{'Proline-MUTANT-ENERGY'} ne "") && ($$pResults{'Proline-NATIVE-ENERGY'} ne ""))
        {
            $result = $$pResults{'Proline-MUTANT-ENERGY'} - $$pResults{'Proline-NATIVE-ENERGY'}; 
        }
        else 
        {
            $result = '?';
        }
    }
    else
    {
        #if not mutant from Pro we set the diffrent in energy to -10 or -100
        $result = '-100';
    }
    return($result);        
}
#-------------------------------------------------------------------------
sub CisPro
{
    my($pResults) = @_;
    my $result = "";
    
    if ($$pResults{'CisPro-BOOL'} ne "")
    {
        if ($$pResults{'CisPro-BOOL'} eq 'OK')  {$result='0';}
        elsif ($$pResults{'CisPro-BOOL'} eq 'BAD') {$result='1';}    
    }
    else {$result = '?';}
    return($result);
}
#-------------------------------------------------------------------------

