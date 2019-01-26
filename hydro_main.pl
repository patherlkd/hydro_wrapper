#!/usr/bin/perl -w

use v5.10;
use Math::Trig;
use lib "./";
use hydro::hydroxyz;
use strict;


my $debug = $ARGV[0]; # 1 to show calculations on command line  0 to show nothing
my $dir = $ARGV[1]; # dir containing xyz files
my $resdir = $ARGV[2]; # dir to output results too
my @xyzfiles= glob "${dir}/*xyz"; # sim trajectories

foreach my $file (@xyzfiles){

    my @ptypes = (1,2,4); # particle types from lammps
    my %diameters = ( 1 => 0.38*1.5, # e.g. 1.5 times size of amino acid
		      2 => 0.38*1.5,
		      4 => 0.38*1.5,
	); # bead diameters
    my %molweights = ( 1 => 110*1.5,
		       2 => 110*1.5,
		       4 => 110*1.5,
	); # molecular weights in Da (ave amino acid
    my %aminosperbeads = ( 1 => 1.5,
			   2 => 1.5,
			   4 => 1.5,
	); # amino acids in a bead

    my $sim = hydro::hydroxyz->new( # create the class to analyse sim trajectories
	{exe => "./hydro++", # hydro++ executable name to run in terminal
	    xyzfile => $file,
	 res_dir => $resdir,
	 diameters => \%diameters, # diameter * length_scale = real diameter
	 molweights => \%molweights,
	 aminosperbeads => \%aminosperbeads,
	 skip => 400,
	 debug => $debug, # 1 is to debug
	 length_scale => "1.E-07", # in cm
	 temperature => 20, # in degrees centigrade
	 viscosity => 0.01, # in Poise
	 density => 1.0, # density of solution g/cm^3
	 icase => 12, # ICASE in hydro++ program
	 idif => 1, # for full diffusion tensors
	 ptypes => \@ptypes,
	}
	);

    $sim->openfile();

    my $N=1;
    my $Rs=0;
    my $Rs2=0;

   my $D=0;
    my $Rg=0;

    while($sim->dotrajectory)
    {

	$D+=$sim->trans_D; 
	$Rs+=$sim->stokes;
	$Rs2+=($Rs*$Rs);
	$Rg+=$sim->gyration;
	$N++;
	$sim->fastforward(10); #skip N trajs	

    }  

    $D /= $N;
    $Rs /= $N;
    $Rs2 /= $N;
    $Rg /= $N;
    
    my $variance = sqrt($Rs2 - $Rs*$Rs);

    $sim->saveresult("transD",$D);
    $sim->saveresult("stokes",$Rs);
    $sim->saveresult("Rg",$Rg);
    $sim->saveresult("variance",$variance);


    $sim->clean;
}

