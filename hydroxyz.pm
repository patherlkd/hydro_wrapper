#!/usr/bin/perl

package hydro::hydroxyz;

use v5.10;
use strict;
use Scalar::Util qw{ openhandle };
use Math::Trig ':pi';
use lib "./";
use Exporter qw(import);
use _Perl::List qw{ contains };
use Scalar::Util qw(looks_like_number); 
use File::Basename;
use Cwd;


use constant DaToGram => 1.6605e-24; # factor (when multiplied) converts Da to g
use constant kB => 1.38064e-23; # Boltzmann's constant  (m^2 kg s^-2 k^-1)
use constant PoiseToPas => 0.1; # factor (when multiplied) converts poise to pascal x seconds = kg m^-1 s^-1
use constant DcmToDm => 1e-4; # factor (when multiplied) converts diffusion coefficient in cm^2/s to m^2/s
use constant DegToKel => 273.15; # To convert degrees to kelvin (add)

sub new
{
    my ($class,$args) = @_;
    
    my $self = bless {
	_exe => $args->{exe},
	_xyzfile => $args->{xyzfile},
	_res_dir => $args->{res_dir},
	_debug => $args->{debug},
	_atoms => 0,
	_traj => 0,
	_skip => $args->{skip}, 
	_ptypes => $args->{ptypes},
	_diameters => $args->{diameters},
	_length_scale => $args->{length_scale},
	_T => $args->{temperature},
	_eta => $args->{viscosity},
	_MWs => $args->{molweights},
	_rho => $args->{density},
	_apbs => $args->{aminosperbeads},
	_icase => $args->{icase},
	_hvals => $args->{hvalues} || 0, 
	_hmax => $args->{hmax} || 0,
	_rvals => $args->{rvalues} || 0, 
	_rmax => $args->{rmax} || 0,
	_trials => $args->{ntrials} || 0,
	_idif => $args->{idif} || 0,
    } ,$class;


    return $self;
}

sub openfile
{
    
    my $self=shift;
    open ($self->{_xyz},'<',$self->{_xyzfile}) || die "Could not open $self->{_xyzfile}";

    say "\nOpened $self->{_xyzfile}" if $self->{_debug};
    set_atoms($self);

}


sub dotrajectory
{ # run hydro++ for one trajectory
    my $self=shift;
    
    my $atomcount = 1;
    my @cols = (1);

    my $mw=0.0; # mol weight in Da
    my $m=0.0; # mass in g
    my $vol=0.0; # vol in cm^3
    my $vbar=0.0; # partial specific vol in cm^3/g

    my $simname = fileparse($self->{_xyzfile},'\.[^\.]*');

    my $namelength = length($simname);
    
    if($namelength>20){
	$simname = substr($simname,$namelength - 20);
    }

    if(!$self->{_name}){
    $self->{_name} = $simname;
    }

    
    
    my $tempcoords = $simname."tempcd";
    my $tempinput = $simname."tempin";

    open(COORDS,'>',$tempcoords) || die;
    open(INPUT,'>',$tempinput) || die;

    say COORDS "$self->{_length_scale},";
    say COORDS "$self->{_atoms},";

    while($cols[0]!=$self->{_atoms})
    {

	if(eof($self->{_xyz})){ return 0; };
	chomp(my $line = readline($self->{_xyz}));
	
	@cols = split(/\s+/,$line);
       
	if(!looks_like_number($cols[0]))
	{ $cols[0] = 0;  
	}
	elsif(contains($self->{_ptypes},$cols[0]))
	{
	
	    # cartesian coords
	    my $ptype=$cols[0];
	    my $x = sprintf("%.4f",$cols[1]);
	    my $y = sprintf("%.4f",$cols[2]);
	    my $z = sprintf("%.4f",$cols[3]);
	    
	    my $dia=${$self->{_diameters}}{$ptype};
	    my $molweight=${$self->{_MWs}}{$ptype};
	    my $apb=${$self->{_apbs}}{$ptype};

	    $vol += 4.0/3.0*pi*(${dia}*0.5*$self->{_length_scale})**3;
	    
	    $mw += $molweight*$apb; 
	    
	    say COORDS "  $x  $y  $z  $dia";
	}
       
    }



    say " volume of mol = $vol (cm^3) " if $self->{_debug};
    say " mol weight of mol = $mw (Da) " if $self->{_debug};
    $m = $mw * DaToGram;
    say " mass of mol = $m (g) " if $self->{_debug};
    $vbar = $vol/$m;
    say " partial specific volume of mol = $vbar (cm^3/g) " if $self->{_debug};

    my $output = $simname."$self->{_traj}";

    if(length($simname)>20){ say "simname > 20 chars long";}
    if(length($output)>30){ say "results filename > 30 chars long";}
    if(length($tempcoords)>30){ say "coords filename > 30 chars long";}

    

    say INPUT $simname;
    say INPUT $output;
    say INPUT $tempcoords;
    say INPUT $self->{_icase};
    say INPUT $self->{_T};
    say INPUT $self->{_eta};
    say INPUT $mw;
    say INPUT sprintf("%.3f",$vbar);
    say INPUT $self->{_rho};
    say INPUT $self->{_hvals}.",";
    say INPUT $self->{_hmax}.",";
    say INPUT $self->{_rvals}.",";
    say INPUT $self->{_rmax};
    say INPUT $self->{_trials}.",";
    print INPUT $self->{_idif};
    print INPUT "                   lol\n";
    say INPUT "\*           End of file";
    
   

    close(COORDS);
    close(INPUT);
    
    # Run hydro++ 


   `$self->{_exe} <<< $tempinput 2>&1`;

    chomp($self->{_results}=`ls ${simname}*-res.txt`);

    $self->{_traj}++;    
    return 1; # succeeded
}


sub stokes
{ # computes the stokes radius
    
    my $self = shift;
    say "Computing stokes radius" if $self->{_debug};

    my $D=$self->{_trD};
    my $Rs = kB*($self->{_T}+DegToKel);
    $Rs = $Rs/(6.0*pi*($self->{_eta}*PoiseToPas)*($D*DcmToDm));
    $Rs=$Rs*100;

    $Rs/= $self->{_length_scale}; # convert to normal length scale from user


    return $Rs;
}

sub gyration
{ # computes radius of gyration
    my $self = shift;
    return quantity($self,"Radius of gyration")/$self->{_length_scale};
}

sub trans_D
{ # get the computed translational diffusion constant cm^2/s
    
    my $self = shift;
    return quantity($self,"Translational diffusion coefficient:");  
}

sub quantity
{
    my ($self,$pat) =@_;

    my $res=$self->{_results};
    
    if(!$res){ return 0; }

    open (IN,'<',$res) || die;
    my @lines = <IN>;
    close(IN);

    
    say "Computing $pat" if $self->{_debug};

    my $D = 0;

    for my $line (@lines)
    {
	if($line =~ /$pat/)
	{
	    my @parts = split(/\s+/,$line);
	    $D = getvariable(\@parts);
	}
    }

    $self->{_trD}=$D;
    return $D;
}

sub saveresult
{
    my $self = shift;

    my ($resname,$res) = @_;
    
    my $file = $self->{_res_dir}.$self->{_name}."_temp$self->{_T}_visc$self->{_eta}"."_".$resname.".txt";

    open(OUT,'>',$file) || die;

    say OUT $res;

    close(OUT);
    
    say "Saved $resname = $res to $file" if $self->{_debug};

}

sub getvariable
{
    my ($ref) = @_;
    
    my @array = @{$ref};
    my $var =0;

    
    for my $x (@array)
    {
	if(looks_like_number($x)){$var =$x; }
    }

    return $var;

}

sub fastforward
{ # skip unequilibrated configurations

    my $self=shift;

    my ($add) = @_;

    if($add){ $self->{_skip}=$self->{_traj}+$add; }

    say "Fast forwarding" if $self->{_debug};

    while($self->{_traj} < $self->{_skip})
    {
	if(eof($self->{_xyz})){ return 0; }
	chomp(my $line=readline($self->{_xyz}));
	
	my @cols = split(/\s+/,$line);
	if($cols[0] eq "$self->{_atoms}"){ $self->{_traj}++;  };
    }

    say "Traj is now: $self->{_traj}" if $self->{_debug};

}


sub set_atoms
{
    my $self=shift;

    chomp($self->{_atoms} = readline($self->{_xyz}));
    say "Atoms = $self->{_atoms}" if $self->{_debug};
}

sub clean
{
    my $self=shift;
    my @files = glob("*");

    foreach my $file (@files){
	
	if( -d $file ){ ; }
	elsif( $file !~ /hydro_main/ )
	{
	    unlink $file or warn "Could not unlink $file";
	}
    }
}

sub DESTROY
{
    my $self=shift;
    
    close(openhandle($self->{_xyz}));    

    $self->{handle}->close() if $self->{handle};
}


1;
