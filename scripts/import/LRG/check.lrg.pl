=pod

SYNOPSIS

    Script to run some healthchecks on a LRG XML record

DESCRIPTION

    Will run a number of healthchecks on a LRG record and check for inconsistencies or obvious mistakes
  
EXAMPLE
  
    List all available checks:
        perl check.lrg.pl -list_checks
        
    Run all available checks on LRG_1:
        perl check.lrg.pl -xml_file LRG_1.xml
  
    Validate the schema of LRG_1 XML file:
        perl check.lrg.pl -xml_file LRG_1.xml -check schema
          
=cut

#!perl -w

use strict;

use Getopt::Long;
use LRG::LRGHealthcheck;

usage() if (!scalar(@ARGV));

my $xml_file;
my @checks;
my $java_executable;
my $jing_jar;
my $rnc_file;
my $list;
my $verbose;
my $help;

# get options from command line
GetOptions(
  'xml_file=s'		=> \$xml_file,
  'check=s'		=> \@checks,
  'java=s'		=> \$java_executable,
  'jing=s'		=> \$jing_jar,
  'rnc=s' 		=> \$rnc_file,
  'list_checks!'        => \$list,
  'verbose!' 		=> \$verbose,
  'help!'               => \$help
);

usage() if (defined($help));

if ($list) {
    print "\nAvailable checks:\n\n";
    foreach my $check (@LRGHealthcheck::CHECKS) {
        print "\t$check\n";
    }
    print "\n";
    exit(0);
}

@checks = @LRGHealthcheck::CHECKS if (!defined(@checks));
$LRGHealthcheck::JAVA = $java_executable if (defined($java_executable));
$LRGHealthcheck::JING_JAR = $jing_jar if (defined($jing_jar));
$LRGHealthcheck::RNC_FILE = $rnc_file if (defined($rnc_file));

my $hc = LRGHealthcheck::new($xml_file);
foreach my $check (@checks) {
    if (!grep(/^$check$/,@LRGHealthcheck::CHECKS)) {
        print "Unknown healthcheck '$check'\n";
        next;
    }
    $hc->$check();
}

foreach my $check (@checks) {
    print "$check\t" . ($hc->{'check'}{$check}{'passed'} ? "PASSED" : "FAILED") . "!\n" if ($verbose || !$hc->{'check'}{$check}{'passed'});
    if (exists($hc->{'check'}{$check}{'message'})) {
        print "\t" . join("\n\t\t",split(/\/\//,$hc->{'check'}{$check}{'message'})) . "\n" if ($verbose || !$hc->{'check'}{$check}{'passed'});
    }
}

sub usage {
    
  print qq{
  Usage: perl check.lrg.pl [OPTION]
  
  Run a series of healthchecks on a LRG XML record
	
  Options:
    
        -xml_file       Path to LRG XML record to be checked (required)
        -check          Name of the check(s) to run (multiple checks can be specified). By default,
                        all available checks are run. To see a list of available checks, use parameter
                        -list_checks
        
    The check to validate the schema launches the jing application. This needs correct paths to the java
    executable, the jing jar file and the RelaxNG Compact schema definition file. To override the default
    settings, use the options below:
    
        -java           Path to java executable. Default is to assume 'java' is in the system's path
        -jing           Path to the jing executable jar file. Default is to look in the current working directory.
        -rnc            Path to RelaxNG Compact schema definition file. Default is to look in the current
                        working directory.
                        
        -verbose        Print some progress information
        -list_checks    Print a list of the available healthchecks and exit
        -help           Print this message
        
  };
  exit(0);
}