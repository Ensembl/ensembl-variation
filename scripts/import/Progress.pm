use strict;
use warnings;

# Package for outputting progress information

package Progress;

sub location {
    my ($pack,$file,$line) = caller();
    my $str = localtime() . "\t\tAt " . $pack . "::" . $file . ", line $line\n";
    return $str;
}


sub time_format {
    my $time = shift;
    
    my $minute = 60;
    my $hour = 60*$minute;
    my $day = 24*$hour;
    my $week = 7*$day;
    
    my $weeks = int($time/$week);
    $time -= $weeks*$week;
    
    my $days = int($time/$day);
    $time -= $days*$day;
    
    my $hours = int($time/$hour);
    $time -= $hours*$hour;
    
    my $minutes = int($time/$minute);
    $time -= $minutes*$minute;
    
    my $seconds = $time;
    
    my %formatted = (
        'weeks' => $weeks,
        'days' => $days,
        'hours' => $hours,
        'minutes' => $minutes,
        'seconds' => $seconds
    );
    
    return \%formatted;
}

1;