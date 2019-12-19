# Script to run several configuration files in parallel.
# Useful for running on a multi-core machine.
#
# Example: perl multirun.pl config.1 config.2
# should run both config.1 and config.2 (can use for N files: config.1 ... config.N)
#
# Make sure the output directories don't collide!
# (This script doesn't check.)

@info = @ARGV;

$N = scalar(@info);

print "$N runs\n";

for $i (0..$N-1) {
  if (-e "tempresults-$info[$i]") {
    print "Error: need to write to tempresults-$info[$i], but file already exists.\n";
    exit;
  }
}

my $i=0;
for $i (0..$N-1) {
  my $pid=fork;
  if (not defined $pid) {
    print STDERR "Error: fork failed\n";
    exit;
  }
  if (not $pid) {
    system "date > tempresults-$info[$i]";
    system "python test_run.py $info[$i] >> tempresults-$info[$i]";
    system "date >> tempresults-$info[$i]";
    exit;
  }
}

# Wait for children
my $k;
for $k (1..$N) {wait();}

for $i (0..$N-1) {
  print "=== Results from configuration file $i/$N -> $info[$i] ===\n";
  system "cat tempresults-$info[$i]";
  print "\n";
  system "rm tempresults-$info[$i]";
}
