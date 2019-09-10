($file) = @ARGV;

print STDERR "<-- $file\n";

open(IN, $file) or die;
$head = <IN>;
$count = 0;
for $fields (split ' ', $head) {
  if ($fields =~ m/^[\d\.Ee]+$/) {
    $current[$count] = $fields;
    $count++;
  }
}
#print "Currents: @current ($count values)\n";
$N=$count;
$Nrow=0;
while ($line=<IN>) {
  if ($line!~m/^\#/) {
    @data = split ' ', $line;
    for $i (0..$N-1) {$ID[$i][$Nrow]=$data[$i+1]-($data[0]+$data[$N+1])/2.;}
    $Nrow++;
  }
}
close IN;

#print "$Nrow rows of data\n";

# Print outputs
for $i (0..$N-1) {
  $sum=$sum2=0;
  for $j (0..$Nrow-1) {
    $sum += $ID[$i][$j];
    $sum2 += $ID[$i][$j]**2;
  }
  $sum/=$Nrow;
  $sum2/=$Nrow;
  $stdev = sqrt(($sum2-$sum**2)/($Nrow-1));
  print (sprintf "%11.5E %11.5E %11.5E\n", $current[$i], $sum, $stdev/$sum);
}
