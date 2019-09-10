($dir, $nchan) = @ARGV;

@files = split "\n", `cat currents.txt`;
@IPD = split ' ', `cat leds.txt`;

$j_use_for_beta = 7;

@color = ("#000000", "#600000", "#ff0000", "#a06000", "#30c000", "#00c070", "#0000ff", "#200090", "#808080",
  "#000000", "#600000", "#ff0000", "#a06000");

@color = ("#808080", "#404040", "#000000", "#600000", "#ff0000", "#ff6000", "#a06000", "#30c000",
   "#00c070", "#0000ff", "#200090", "#800080", "#ff00ff");

for $ichan (1..$nchan) {
  $C1 = 2*$ichan-1;
  $C2 = $C1+1;
  open(G, "| tee temp.script | gnuplot") or die;
  print G qq|set term postscript enhanced 22 color eps\n|;
  print G qq|set output "$dir/crnl-$ichan.eps"\n|;
  print G qq|set xlabel "log_{10} Signal level [DN]"\n|;
  print G qq|set ylabel "ln (observed signal/expected signal) + offset"\n|;
  print G qq|set size 1.8,1.8\n|;
  print G qq|set title "Channel $ichan"\n|;
  print G qq|set xrange [1.5:4.5]; set xtics .25; set yrange [-.035:.035]; set ytics .005\n|;
  print G qq|set grid\n|;
  for $j (0..((scalar @files)-1)) {
    $jj = $j+1;
    print G qq|set style line $jj lt 1 lw .5 pt 1 ps .5 lc rgb "$color[$j]"\n|;
  }
  print G qq|plot|;

  $bg_use = 0;
  if(1) {
    ($f, $dt, $ledsline) = split ' ', $files[$j_use_for_beta];
    $rampdata = (split "\n", `cat $f\_offsets.txt`)[$ichan-1];
    ($offset, $offset_final, $bg_use) = @coefs = split ' ', $rampdata;
    #$bg_use = 0;
    @nlc = split ' ', `head -n $ichan nl.txt | tail -n 1`;
    $dnl = 1 + scalar @nlc;
    print 'NL: ', @nlc, " [$dnl]\n";

    $nlines = (split ' ', `wc $f\_mean.txt`)[0];
    @lastline = split ' ', `tail -n 1 $f\_mean.txt`;
    $voffset = log($lastline[$C1-1]/($nlines-1)/$dt/$IPD[$ledsline * 3 + 1]);
    if (0) {
      $voffset -= log(1. - $coefs[3]*($offset*2+$lastline[$C1-1]));
    } else {
      $voffset -= log(1. - $bg_use*($offset*2+$lastline[$C1-1]));
    }
  }
  print "channel $ichan --> beta*g used = $bg_use DN^{-1}; voffset = $voffset\n";

  for $j (0..((scalar @files)-1)) {
    $jj=$j+1;
    ($f, $dt, $ledsline) = split ' ', $files[$j];
    $k = $ledsline * 3 + 1;

    $rampdata = (split "\n", `cat $f\_offsets.txt`)[$ichan-1];
    ($offset, $bg) = split ' ', $rampdata;

    $title = sprintf "%9.3E A", $IPD[$k];
    #print G qq| "$f\_mean.txt" using (log10(\$$C1)):(log(\$$C1/\$0/$dt/$IPD[$k]/(1.-($bg_use)*($offset*2+\$$C1)))-$voffset):(\$$C2/\$$C1) with yerrorbars title "$title" ls $jj|;
    print G qq| "$f\_mean.txt" using (log10(\$$C1)):(log(\$$C1/\$0/$dt/$IPD[$k]*(1.|;
    for $o (2..$dnl) {
      $nc = $nlc[$o-2];
      print G qq|+($nc)*(($offset+\$$C1)**$o - ($offset)**$o)/(\$$C1)|;
    }
    print G qq|))-$voffset):(\$$C2/\$$C1) with yerrorbars title "$title" ls $jj|;
    if ($j<(scalar @files)-1) {print G qq|,|;}
  }
  print G qq|\n|;
  close G;
}
