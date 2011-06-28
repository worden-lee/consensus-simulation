#!/usr/bin/perl

my $usage =
"Usage:  ./run-batch2.pl [options] exptname
options:
    -r | remove all data files, leave outcomes, statistics, and graphs
    -d | display debugging information
    -c | collect outcomes only, don't run simulations
    -s | make statistics from outcomes only
    -D dir  | directory for data sets
    -O dir  | directory for outcomes.dat
    -S dir  | directory for statistics.dat
    -f file | settings file to be passed on to simulation program
    --x=y   | argument to be passed on to simulation program
    --compare=A,B  | run with settings file A and with settings file B
    --vary X=A,B,C | run with --X=A, with --X=B, and with --X=C
    --vary X=A-Z   | run with --X=Y for all Y in A, A+1, ..., Z
    --vary X=A-e-Z | run with --X=Y for all Y in A, A+e, ..., Z
    --executable=E | the program executable is E
    perlexpr=value | perl code to evaluate in run-batch.pl
";

sub usage_err {
  print $usage;
  exit;
}

# package for mkpath
use File::Path;

use Tie::Handle::CSV;

$[ = 0;

# independent vars and how they loop
my @independent_vars, @independent_vars_values;
my $includetxt = "";
my @evals;
my $nreps = "";
$intro = 0;

# dependent vars
my @dependent_vars = ("n.individuals", "min.value", "mean.value", "max.value",
  "n.satisfied");
# indexes of equil, ess, final_ess in the above list
#my @dep_vars_to_plot = (0, 2, 5);

my $pcomm;

# default values
$dbg = 0;
$executable = './consensus';
$delete = 0;
$runsims = 1;
$collect = 1;
$dostats = 1;
$makeplots = 1;
$customvalues = "";
@settings = ();
#$exptname = "random-community";

# process command line
while (@ARGV)
{
  my $arg = shift ARGV;
  if ($arg =~ /^-r$/) {
    $delete = 1;
    print "Received flag $arg: removing output directories\n";
  }
  elsif ($arg =~ /^-d$/)
  { $dbg=1;
    $dbg && print "Saw -d option\n";
  }
  elsif ($arg =~ /^-c$/)
  { $runsims=0;
    $dbg && print "Saw -c option\n";
  }
  elsif ($arg =~ /^-s$/)
  { $runsims=0;
    $collect=0;
    $dbg && print "Saw -s option\n";
  }
  elsif ($arg =~ /^-D$/)
  { $data_dir = shift ARGV; }
  elsif ($arg =~ /^-O$/)
  { $outcomes_dir = shift ARGV; }
  elsif ($arg =~ /^-S$/)
  { $statistics_dir = shift ARGV; }
  elsif ($arg =~ /^-f$/)
  { push @settings, shift ARGV; }
  elsif ($arg =~ /^--compare=(.*)$/)
  { push @independent_vars, 'include';
    push @independent_vars_values, [ split(/,/,$1) ];
  }
  elsif ($arg =~ /^--executable=(.*)$/)
  { $executable = $1; }
  elsif ($arg =~ /^--vary$/)
  { $arg = shift ARGV;
    if ($arg =~ /^(.*)=(.*)-(.*)-(.*)$/)
    { push @independent_vars, $1;
      my @vals = map { $2 + $3*$_ } (0..(($4-$2)/$3));
      push @independent_vars_values, [ @vals ];
    }
    elsif ($arg =~ /^(.*)=(.*)-(.*)$/)
    { push @independent_vars, $1;
      push @independent_vars_values, [$2..$3];
    }
    elsif ($arg =~ /^(.*)=(.*)$/)
    { push @independent_vars, $1;
      my @vals = split(/,/,$2);
      push @independent_vars_values, [ @vals ];
    }
    else
    { &usage_err(); }
  }
  elsif ($arg =~ /^--(.*)=(.*)$/)
  { $customvalues .= "$1 $2\n";
  }
  elsif ($arg =~ /=/)
  { push @evals, $arg; }
  elsif ($arg =~ /^-/)
  { &usage_err(); }
  else
  {
    $exptname = $arg;
    print "replaced experiment name with '$arg'\n";
  }
}

# eval the commandline specifications before and after
for (@evals)
{ eval $_; }

#$dbg && print "experiment name is $exptname\n";

if (0) {
#
# evaluate the experiment name
#
# these are more or less ordered from old to new, some of the old
# ones may not work any more
#
if ($exptname eq "one") {
  # 2 resources, 1 species
  @independent_vars = ("replicate", "resources[0].heating",
		       "resources[1].heating", "species[0].tau" );
  @mins =  (1, 0, 0, 0);
  @maxs =  (2, 1, 1, 100);
  @steps = (2, 1, 1, 2);  # number of steps, counting both min and max
#  @steps = (3, 3, 5, 2);
  $includetxt .= "include experiments/$exptname.settings\n";
}
}

# eval the commandline specifications before and after
for (@evals) {
  eval $_;
}

my @independent_vars_norep; # like @independent_vars but without replicate
my $nkeys = $#independent_vars + 1; # including replicate
my $replicate_column;
foreach $i (0 .. $#independent_vars)
{ if ($independent_vars[$i] ne "replicate")
  { push @independent_vars_norep, $independent_vars[$i]; }
  else
  { $replicate_column = $i; }
}


if (!defined $data_dir){
  $data_dir = "batch-out/$exptname";

  # special case -- if running a new simulation from scratch make sure
  # it's a fresh directory
  if ($run_sims && -e $data_dir) {
    for ($i=1; $i<10000; $i++) {
      $try = $data_dir . "-$i";
      if (!-e $try) {
        $data_dir = $try;
        last; # i=100;
      }
    }
  }
}
$data_dir =~ s|/*$||;

if (!defined $outcomes_dir) {
  $outcomes_dir = $data_dir;
}
$outcomes_dir =~ s|/*$||;

if (!defined $statistics_dir) {
  $statistics_dir = $outcomes_dir;
}
$statistics_dir =~ s|/*$||;

if ($runsims)
{
  $settingsfile = "batch.settings";
  $batchlogfile = "batch.log";

  # delete files if we have -r activated
  if ($delete) {
    #print "Deleting files!\n";
    $comm = "rm -R $data_dir/*";
    print "$comm\n";
    system("$comm") || die "couldn't delete $data_dir/*";
  }

  if (!-e $data_dir)
  {
    print "Create $data_dir \n";
    mkpath($data_dir) || die "couldn't create $data_dir/";
  }

  open VARS, ">$data_dir/vars"
    or die "couldn't open $data_dir/vars";
  print VARS join(" ", @independent_vars), "\n";
  close VARS;

  my $maxcounter = 1;
  for (@independent_vars_values)
  { $maxcounter *= scalar @$_; }
  $dbg && print "maxcounter is $maxcounter\n";
  my $counter = 0;
  while ($counter < $maxcounter)
  {
    my $c = $counter;
    my @vals;
    for my $i (reverse (0 .. $#independent_vars))
    { my $ni = scalar @{$independent_vars_values[$i]};
      $vals[$i] = $independent_vars_values[$i][$c % $ni];
      $c /= $ni;
    }
    print "@vals\n";

    my $tmpdir = "$data_dir/running";
    if (!-e $tmpdir)
    { mkpath($tmpdir); }
    system("rm -rf $tmpdir/*");

    open SETTINGS, ">$tmpdir/$settingsfile"
      or die "couldn't open $tmpdir/$settingsfile";
    #print SETTINGS "include settings/fast.settings\n";
    #print SETTINGS "include settings/precise.settings\n";
    print SETTINGS map { "include $_\n"; } @settings;
    print SETTINGS $includetxt;

    print SETTINGS "disableDisplaying true\n";
#    print SETTINGS "displayFlowsInNetwork true\n";

    print SETTINGS $customvalues;

    print SETTINGS "outputDirectory $tmpdir\n";

#    system("rm -rf out");

    for (0 .. $#independent_vars)
    { print SETTINGS "$independent_vars[$_] $vals[$_]\n"; }

    close SETTINGS;

    my $comm = "$executable -f $tmpdir/$settingsfile 2>&1 | /usr/bin/tee $tmpdir/$batchlogfile";
    print "$comm\n";
    my $t0 = time;
    if (system($comm))
    { if ($? == -1) {
        print "failed to execute: $!\n";
      }
      elsif ($? & 127) {
        my $sig  = ($? & 127);
        printf "child died with signal %d, %s coredump\n",
          $sig,  ($? & 128) ? 'with' : 'without';
        if ($sig == 2) {
          &cleanup;
          die "terminating";
        }
      }
      else {
        my $ret = $? >> 8;
        printf "child exited with value %d\n", $ret;
        if ($ret == 255)
        { &cleanup;
          die "terminating";
        }
      }
      &cleanup;
      die "simulation returned nonzero";
    }
    my $t1 = time;
    my $runtime = $t1 - $t0;

    open LOG, ">>$tmpdir/log.0"
      or die "can't append to $tmpdir/log.0";
    print LOG "runtime: $runtime sec\n";
    close LOG;
    print "runtime: $runtime sec\n";

    #calculate here whether we've finished a replicate set, if so, record?
    # this would be where $i==$#steps and $c%$steps[$i]==0

#     $comm = "mv $settingsfile out/";
#     print "$comm\n";
#     system("$comm") && die "couldn't relocate $settingsfile";

#     $comm = "mv $batchlogfile out/";
#     print "$comm\n";
#     system("$comm") && die "couldn't relocate $batchlogfile";

#     $comm = "mv out/ $dir";
#     print "$comm\n";
#     system("$comm") && die "couldn't relocate out/";

    my $outdir = "$data_dir/".join('/', @vals);
    #print "\n$dir\n";
    if (!-e $outdir) # make the directory AND all its parent dirs
    { mkpath($outdir); }
    system("rm -rf $outdir");
    $comm = "mv $tmpdir/ $outdir";
    print "$comm\n";
    system("$comm") && die "couldn't relocate $tmpdir to $outdir";

    if($delete && $vals[$#vals] == $maxs[$#maxs]) {
      if ($collect) {
        $dbg && print "About to record results!\n";
        record_results();
        $dbg && print "Done recording results!\n";
      }

      $dbg && print "Deleting files!\n";
      $comm = "rm -R $data_dir/$vals[0]";
      print "$comm\n";
      system("$comm") && die "couldn't delete $data_dir/$vals[0]";
    }

    ++$counter;
  }
}

if ($collect && !$delete) {
  $dbg && print "About to record results!\n";
  record_results();
  $dbg && print "Done recording results!\n\n";
}
if ($dostats) {
  $dbg && print "About to record statistics!\n";
  record_statistics();
  $dbg && print "Done recording statistics!\n\n";
}

if ($makeplots) {
  $dbg && print "About to make plots!\n";
  make_plots();
  $dbg && print "Done making plots!\n\n";
}

sub sort_by_columns
{
  my($a,$b) = @_;
  return 0 if $a !~ /^\s*(\S+)\b(.*)$/;
  my($a0,$a1) = ($1,$2);
  return 0 if $b !~ /^\s*(\S+)\b(.*)$/;
  my($b0,$b1) = ($1,$2);
  my $comp = ($a0 <=> $b0);
  return $comp if ($comp != 0);
  $comp = ($a0 cmp $b0);
  return $comp if ($comp != 0);
  return sort_by_columns($a1,$b1);
}

sub key_sort_sub
{
  sort_by_columns($a,$b);
}

sub compare_arrays {
  my ($first, $second) = @_;
  my @a1 = @{$first};
  my @a2 = @{$second};
  if (scalar(@a1) ne scalar(@a2)) {
    return 0;
  }
  my $n = 0;  #index counter
  while ($a1[$n]) {
    if ($a1[$n] ne $a2[$n]) {
      return 0;
    }
    $n++;
  }
  return 1
}

sub record_results {

  if (!-e $outcomes_dir)
  { mkpath($outcomes_dir); }

  my $outcomes_file = "$outcomes_dir/outcomes.out";

  my %data;
  for $csv (`find $data_dir -name outcome.csv -print`)
  { chomp $csv;
    my $key = $csv;
    $key =~ s|$data_dir/(.*)/outcome.csv|$1|;
    $key =~ s|/| |g;
    $dbg && print "opening $csv - key is $key\n";
    my %outcome;
    my $fh = Tie::Handle::CSV->new($csv, header=>1);
    while (my $csv_line = <$fh>) # expect only one line
    { my @names = keys %$csv_line;
      $dbg && print "names: ". join(' | ',@names). "\n";
      $dbg && print "values: " . join(' ', @$csv_line{@names})."\n";
      @outcome{@names} = @$csv_line{@names};
    }
    $data{$key} = join(' ',@outcome{@dependent_vars});
    print "$key $data{$key}\n";
  }

  if (0) {
  for $dir (`find $data_dir -name log.0 -print`)
  { my %current_community = (), %last_community = ();

    chomp $dir;
    $dir =~ s|/log.0||;
    my $logfile = "$dir/log.0"; # need to check log.1 etc.?
    $dbg && print "reading $logfile\n";
    open LOGFILE, "$logfile" or die "couldn't open $logfile";

    my %outcome;
    @outcome{@dependent_vars} = (0) x scalar(@dependent_vars);

    # list of keys is in the name of the directory
    my $td = $dir;
    $td =~ s|$data_dir/||;
    my $key = $td; $key =~ s|/| |g;
    $dbg && print "key: $key\n";

    while (<LOGFILE>)
    {
      if (/[0-9] equilibrium/)
      { $outcome{"equil"} = 1;
      }
      #if (/equilibrium/) # includes 'tired of waiting'
      if (/^([0-9.e+-]+) equilibrium/)
      { $read_after_equil = 1;
	$dbg && print "at equilibrium t=$1: ";
      }
      if (/ESS/) { $outcome{ess} = 1; ++$outcome{ess_count}; }
  #    if (/No more species/) { $ess_nsp = 0; }
      if (/T -> ([\d\.e+-]+)/) { $outcome{final_T} = $1; }
      if (/couldn't find a viable new species type/)
      { $outcome{final_ess} = 1; }
      if ($outcome{final_ess} && /R\d+ ->/) { $outcome{final_nr}++; }
      if (/lineages -> (\d+)/)
      { my $nl = $1;
        if (($outcome{eq_nsp} == 0) && $outcome{equil})
        { $outcome{eq_nsp} = $nl; }
        if (($outcome{ess_nsp} == 0) && $outcome{ess})
        { $outcome{ess_nsp} = $nl; }
        if ($outcome{final_ess}) { $outcome{final_nsp} = $nl; }
      }
      if (/runtime: (\d+) sec/)
      { $outcome{runtime} = $1; }
      if (/speciation \((N[0-9]+) => (N[0-9]+)\)/)
      { $lineage{$2} = $lineage{$1}; }
      if (/^new species type/)
      { $new_species_type = 1; }
      elsif ($new_species_type)
      { if (/^(N[0-9]+) /)
        { $lineage{$1} = $1; }
        else
        { print "New species type not found ($logfile line $.)!\n"; }
        $new_species_type = 0;
      }
      if ($read_after_equil && /(N[0-9]+) ->/)
      { ++$current_community{$lineage{$1}};
        $dbg && print "$1:",$lineage{$1}," ";
      }
      if ($read_after_equil && /^\}/)
      { $read_after_equil = 0;
#	$dbg &&
#	  print "\n", join(" ",sort keys %current_community), " =? ",
#	    join(" ",sort keys %last_community), "\n";
        if (join(" ",sort keys %current_community)
            ne join(" ",sort keys %last_community))
        { ++$outcome{community_structures};
          %last_community = %current_community;
          $dbg && print " DIFFERENT COMMUNITY";
        }
        $dbg && print "\n";
        %current_community = ();
      }
    }
    close LOGFILE;

#    $dbg && print "outcome{equil} = $outcome{\"equil\"}\n";
#    $dbg && print "outcome{@dependent_vars} = @outcome{@dependent_vars}\n";
    $data{$key} = join(' ',@outcome{@dependent_vars});
    print "$key $data{$key}\n";
  }
  }

  my @keys = sort key_sort_sub keys %data;
  print "\n";

  if($intro==0) {
    open OUTCOME, ">$outcomes_file"
      or die "couldn't open $outcomes_file for writing";
    #print OUTCOME "# $independent_vars $dependent_vars\n";
    print OUTCOME "# ",
      join(' ',(@independent_vars,@dependent_vars)),
      "\n";
    $intro=1;
  }
  else {
    open OUTCOME, ">>$outcome"
      or die "couldn't open $outcome for writing";
  }

  my $lasth0 = "";
  for $k (@keys)
  {
    if (($k =~ /(\S+)\s/) and ($1 ne $lasth0))
    {
      if ($lasth0 ne "")
      {
        print "\n";
        print OUTCOME "\n";
      }
      $lasth0 = $1;
    }
    print "$k $data{$k}\n";
    print OUTCOME "$k $data{$k}\n";
  }
  close OUTCOME;
}

sub record_statistics {

  my $outcome = "$outcomes_dir/outcomes.out";

  if (!-e $statistics_dir)
  { mkpath( $statistics_dir ); }

  my $statistics = "$statistics_dir/statistics.out";
  open STATISTICS, ">$statistics"
    or die "couldn't open $statistics for writing";

  open READFILE, "$outcome" or die "Couldn't open $outcome";

  # read the header line and don't use it, use internal data
  $line = <READFILE>;
  chomp $line;

  $dbg && print "\n\nDebugging Statistics\n============================\n\n";
  #my $sigmaNumber = -1;
  my $count = 0;

  $dbg && print "keys: @independent_vars_norep\n";
  $dbg && print "replicate_column: $replicate_column\n";
  $dbg && print "dependent vars: @dependent_vars\n";

  print STATISTICS
    join(' ',(@independent_vars_norep,"replicates",@dependent_vars)),"\n";

  # now read each line of outcomes and collect in %stats

  while ($line = <READFILE>)
  {
    chomp $line;
    next if ($line =~ /^\s*$/);

    my @fields = split(/ +/, $line);
    my @keys;
    foreach $i (0 .. $nkeys-1)
    { if ($i != $replicate_column)
      { push @keys, $fields[$i]; }
    }
    my $replicate = $fields[$replicate_column];
    my @values = @fields[$nkeys .. ($#fields)];

    my $key_string = join(' ',@keys);

    $dbg && print "====================\n\n";
    $dbg && print "Read line: @fields\n";
    $dbg && print "  keys: @keys\n";
    $dbg && print "  replicate: $replicate\n";
    $dbg && print "  values: @values\n";

    # $key_string is the independent variables' values
    #  (not counting replicate) joined by ' '
    # $stats{$key_string} is (by reference) an array containing
    #  first the # of replicates, and then the sums of each of the
    #  dependent variables.
    if (!defined($stats{$key_string}))
    { $stats{$key_string} = [1, @values];
#      $dbg &&
#	print "created stats{$key_string}: @{$stats{$key_string}}\n";
    }
    else
    { ++$stats{$key_string}[0];
      foreach $i (0 .. $#dependent_vars)
      {
#	$dbg && print "stats{$key_string}[",$i+1,"] == $stats{$key_string}[$i+1]; ";
#	$dbg && print "+= $values[$i]: ";
	$stats{$key_string}[$i+1] += $values[$i];
#	$dbg && print "$stats{$key_string}[$i+1]\n";
      }
    }
    $dbg && print "  Current totals (@keys): @{$stats{$key_string}}\n\n";
  }

  $dbg && print "Done! statistics: \n";
  foreach $k (sort keys %stats)
  { my $nreps = $stats{$k}[0];
    print STATISTICS "$k $nreps";
    $dbg && print "$k $nreps";
    foreach $i (0 .. $#dependent_vars)
    { print STATISTICS " ",($stats{$k}[$i+1]/$nreps);
      $dbg && print " ",($stats{$k}[$i+1]/$nreps);
    }
    print STATISTICS "\n";
    $dbg && print "\n";
  }

  close READFILE;
  close STATISTICS;
}

sub make_plots {
  # make plots of each independent_vars -- except replicate -- vs
  #    equil, ess, final_ess
  foreach $i (0 .. $#independent_vars_norep) {
    foreach $d (@dep_vars_to_plot) {
      my ($ind, $dep)
	= ($independent_vars_norep[$i], $dependent_vars[$d]);

      $dbg && print "\nPlotting $dep vs $ind\n";
      my $title = $dep."_vs_".$ind.".gp";
      my $depcol = $d + $#independent_vars_norep + 2;
      # the column numbers used here are in statistics.out, and starting
      # from 1
      my $pcomm = "../bin/plot-cols -x ".($i+1)." ".($depcol+1).
	" $statistics_dir/statistics.out ".($dbg?"-d ":"").
	  "-nokey -t \"$dep vs. $ind\" -w p -g $statistics_dir/$title";
      $dbg && print "Plotting data: $pcomm\n\n";
      system("$pcomm") == 0
	or die "couldn't run graphing script!\n";
    }
  }
}

sub cleanup {
  unlink($settingsfile);
  unlink($batchlogfile);
}

