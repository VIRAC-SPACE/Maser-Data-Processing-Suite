#!/usr/bin/perl
use POSIX;
use File::Copy;
use List::Util qw( min max );
use Data::Dumper qw(Dumper);
 
my @source = @ARGV[0];
my @outputDir = @ARGV[1];

my $directory = "@source/@outputDir";

open(DB,">DB.csv");
open(PYTHON,">python.txt");
print DB "name,min_date,max_date,JD_min,JD_max,days,min_v,max_v,n_spectra\n";

 opendir( my $DIR_list, $directory );
#while ( my $dir_entry = readdir $DIR_list ) {
    #next unless -d $directory . '/' . $dir_entry;
    #next if $dir_entry eq '.' or $dir_entry eq '..';
    #print "Found directory $dir_entry\n";

	open(F,">$directory/3d.txt");
	open(list_days,">$directory/days.txt");
	open(labels,">$directory/labels.txt");

    opendir (DIR, "$directory/") or die $!;

	my @dates = (); undef @dates;
	my @jds = (); undef @jds;
	my @diffs = (); undef @diffs;
	my @files = (); undef @files;

	  while (my $file = readdir(DIR)) {
		  next if not -f "$directory/$dir_entry/".$file;
		  next if not $file =~ /.dat$/;

        #print "$file\n";
		my @array = split('_',$file);
		my $name = $array[0];
		my $h = $array[1];
		my $m = $array[2];
		my $s = $array[3];
		my $day = $array[4];
		my $month = $array[5];
		my $year = $array[6];

		#print "$name $h $m $s $day $month $year\n";
		

		my @months = qw(jan feb mar apr may jun jul aug sep oct nov dec);
		my $ind = indexArray(lc $month,@months);
		if ($ind >=0 ) {
			$ind++;
			$ind = "0$ind" if $ind < 10;
			my $jd = JulianDate("$year-$ind-$day");
			$jd = sprintf("%.4f",$jd + ($h+$m/60.0)/24.0);
			#print "file=$file, jd=$jd\n";
			#copy("$directory/$file","data_txt/$jd-$year-$ind-$day.txt") or die "Copy failed: $!"; 
			#print "$directory/$file data_txt/$year-$ind-$day-$jd.txt\n";
			push(@dates,"$year-$ind-$day $h:$m");
			push(@jds,"$jd");
			push(@files,"$file");

		} else {
			print "WARNING: file is in wrong format, month $month unknown\n";
		}

    }

	  closedir(DIR);

my $n_spectra = @jds;
my $jd_first = min @jds;
my $jd_last = max @jds;

#print "First JD is $jd_first\n\n";
print list_days "First JD is $jd_first\n";
#print labels "First JD is $jd_first\n";


foreach my $jd (@jds) {
	push(@diffs,sprintf("%.4f",$jd-$jd_first));
}

my @idx = sort { $jds[$a] <=> $jds[$b] } 0 .. $#jds;

@jds = @jds[@idx];
@diffs = @diffs[@idx];
@dates = @dates[@idx];
@files = @files[@idx];


# Find gaps
my $dif_pre;
my $j = 0;
foreach my $dif (@diffs) {

	#if ($dif_pre and ($dif-$dif_pre > 6)) {
		#print "Big difference from $dates[$j-1] to $dates[$j]: ".($dif-$dif_pre)." days\n";
		#print labels "$dif_pre,$dif\n";
	#}
	$dif_pre = $dif;
	$j++;
}

my $count = @diffs;
my $nth = 1;
$nth = 2 if $count > 20;
$nth = 3 if $count > 30;
$nth = 4 if $count > 40;
$nth = 5 if $count > 50;
$nth = 6 if $count > 60;
$nth = 7 if $count > 70;
$nth = 8 if $count > 80;

my @diffs_short = @diffs[grep { ! (($_+1) % $nth) } 0..$#diffs]; # Pick every nth element
my @dates_short = @dates[grep { ! (($_+1) % $nth) } 0..$#dates];
my $list_diffs_short = join(",",@diffs_short);
my $list_dates_short = join(",",@dates_short);
my $list_diffs = join(",",@diffs);

print labels $list_diffs;
print labels "\n";
print labels $list_dates_short;
print labels "\n";
print labels $list_diffs_short;



my $min_v = 999;
my $max_v = -999;

for (my $i = 0; $i < $#jds ; $i++) 
{
	#print "File $files[$i] with date diff $diffs[$i]\n";
	print list_days "File $files[$i], diff: $diffs[$i]\n";
		open(D,"$directory/$dir_entry/".$files[$i]);
		my ($vel,$flux);

		while (<D>) {
			chomp $_;
			my @split_arr = split(/\s/,$_);
			$vel = sprintf("%.3f",$split_arr[0]);
			$min_v = $vel if $vel < $min_v;
			$max_v = $vel if $vel > $max_v;
			$flux = sprintf("%.3f",$split_arr[3]);
			# if ($vel > 55.0 and $vel < 65.0) {
				print F "$vel $flux $diffs[$i]\n";
			# }
			
		}	
		print list_day "$files[$i]  $diffs[$i]  $jds[$i]\n";	
}

#print "min_v = $min_v\n";
#print "max_v = $max_v\n";

close(labels);
close(F);
close(list_days);

my $days = sprintf("%.1f",$jds[$n_spectra-2] - $jds[0]);
print DB "$dir_entry,$dates[0],$dates[$n_spectra-1],$jds[0],$jds[$n_spectra-1],$days,$min_v,$max_v,$n_spectra\n";

if ($min_v == 999 or $max_v == -999) {
	$min_v = -20;
	$max_v = 20;
}

#my $cmd = "python3 code/plot_multiple.py $directory/3d.txt $directory/labels.txt @source.png 0 $days $min_v $max_v $jd_first @source";
my $cmd = "{\n\"days\":$days,\n\"min_v\":$min_v,\n\"max_v\":$max_v,\n\"jd_first\":$jd_first\n}";
print PYTHON $cmd;
#system($cmd);


#}

close(DB);
close(PYTHON);

sub indexArray($@)
{
 my $s=shift;
 $_ eq $s && return @_ while $_=pop;
 -1;
}


sub JulianDate {

	my $date = $_[0];

	my ($yr,$mn,$mday) = split('-',$date);
	my $hr = 0;

	# In leap years, -1 for Jan, Feb, else 0
	my $L = POSIX::ceil(($mn - 14) / 12);

	my $p1 = $mday - 32075 + POSIX::floor (1461 * ($yr + 4800 + $L) / 4);
	my $p2 = POSIX::floor (367 * ($mn - 2 - $L * 12) / 12);
	my $p3 = 3 * POSIX::floor (POSIX::floor (($yr + 4900 + $L) / 100) / 4);
	my $julian = $p1 + $p2 - $p3;        
	$julian = $julian + ($hr / 24.0) - 0.5;

	my $tjd = $julian - 2440000.5;

	return $julian;
}
