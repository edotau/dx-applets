#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2007-2011 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software 
License Agreement" and the "Illumina Source Code License Agreement", 
and certain third party copyright/licenses, and any user of this 
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party 
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=head1 NAME

configureValidation.pl - Create and configure the test harness.

=head1 SYNOPSIS

configureValidation.pl [options]

Options:

=over 4

=item B<--generate> Produce the chain of validation datasets instead of comparing

Z<>

=item B<--help> brief help message

Z<>

=item B<--man> full documentation

Z<>

=item B<--source-dir=<path>>

Overrides default location of validation data and configuration. 

 Default is: /usr/local/share/bcl2fastq-1.8.4/examples/Validation/Default

=item B<--base-calls-dir=<path>>

Overrides the default for RTA base calls folder.

 Default is: /usr/local/share/bcl2fastq-1.8.4/examples/Validation/Default/Data/Intensities/BaseCalls

=item B<--output-dir=<path>>

path to the output directory.

 Default is current path with name derived from the --source-dir folder name

=item B<--sample-sheet=<path>>

Overrides the default sample sheet path.

 Default is: /usr/local/share/bcl2fastq-1.8.4/examples/Validation/Default/Data/Intensities/BaseCalls/SampleSheet.csv

If the default sample sheet is not present, the dataset is considered non-multiplexed

=item B<Running end-to-end validation on the default test data set included with bcl2fastq:>

 /usr/local/bin/configureValidation.pl --output-dir ./ValidationDefault

 make -C ValidationDefault  # use make -j <parallel jobs> to speedup the conversion

=item B<Running end-to-end validation on the bigger test data set included with bcl2fastq (contains multiplexed RNA data):>

 /usr/local/bin/configureValidation.pl \
 --source-dir /usr/local/share/bcl2fastq-1.8.4/examples/Validation/110120_P20_0993_A805CKABXX \
 --output-dir ./ValidationMultiplexed

 make -C ValidationMultiplexed  # use make -j <parallel jobs> to speedup the conversion

=back

=head1 DIAGNOSTICS

=head2 Exit status

=over 4

=item B<0:> successful completion

=item B<1:> abnormal completion

=item B<2:> fatal error

Z<>

=item B<Errors:> All error messages are prefixed with "ERROR: ".

Z<>

=item B<Warnings:> All warning messages generated by CASAVA are prefixed with "WARNING: ".

Z<>

=back

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Roman Petrovski

=cut

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install

use warnings "all";
use strict;

use File::Path qw(mkpath);
use Getopt::Long;
use Pod::Usage;

use Casava::Common::Log qw(initLog logInfo errorExit);
use Casava::Common::Utils qw(reallyRealPath);

$ENV{PATH} = $1 if $ENV{PATH}=~/^(.*)$/;

initLog( undef, 3, 0, 1, 1);

# load the default configuration parameters before parsing the command line
my $CASAVA_DATADIR = '/usr/local/share/bcl2fastq-1.8.4'; # substituted during the install

my $man = 0;
my $help = 0;
my $make = 0;
my $force = 0;
my $generate = 0;
my $outputDir;
my $validationRunFolder = File::Spec->catfile( $CASAVA_DATADIR, 'examples', 'Validation', 'Default');
my $sampleSheetFile;
my $baseCallsDir;

my $result = GetOptions(
                     'source-dir=s' => \$validationRunFolder,
                     'base-calls-dir=s' => \$baseCallsDir,
                     'output-dir=s' => \$outputDir,
                     'sample-sheet=s' => \$sampleSheetFile,
                     'generate!' => \$generate,
                     'force!' => \$force,
                      help => \$help,
                      man => \$man,
                     ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2,  -input => $1) if ($man and $0 =~ /(.*)/);

errorExit("ERROR: Can't obtain absolute path for source-dir: '$validationRunFolder'.")
    if reallyRealPath($validationRunFolder);

$baseCallsDir = File::Spec->catfile($validationRunFolder, 'Data', 'Intensities', 'BaseCalls') unless $baseCallsDir;
errorExit("ERROR: Can't obtain absolute path for BaseCalls path: '$baseCallsDir'.")
    if reallyRealPath($baseCallsDir);

$sampleSheetFile = File::Spec->catfile( $baseCallsDir, 'SampleSheet.csv') unless $sampleSheetFile;

if (!$outputDir)
{
    my ($vol, $path, $file) = File::Spec->splitpath($validationRunFolder);
    $outputDir = $generate ? "./$file.generation" : "./$file.validation";
    reallyRealPath($outputDir);
    logInfo "--output-dir is unspecified. Will configure the validation in: $outputDir\n";
}
else
{
    errorExit("ERROR: Can't obtain absolute path for --output-dir directory argument: '$outputDir'.")
        if reallyRealPath($outputDir);
}

if (-e $outputDir)
{
    errorExit "ERROR: $outputDir already exists. Please specify a different path." unless $force;
}
else
{
    logInfo "Creating directory '$outputDir'";
    File::Path::mkpath( $outputDir );
    errorExit "ERROR: failed to create directory '$outputDir'" unless -d $outputDir;
}

my @configuration;
push @configuration, "BASE_CALLS_ROOT:=$baseCallsDir";
push @configuration, "ALIGNED_ROOT:=" . File::Spec->catfile( $validationRunFolder, 'Aligned');
push @configuration, "UNALIGNED_ROOT:=" . File::Spec->catfile( $validationRunFolder, 'Unaligned');
my $demuxConfigMk = File::Spec->catfile( $validationRunFolder, 'Unaligned', 'ValidationConfig.mk');
push @configuration, "DEMUX_CONFIG_MK:=" . $demuxConfigMk if -e $demuxConfigMk;
push @configuration, "SAMPLE_SHEET:=$sampleSheetFile" if -e $sampleSheetFile;
push @configuration, "ALIGNMENT_CONFIG_TXT:=" . File::Spec->catfile( $validationRunFolder, 'Aligned', 'ValidationConfig.txt');
my $alignmentConfigMk = File::Spec->catfile( $validationRunFolder, 'Aligned', 'ValidationConfig.mk');
push @configuration, "ALIGNMENT_CONFIG_MK:=" . $alignmentConfigMk if -e $alignmentConfigMk;
push @configuration, "BUILD_CFG_ROOT:=" . File::Spec->catfile($validationRunFolder, 'Build');
push @configuration, "PROJECT_SAMPLE_CONFIG_FILES:=" . join (' ', glob(File::Spec->catfile($validationRunFolder, 'Build', 'Project_*', 'Sample_*', 'ValidationConfig.mk')));
push @configuration, "VALIDATE_ACTION:=generate" if $generate;
push @configuration, "VALIDATE_ACTION:=validate" unless $generate;
my $inputMakefile = File::Spec->catfile( '/usr/local/share/bcl2fastq-1.8.4', "makefiles", 'Validation', 'Makefile');
my $outputMakefile = File::Spec->catfile( $outputDir, 'Makefile');
use File::Copy;
copy($inputMakefile, $outputMakefile) or errorExit "ERROR: Failed to copy $inputMakefile to $outputMakefile: $!";

my $outputConfigMk = File::Spec->catfile( $outputDir, 'config.mk');
open (CONFIG_MK, ">$outputConfigMk") or errorExit "ERROR: Unable to open $outputConfigMk for write";
print CONFIG_MK join("\n", @configuration, '') or errorExit("ERROR: Failed to write the configuration: $!");

__END__
