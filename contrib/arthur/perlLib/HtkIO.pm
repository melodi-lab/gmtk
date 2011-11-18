#A PERL module to read and write HTK format observation files
#
#Released under the MIT License. (http://www.opensource.org/licenses/mit-license.php)
#Permision to use, copy, modify, and create derivative works is granted.  NO WARRANTY.
#
#
#Arthur Kantor
#6/24/05 Initial release

package HtkIO;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(&packHtk  unpackHtk);



# Writes out an observation sequence of n k-sized vectors in HTK feature format. 
#
# parameters:
# (my $filename, my $parmKind, my $sampPeriod, my $vectorSize, my @data) = @_;
#
# my $parmKind: see below
# my $sampPeriod: sample Period in 100ns units 
# my $vectorSize: number of elements per sample vector
# my @data: the observation data stored, vector after vector.  Number of elements should be a multiple of $vectorSize
#
#
#  The parmKind can be one of 
# 0 		 WAVEFORM 		 sampled waveform 								
# 1 		 LPC 		     linear prediction filter coefficients 
# 2 		 LPREFC 		 linear prediction reflection coefficients 
# 3 		 LPCEPSTRA 		 LPC cepstral coefficients 
# 5 		 IREFC 		     LPC reflection coef in 16 bit integer format  
# 6 		 MFCC 		     mel-frequency cepstral coefficients 
# 7 		 FBANK 		     log mel-filter bank channel outputs 
# 8 		 MELSPEC 		 linear mel-filter bank channel outputs 
# 9 		 USER 		     user defined sample kind 
# 10 		 DISCRETE 		 vector quantised data 
# 
# modified by ANDing it with any number of the flags (shown in octal)
# N 		 000200 		 absolute energy suppressed 
# _D 		 000400 		 has delta coefficients 
# _A 		 001000 		 has acceleration coefficients
# _Z 		 004000 		 has zero mean static coef. 
# _O 		 020000 		 has 0'th cepstral coef. 
# 
# Every vector element can be a 32-bit float, or a 16-bit integer, depending on the parmKind.
# if parmKind is DISCRETE, then data is first converted to a 16-bit integer,
# otherwise data is converted to a 32-bit float.  
# See section "HTK Format Parameter Files " in the HTK book for how the the flags affect the layout of the sample vector
# All numbers are stored in big-endian order (default in HTK), regardless of the architecture.

sub packHtk{
	(my $parmKind, my $sampPeriod, my $vectorSize, my @data) = @_;
	
	$elCount = @data;
	my $nSamples = $elCount/$vectorSize;
	die "\@data is not a multiple of vectorSize $vectorSize" if ($nSamples != int($nSamples));
	#die "paramKind is DISCRETE requires that vectorSize $vectorSize == 1" if (($parmKind & 077)==10 && $vectorSize != 1);
	
	my $dtype = ($parmKind & 077)==10 ? 2 : 4; #number of bytes
	my $header = pack("NNnn", $nSamples, $sampPeriod, $vectorSize*$dtype, $parmKind);
	my $data;
	if($dtype == 2){ #short int data
		$data = pack("n[$elCount]",@data);
	}
	else{ #float data
		$data = pack("f[$elCount]",@data);
		#unpack as nativeInt and repack as big-endian int, to do the architcture-dependent conversion
		$data= pack("N[$elCount]", unpack("L[$elCount]", $data));
	}
	
	return $header.$data;

}

#the inverse of packHtk
sub unpackHtk{
	my $obs = shift;
	(my $nSamples, my $sampPeriod, my $bytesPerVector, my $parmKind) = unpack("NNnn", $obs);
	my $dtype = ($parmKind & 077)==10 ? 2 : 4; #number of bytes
	my $vectorSize = $bytesPerVector/$dtype;
	my $elCount = $vectorSize*$nSamples;
	my @data;
	if($dtype == 2){ #short int data
		@data = unpack("x[NNnn] n[$elCount]",$obs);
	}
	else{ #float data
		#unpack as big-endian int and repack as nativeInt to do the architcture-dependent conversion
		my $data= pack("L[$elCount]", unpack("x[NNnn] N[$elCount]", $obs));
		@data = unpack("f[$elCount]",$data);
	}
	
	return ($parmKind, $sampPeriod, $vectorSize, @data);
	



}