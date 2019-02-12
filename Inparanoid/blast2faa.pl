#!/usr/bin/perl

# Version 1.1

# Supports BLAST up to version 2.2.12. If you are using a higher
# version you can try to add the option -VT to emulate the old
# behaviour.

# BLAST2 and PSI-BLAST output parser
# Extracts the query  names, matching sequence names, scores,
# E-values, lengths from BLAST2 or PSI-BLAST output files.
# Creates 1 file for each match, containing 2 aligned sequences

$S_cutoff = shift;# first argument defines bitscore for cutoff 

while (<>)  { # Loops for each Query sequence
   next if /^\s*$/;
   next if /Sequences not found previously/;
   next if /Sequences used in model/;
   if (/^Query=/) { # New query
      @tmp = split(/\s+/);
      $query = $tmp[1];
      # Get length of query as well
      until(/letters/){
         $_ = <>; # Take next line
	 m/\s+\((\d+)/;
	 $query_length = $1;     
      }
   }


   if (/Results from round/){  # PSI-Blast starts another round
      m/\D+(\d+)/;
      $round = $1;
      $j=0;$sorting_necessary=0; # Reset all values!
   }
     
   if (/^Sequences producin/) { # Has some hits in database!!!
      $_ = <>; # Now on a blank line
      $_ = <> if(/Sequences used in model/);
      $i = 0; # reset match counter
      while (<>) { # Loops until the list of matches
         if (/^\s*$/){# What if blank line? Let's look one more line:
            $_ = <>; # Read next line after the end of match list
	    last if /^\s*$/; 
	    if (/^>/){ # Now on alignment line
      	       $j=0; # alignment count
	       goto ALIGNMENTS; 
	    } 
	    last if /^  Database/; # Stops if no alignments shown
	    last if /^Searching\.\.\./; # Stops if no alignments shown (PSI)
	    if (/Sequences not found previously/){ # This is for PSIBLAST
	       $_ = <>; # Next blank line
	       $sorting_necessary = 1;
	       next;
	    }	
         } 
	 ++$i;
	 
	 chomp;      
	 @tmp=split(/\s+/,$_);    # Take matching name, score and E
	 $match[$i] = $tmp[0];   # match name 
	 $score[$i] = substr ($_,68,5);
	 $E_value[$i] = substr($_,74,5);
	 @tmp = split(//,$E_value[$i]);# Check for correct writing of E-values
	 $E_value[$i] = "1" . $E_value[$i] if $tmp[0] eq "e";
      }# End of match list, continue analyzing alignments
      
      $_ = <>;
      if (/Results from round/){  # PSIblast starts another round
         m/\D+(\d+)/;
	 $round = $1;
	 $j=0;$sorting_necessary=0; # Reset all values!
      }
      
      # Now continue to find length of each matching protein/NA 
      $j=0; 
      while (<>){ # Loops for each query sequence
         last if (/Searching\.\.\./); # New PSI-BLAST round started!
	 ALIGNMENTS:
	 if(/^>/){ # Found lines with alignments
	    ++$j; $HSPcount[$j] = 0;
	    @qstart[$j] = @qend[$j] = ();
	    @mstart[$j] = @mend[$j] = ();
	    @qline[$j] = @mline[$j] = ();
	    chomp;
	    s/\>//;
	    @tmp = split(/\s+/);
	    $match[$j] = $tmp[0];
	    until(/Length =/){ # repeat until target length is found	
	       $_ = <>; # Next line	           
	       @tmp = split(/\s+/); 
	       $match_length[$j] = $tmp[3];
	       next if($tmp[1] ne "Length"); # Some checks:
	       next if($match_length[$j]==NULL or $match_length[$j] =~ /\D/);
	       last;     # Yes, got the length!		  
	    }
	    ##############
	    $_ = <>;$_ = <>; # Score line
	    chomp;
	    s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
	    @tmp = split(/\/|\s+/);
	    $score[$j] = $tmp[2]; # in bits. Overrides the value that was read from initial list
	    
	    $_ = <>; # Now on identity line
	    chomp;
	    s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
	    @tmp = split(/\/|\s+/); 		
	    $HSP_length[$j] = $tmp[3];
	    $id[$j] = $tmp[2];     # Number of identities in match region
	    $pos[$j] = $tmp[6];   # Number similarities in match region
	    $ngaps[$j] = $tmp[10]; # Number of gaps
	    
	    $n = $HSP_count[$j] = 1; #segment number
	    $_ = <>;$_ = <>; # Get query start
	    chomp;
	    @tmp = split(/\s+/);
	    $qstart[$j][$n] = $tmp[1]; # Query start	    
	    $qend[$j][$n] = $tmp[3];   # Query end
	    $qline[$j][$n] = $tmp[2];  # Query line
	    $_ = <>;$_ = <>;
	    chomp;
	    @tmp = split(/\s+/);
	    $mstart[$j][$n] = $tmp[1]; # Match start
	    $mend[$j][$n] = $tmp[3];   # Match end
	    $mline[$j][$n] = $tmp[2];  # Match line
	    $_ = <>;$_ = <>;	
	    for(;;){	   
	       last if /^\s*$/; # empty line after alignment
	       chomp;
	       @tmp = split(/\s+/);
	       $qend[$j][$n] = $tmp[3];
	       $qline[$j][$n] .=  $tmp[2];  # Query line
	       $_ = <>;$_ = <>;
	       chomp;@tmp = split(/\s+/);
	       $mend[$j][$n] = $tmp[3];
	       $mline[$j][$n] .= $tmp[2];  # Match line
	       $_ = <>;$_ = <>; # Get next query start
	    }	
	 } 
	 if((/^ Score =/) and ($n < 5)){ # This happens if we have more than 1 segment reported
 	    ++$n; ++$HSP_count[$j];# HSP numbering 
	    chomp;
	    s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
	    @score_line = split(/\/|\s+/);
	    
	    $_ = <>; # Now on identity line
	    chomp;
	    s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
	    @id_line = split(/\/|\s+/);
	    
	    $_ = <>;$_ = <>; # Get query start 
	    chomp;@tmp = split(/\s+/);
	    $qstart[$j][$n] = $tmp[1];
	    $qend[$j][$n] = $tmp[3];
	    $qline[$j][$n] =  $tmp[2];  # Query line
	    $_ = <>;$_ = <>;
	    chomp;@tmp = split(/\s+/);
	    $mstart[$j][$n] = $tmp[1]; # Match start
	    $mend[$j][$n] = $tmp[3];   # Match end
	    $mline[$j][$n] = $tmp[2];  # Match line
	    $_ = <>;$_ = <>;
	    for(;;){
	       if (/^\s*$/){ # 2 empty lines after alignment
	          #Decide if we want to accept this HSP 
		  # Check for overlap
		  $accepted = 1;
		  for $c(1..($n-1)){
	             # Which HSP is N-terminal?
		     if ($qstart[$j][$n] < $qstart[$j][$c]){
	                $qstart1 = $qstart[$j][$n];
			$qend1   = $qend[$j][$n];
			$mstart1 = $mstart[$j][$n];
			$mend1   = $mend[$j][$n];
			$qstart2 = $qstart[$j][$c];
			$qend2   = $qend[$j][$c];
			$mstart2 = $mstart[$j][$c];
			$mend2   = $mend[$j][$c];
		        $qline1  = $qline[$j][$n]; 
		        $qline2  = $qline[$j][$c];
		        $mline1  = $mline[$j][$n];
			$mline2  = $mline[$j][$c];
		     }
		     else{
		     #if ($qstart[$j][$n]+$mstart[$j][$n] > $qstart[$j][$c]+$mstart[$j][$c]){
	                $qstart2 = $qstart[$j][$n];
			$qend2   = $qend[$j][$n];
			$mstart2 = $mstart[$j][$n];
			$mend2   = $mend[$j][$n];
			$qstart1 = $qstart[$j][$c];
			$qend1   = $qend[$j][$c];
			$mstart1 = $mstart[$j][$c];
			$mend1   = $mend[$j][$c];
	                $qline2  = $qline[$j][$n]; 
		        $qline1  = $qline[$j][$c];
		        $mline2  = $mline[$j][$n];
			$mline1  = $mline[$j][$c];
		     }
#		     print "Working with $match[$j], HSP $n, checking against HSP $c\n";
#		     print "$qstart1-$qend1 $mstart1-$mend1   $qstart2-$qend2 $mstart2-$mend2\n";
		     $qlength1 = $qend1 - $qstart1;
		     $qlength2 = $qend2 - $qstart2;
		     $qlength  = ($qlength1 < $qlength2)?$qlength1:$qlength2;
		     
		     $mlength1 = $mend1 - $mstart1;
		     $mlength2 = $mend2 - $mstart2;
		     $mlength  = ($mlength1 < $mlength2)?$mlength1:$mlength2;
		     
		     $qoverlap = $qstart2 - $qend1;
		     $moverlap = $mstart2 - $mend1;
		     
#		     print "$qoverlap\t$qlength\t$moverlap\t$mlength\n";
		     if (!($qlength * $qlength)){ # Zero length - why?
		        --$n; # Discard this HSP
			--$HSP_count[$j];
			$accepted = 0;
			last;
		     }	
		     # 5% overlap is allowed:
		     if (($qoverlap/$qlength < -0.05) or ($moverlap/$mlength < -0.05)){ 
		        --$n; # Discard this HSP
			--$HSP_count[$j];
#			print "########### Rejected #########\n";
			$accepted = 0;
		        last;
		     }
		  }
		  $score[$j]      += $score_line[2] if($accepted); # add score of this HSP
	          $HSP_length[$j] += $id_line[3]    if($accepted);
		  $id[$j]         += $id_line[2]    if($accepted);
		  $pos[$j]        += $id_line[6]    if($accepted);
		  $ngaps[$j]      += $id_line[10]   if($accepted);
		  last; 
	       }
	       chomp;@tmp = split(/\s+/);
	       $qend[$j][$n] = $tmp[3];
	       $qline[$j][$n] .= $tmp[2];
	       $_ = <>;$_ = <>;
	       chomp;@tmp = split(/\s+/);
	       $mend[$j][$n] = $tmp[3];
	       $mline[$j][$n] .= $tmp[2];
	       $_ = <>;$_ = <>; # Get next query start
	    }
         }  
##########################################################################################
         if(/^effective length of query:/){ # This is query length corrected for edge effects
            chomp;
	    @tmp = split(/\s+/);
	    $M = $tmp[4];
            $M =~ s/,//g;
	 }
	 if(/^effective length of database:/){ # End of this query..., Print results:
            chomp;
	    @tmp = split(/\s+/);
	    $N = $tmp[4]; # Get database length
	    $N =~ s/,//g;
####### Sort HSP on query:
      	    for $a(1..$j){
	       for $b(1..($HSP_count[$a]-1)){
	          while($qstart[$a][$b] > $qstart[$a][$b+1]){
	             $qs = $qstart[$a][$b];
		     $qe = $qend[$a][$b];
		     $ms = $mstart[$a][$b];
		     $me = $mend[$a][$b];
		     $ql = $qline[$a][$b];
		     $ml = $mline[$a][$b];
		     $qstart[$a][$b] = $qstart[$a][$b+1];
		     $qend[$a][$b]   = $qend[$a][$b+1];
		     $mstart[$a][$b] = $mstart[$a][$b+1];
		     $mend[$a][$b]   = $mend[$a][$b+1];
		     $qline[$a][$b]  = $qline[$a][$b+1];
		     $mline[$a][$b]  = $mline[$a][$b+1];
		     $qstart[$a][$b+1] = $qs;
		     $qend[$a][$b+1]   = $qe;
		     $mstart[$a][$b+1] = $ms;
		     $mend[$a][$b+1]   = $me;
		     $qline[$a][$b+1]  = $ql;
		     $mline[$a][$b+1]  = $ml;
		     --$b if ($b > 1);	   
	          }
	       }
	    }	  
	    for $k(1..($j-1)){ # Resort everything by score
	       while($score[$k] < $score[$k+1]){
	          $tempM   = $match[$k];
		  $tempML  = $HSP_length[$k];
		  $tempTL  = $match_length[$k];
		  $tempS   = $score[$k];
		  $tempID  = $id[$k];
		  $tempSIM = $pos[$k];
		  $tempGAP = $ngaps[$k];
		  $tempHC  = $HSP_count[$k];
		  @tempQS  = @qstart[$k];
		  @tempQE  = @qend[$k];
		  @tempMS  = @mstart[$k];
		  @tempME  = @mend[$k];
		  @tempQL  = @qline[$k];
		  @tempML  = @mline[$k];
		  
		  $match[$k] = $match[$k+1];
		  $HSP_length[$k] = $HSP_length[$k+1];
		  $match_length[$k] = $match_length[$k+1];
		  $score[$k] = $score[$k+1];
		  $id[$k] = $id[$k+1];
		  $pos[$k] = $pos[$k+1];
		  $ngaps[$k] = $ngaps[$k+1];
		  $HSP_count[$k] = $HSP_count[$k+1];
		  @qstart[$k] = @qstart[$k+1];
		  @qend[$k] = @qend[$k+1];
		  @mstart[$k] = @mstart[$k+1];
		  @mend[$k] = @mend[$k+1];
		  @qline[$k] = @qline[$k+1];
		  @mline[$k] = @mline[$k+1];
		  
		  $match[$k+1] = $tempM;
		  $HSP_length[$k+1] = $tempML;
		  $match_length[$k+1] = $tempTL;
		  $score[$k+1] = $tempS;
		  $id[$k+1] = $tempID;
		  $pos[$k+1] = $tempSIM;
		  $ngaps[$k+1] = $tempGAP;
		  $HSP_count[$k+1] = $tempHC;
		  @qstart[$k+1] = @tempQS;
		  @qend[$k+1] = @tempQE;
		  @mstart[$k+1] = @tempMS;
		  @mend[$k+1] = @tempME;
		  @qline[$k+1] = @tempQL;
		  @mline[$k+1] = @tempML;
		  
		  --$k if ($k > 1);
	       }
	    }  
	    $count = 0;
            for $k(1..$i){ # Print things out  
  	       next if ($score[$k] < $S_cutoff); # if they make it over cutoff
	       next if ($query eq $match[$k]);   # ignore self
	       ++$count;
	    }
	    last if ($count < 2);
#	    $filename = $query . ".faa";
#	    open F, ">$filename";	    
	    for $k(1..$i){ # Print things out 
	       next if ($score[$k] < $S_cutoff); # if they make it over cutoff
	       next if ($query eq $match[$k]);   # ignore self
	       printf (">%-25.25s", $query);
	       printf ("\t%-25.25s" , $match[$k]);
	       printf ("\t%.1f", $score[$k]);
	       $score[$k] = 25 if ($score[$k] < 25); # To avoid underflow of the E-value
	       $E = $M*$N*(2**-$score[$k]);
	       printf ("\t%3.1g", $E);
	       print  "\t$query_length";
	       	       
	       # Alignments were found, print match_length, HSP_length,
	       # identity %, match % and query-target ratio
	       if ($HSP_length[$k]){
	          print "\t$match_length[$k]";
		  $HSP_length[$k] -= $ngaps[$k] if ($ngaps[$k]);
		  printf ("\t%.0f", $HSP_length[$k]);  # The length of actually aligned amino  acids
                  $match_area = $qend[$k][$HSP_count[$k]] - $qstart[$k][1] + 1;
		  printf ("\t%.0f", $match_area); # The length of matching area on query sequence
		  $id = int(100*($id[$k]/$HSP_length[$k])) if ($HSP_length[$k]);
		  $sim = int(100*($pos[$k]/$HSP_length[$k])) if ($HSP_length[$k]);
		  print "\t$id%";
		  print "\t$sim%";
		  for $hsp (1..$HSP_count[$k]){ # Query HSPs
		     print "\t$qstart[$k][$hsp]-$qend[$k][$hsp]" if ($qstart[$k][$hsp]);
		  }
#		  print "\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
#		  for $hsp (1..$HSP_count[$k]){ # Match HSPs
#                    print "\t$mstart[$k][$hsp]-$mend[$k][$hsp]" if ($mstart[$k][$hsp]);
#	          }
	       }
	       print "\n";
	       $qline = "";
	       $mline = "";
	       for $hsp (1..$HSP_count[$k]){ # Query HSPs
	          $qline .= $qline[$k][$hsp];
		  $mline .= $mline[$k][$hsp];
	       }
	       print "$qline\n";
	       print "$mline\n";
	    }
#	    close F;
	    last; # End of this query
	 }
      }
   }
}