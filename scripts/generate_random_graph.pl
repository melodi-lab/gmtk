#!/usr/nikola/bin/perl -w 

use strict;
use Getopt::Long;
use Graph::Directed;

my $valid_options; 
my $nodes_per_frame; 
my $nmbr_verticies; 
my $max_allowable_edges;
my $min_edges_in_graph;
my $max_edges_in_graph;
my $max_in_degree;
my $nmbr_frames;
my $total_iterations;
my $deterministic;
my $observed;
my $max_rndm_card;
my $max_disc_card;
my $master_file_name;
my $use_discrete;
my $continuous_obs;
my $allow_disconnected;
my @V;
my $i;
my $j;
my $G;

##############################################################################
# Get command line parameters 
##############################################################################

$valid_options = &GetOptions( 
  "nodes:i"         => \$nodes_per_frame, 
  "frames:i"        => \$nmbr_frames, 
  "max_edges:i"     => \$max_allowable_edges, 
  "max_in_degree:i" => \$max_in_degree, 
  "iterations:i"    => \$total_iterations, 
  "deterministic:f" => \$deterministic, 
  "observed:f"      => \$observed, 
  "max_rndm_card:i" => \$max_rndm_card, 
  "max_disc_card:i" => \$max_disc_card, 
  "master_file:s"   => \$master_file_name, 
  "discrete_observed!"  => \$use_discrete,
  "allow_disconnected!" => \$allow_disconnected
  );

##########################################################################
# Get number of nodes and frames 
##########################################################################
if (!defined $nodes_per_frame) {
  $valid_options = 0;
}

if (!$valid_options) 
{
  print "***ERROR, command line arguments are\n";
  print "   -nodes               Number of nodes per frame (required)\n";
  print "   -max_edges           Maximum number of edges\n";
  print "   -max_in_degree       Maximum indgree of any node\n";
  print "   -frames              Number of frames in DBN (1 or 2)\n";
  print "   -iterations          Number of iterations in MCMC\n";
  print "   -deterministic       Probability that a node is deterministic\n";
  print "   -observed            Probability that a node is observed\n";
  print "   -max_rndm_card       Maximum cardinality of random nodes\n";
  print "   -max_disc_card       Maximum cardinality of a discrete node\n";
  print "   -master_file         Name of master file\n";
  print "   -discrete_observed   Use discrete observations (not continuous)\n";
  print "   -allow_disconnected  Allows disconnected graphs\n";
  die "\n";
}

if (!defined $nmbr_frames) {
  $nmbr_frames = 2;
} 

(($nmbr_frames == 1) || ($nmbr_frames == 2)) or die "***ERROR:  Number of frames must be 1 or 2"; 

$nmbr_verticies = $nodes_per_frame*$nmbr_frames; 
($nmbr_verticies > 1) or die "***ERROR:  Must have at least two verticies\n"; 

##########################################################################
# Get number of edges 
##########################################################################
$max_edges_in_graph = ($nmbr_verticies*$nmbr_verticies - $nmbr_verticies)/2; 

if (!defined $max_allowable_edges) {
  if ($nmbr_verticies == 2) {
    $max_allowable_edges = 1; 
  }
  elsif ($nmbr_verticies == 3) {
    $max_allowable_edges = 3; 
  }
  elsif ($max_edges_in_graph/2 < $nmbr_verticies) {
    $max_allowable_edges = $nmbr_verticies+1;
  }
  else { 
    $max_allowable_edges = int($max_edges_in_graph/2);  
  }
}

($max_allowable_edges <= $max_edges_in_graph) or die "***ERROR:  Maximum number of edges with $nmbr_verticies verticies is $max_edges_in_graph\n";  

$min_edges_in_graph = $nmbr_verticies-1;
(($allow_disconnected) || ($max_allowable_edges >= $min_edges_in_graph)) or die
  "***ERROR:  Minimum number of edges with $nmbr_verticies verticies is $min_edges_in_graph\n";  

##########################################################################
# Get other options 
##########################################################################

if (!defined $max_in_degree) {
  $max_in_degree = $nmbr_verticies;
}

if (!defined $total_iterations) {
  $total_iterations = 500;
} 

if (!defined $deterministic) {
  $deterministic = 0.5;
} 

if (!defined $observed) {
  $observed = 0.2;
} 

if (!defined $max_rndm_card) {
  $max_rndm_card = 10;
} 

if (!defined $max_disc_card) {
  $max_disc_card = 5e3;
} 

if (defined $master_file_name)
{
  open MASTER, ">$master_file_name" or die "***ERROR:  Can't open master file \'$master_file_name\' for writing\n";
}


if ((!defined $use_discrete) || (!$use_discrete))
{
  $continuous_obs = 1;
}
else 
{
  $continuous_obs = 0;
}

##############################################################################
# Display comments 
##############################################################################
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
print "% Random DAG with:\n";
if ($nmbr_frames == 1) {
  print "%    1 frame\n"; 
}
else {
  print "%    \n"; 
  print "%    $nmbr_frames frames\n"; 
}

print "%    $nmbr_verticies verticies per frame\n"; 
print "%    Maximum of $max_allowable_edges edges\n";
print "%    Maximum indegree $max_in_degree\n";
print "%    $total_iterations iterations in MCMC chain\n"; 
print "%    $deterministic probability of nodes being deterministic\n"; 
print "%    $observed probability of nodes being observed\n"; 
print "%    Maximum random node cardinality of $max_rndm_card\n"; 
print "%    Maximum discrete node cardinality of $max_disc_card\n"; 
if ($continuous_obs) {
  print "%    Using continuous observations\n"; 
}
else {
  print "%    Using discrete observations\n"; 
}
if ($allow_disconnected) {
  print "%    Graph may be disconnected\n"; 
}
else {
  print "%    Graph must be connected\n"; 
}
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";

##############################################################################
# Initialize graph to be a chain
##############################################################################
my $node_name;

$G = Graph::Directed->new();

for ($i=0; $i<$nmbr_frames; $i++) 
{
  for ($j=0; $j<$nodes_per_frame; $j++)   {
    $node_name = sprintf "%d", ($i*$nodes_per_frame+$j);
    push @V, $node_name; 
    $G->add_vertex( $node_name );
    $G->set_attribute( 'frame', $node_name, $i ); 
    $G->set_attribute( 'name',  $node_name, "v$j" ); 
  }
}

for ($i=0; $i<$nmbr_frames; $i++) 
{
  for ($i=0; $i<$nodes_per_frame-1; $i++) {
    $j = $i + 1;
    $G = $G->add_edge($V[$i], $V[$j]);
  }
}

if ($nmbr_frames == 2) 
{
  $G = $G->add_edge($V[0], $V[$nodes_per_frame]);
}

##############################################################################
# Generate a chain of graphs of length $total_iterations 
##############################################################################
my @changed_edges;
my $moves_accepted;

$moves_accepted = 0;
while($moves_accepted<$total_iterations)
{
  @changed_edges = change_edge($G);

  if ( ((scalar $G->edges) <= $max_allowable_edges)  &&
       (($allow_disconnected) || (is_connected($G))) &&
       (is_DAG($G))                                  &&
       (max_in_degree($G) <= $max_in_degree)         && 
       (($continuous_obs==0) || (scalar get_frame_sinks($G) >= 1)) )
  {
    $moves_accepted++;
  }
  else 
  {
    restore_changed_edges( @changed_edges ); 
  }

}

print  "% $G\n";
printf "%% Edges in graph: %d\n\n", (scalar $G->edges);

##############################################################################
# Randomly assign other graph attributes 
##############################################################################
my $node_index;
my $vertex_b;
my $vertex;
my $rndm_nmbr;
my $crrnt_max_disc_card;
my $self_parent;
my $nmbr_parents; 
my @parents; 
my @parent_indices; 
my @parent_names; 
my $parent_name; 
my $parent; 
my @sorted_vertices; 
my $nmbr_deterministic_CPT; 
my $nmbr_observed; 
my $frame; 
my @sinks;

@sorted_vertices = toposort($G);
@sinks = get_frame_sinks($G);

##############################################################################
# Set the initial cardinalities, and decide if each vertex is random, 
# deterministic, or observed.
##############################################################################
$nmbr_deterministic_CPT = 0;
$nmbr_observed = 0;
for($node_index=0; $node_index<$nodes_per_frame; $node_index++)
{
  $vertex   = $sorted_vertices[$nodes_per_frame*($nmbr_frames-1)+$node_index];
  $vertex_b = $sorted_vertices[$node_index];

  $rndm_nmbr = rand(1);
  $nmbr_parents = scalar $G->predecessors($vertex);

  if ( ($nmbr_parents > 0) && ($rndm_nmbr < $deterministic) ) {
    $G->set_attribute( 'deterministic', $vertex, 1 );    
    $G->set_attribute( 'cardinality',   $vertex, 1 ); 
    $nmbr_deterministic_CPT++;
    if ($nmbr_frames == 2) { 
      $G->set_attribute( 'deterministic', $vertex_b, 1 );    
      $G->set_attribute( 'cardinality',   $vertex_b, 1 ); 
      $nmbr_deterministic_CPT++;
    }
  }
  else {
    $G->set_attribute( 'deterministic', $vertex, 0 );
    $rndm_nmbr = 2 + int(rand($max_rndm_card-2));
    $G->set_attribute( 'cardinality', $vertex, $rndm_nmbr );
    if ($nmbr_frames == 2) { 
      $G->set_attribute( 'deterministic', $vertex_b, 0 );    
      $G->set_attribute( 'cardinality', $vertex_b, $rndm_nmbr );
    }
  }

  ##############################################################################
  # Choose if node is observed
  ##############################################################################
  $rndm_nmbr = rand(1);
  if ($rndm_nmbr < $observed) 
  {
    if (($continuous_obs==0) || (grep /^$vertex$/, @sinks)) {
      set_observed($vertex_b); 
      $nmbr_observed++;
    } 
  }

}

##############################################################################
# Make sure one node is observed 
##############################################################################
if ($nmbr_observed == 0)
{
  if ($continuous_obs)
  {
    ##########################################################################
    # Choose from sink nodes for continuous observation 
    ##########################################################################
    $rndm_nmbr = int(rand(scalar @sinks));
    $vertex = $sinks[$rndm_nmbr];
  }
  else 
  {
    ##########################################################################
    # Choose any node for discrete observation 
    ##########################################################################
    $rndm_nmbr = int(rand($nodes_per_frame));
    $vertex = $V[$rndm_nmbr];
  }

  set_observed($vertex); 
}

##############################################################################
# Possibly increase state space of determistic nodes when parent 
#   cardinalities are large. 
##############################################################################
for($frame=0; $frame<$nmbr_frames; $frame++)
{
  for($node_index=0; $node_index<$nodes_per_frame; $node_index++)
  {
    $vertex   = $sorted_vertices[$frame*$nodes_per_frame+$node_index];
    $vertex_b = $sorted_vertices[$node_index];

    if ( ($G->get_attribute('deterministic', $vertex ) == 1) &&
         !($G->has_attribute('observed', $vertex )) ) 
    { 
      @parents = $G->predecessors($vertex);

      $crrnt_max_disc_card = 1;
      $self_parent = 0;
      foreach $parent (@parents) 
      {
        if ($G->get_attribute('name', $parent) ne 
            $G->get_attribute('name', $vertex))  
        {
          $crrnt_max_disc_card = $crrnt_max_disc_card *
            $G->get_attribute('cardinality', $parent);
        }
        else {
          $self_parent = 1;
        }

        if ($crrnt_max_disc_card > $max_disc_card) {
          $crrnt_max_disc_card = $max_disc_card;
        }
      }

      if ($self_parent)  {
        $crrnt_max_disc_card = $crrnt_max_disc_card*$crrnt_max_disc_card;
      }

      if ($crrnt_max_disc_card > $max_disc_card) {
        $crrnt_max_disc_card = $max_disc_card;
      }

      if ($crrnt_max_disc_card<=3)  {
        $rndm_nmbr = 2;
      }
      else { 
        $rndm_nmbr = 2 + int(rand($crrnt_max_disc_card-3));
      }     
      $G->set_attribute( 'cardinality', $vertex, $rndm_nmbr ); 

      if ($nmbr_frames == 2) { 
        $G->set_attribute( 'cardinality', $vertex_b, $rndm_nmbr );  
      }
    } 
  } 
}

##############################################################################
# Display graph as a structure file 
##############################################################################
my $prnt_index;

print "GRAPHICAL_MODEL Random\n\n";

for($frame=0; $frame<$nmbr_frames; $frame++)
{
  print "frame : $frame\n"; 
  print "{\n"; 

  if ($nmbr_frames == 1)
  {
    print  "   variable : dummy {\n";
    print  "      type: discrete hidden cardinality 2;\n"; 
    print  "      switchingparents: nil;\n"; 
    printf "      conditionalparents: nil\n"; 
    print "         using DeterministicCPT(\"dummy_DtrmnstcCPT\");\n";
    print "   }\n\n";
  }

  for($node_index=0; $node_index<$nodes_per_frame; $node_index++)
  {
    if ($nmbr_frames==1) {
      $vertex = $sorted_vertices[$node_index];
    }
    else {
      $vertex = $sorted_vertices[$frame*$nodes_per_frame+$node_index];
    }

    printf "   variable : %s {\n", $G->get_attribute('name', $vertex);

    if ($G->has_attribute('observed', $vertex)) { 
      if ($continuous_obs) {
        printf "      type: continuous observed 0:%d;\n", 
          ($G->get_attribute('cardinality', $vertex)-1); 
      }
      else {
        printf "      type: discrete observed 0:0 cardinality %d;\n",
          $G->get_attribute('cardinality', $vertex); 
      }
    }
    else {
      printf "      type: discrete hidden cardinality %d;\n", 
        $G->get_attribute('cardinality', $vertex);
    }
 
    print "      switchingparents: nil;\n"; 
    print "      conditionalparents: ";

    @parent_names = (); 
    if (($nmbr_frames == 1) && ($node_index == 0)) {
      @parent_names = ('dummy(0)', @parent_names); 
    }

    @parent_indices = $G->predecessors($vertex);
    for($prnt_index=0; $prnt_index<(scalar @parent_indices); $prnt_index++)
    {
      $parent_name = $G->get_attribute('name', $parent_indices[$prnt_index]);

      if (($nmbr_frames == 1) || 
          ($G->get_attribute('frame', $parent_indices[$prnt_index]) == $frame)) 
      { 
        $parent_name = $parent_name . "(0)";
      }
      else {
        $parent_name = $parent_name . "(-1)";
      }
      push @parent_names, $parent_name;
    }

    if ((scalar @parent_names)==0) {
      print "nil ";
    }
    else { 
      for($prnt_index=0; $prnt_index<(scalar @parent_names); $prnt_index++) 
      {
        printf "%s", $parent_names[$prnt_index];
 
        if ($prnt_index<(scalar @parent_names)-1) {
          print ",";
        }
        print " ";
      }
    } 

    if ($continuous_obs && ($G->has_attribute('observed', $vertex))) {
      print "using mixture\n";
      print "        collection(\"myglobal\")\n";
      print "        mapping(\"${vertex}_DT\");\n"; 
    } 
    elsif ($G->get_attribute('deterministic', $vertex)) {
      print "using DeterministicCPT(\"${vertex}_DeterministicCPT\");\n"; 
    }
    else {
      printf "using DenseCPT(\"%d_DenseCPT\");\n", $vertex; 
    }

    printf "   }\n\n";
  }

  printf "}\n\n";
}

if ($nmbr_frames == 1)
{
  print  "frame : 1\n";   
  print  "{\n";
  print  "   variable : dummy {\n";
  print  "      type: discrete hidden cardinality 2;\n"; 
  print  "      switchingparents: nil;\n"; 
  print  "      conditionalparents: dummy(-1) using DenseCPT(\"dummy_DenseCPT\");\n"; 
  print  "   }\n\n";
  print  "   variable : dummy_2 {\n";
  print  "      type: discrete hidden cardinality 2;\n"; 
  print  "      switchingparents: nil;\n"; 
  print  "      conditionalparents: dummy(0) using DenseCPT(\"dummy_DenseCPT\");\n"; 
  print  "   }\n";
  print  "}\n\n";
}

printf "chunk 1:1\n\n";


##############################################################################
# Write master file 
##############################################################################

if (defined $master_file_name) 
{
  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  print MASTER "% Master file for random DAG\n";
  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";

  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  print MASTER "% Gaussian Mixtures\n";
  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  print MASTER "DPMF_IN_FILE  initialGMParams ascii % dense 1D pmfs\n";
  print MASTER "MEAN_IN_FILE  initialGMParams ascii % means\n";
  print MASTER "COVAR_IN_FILE initialGMParams ascii % variances\n";
  print MASTER "MC_IN_FILE    initialGMParams ascii % gaussian components\n";
  print MASTER "MX_IN_FILE    initialGMParams ascii % mixtures of Gaussians\n";
  print MASTER "\n";

  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  print MASTER "% Decision Trees\n";
  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
  print MASTER "DT_IN_FILE inline\n\n";

  if ($nmbr_frames == 1) {
    $nmbr_deterministic_CPT++; 
  }

  print MASTER "$nmbr_deterministic_CPT  % Number of decision trees\n\n";

  my $DT_count;
  $DT_count = 0;

  if ($nmbr_frames == 1) {
    print  MASTER "$DT_count  % DT number\n";
    $DT_count++; 
    print  MASTER "dummy_DT\n";
    printf MASTER "0 %% Number of parents\n";
    printf MASTER "   -1 0\n\n";
  }

  foreach $vertex (@sorted_vertices)
  {
    if ($G->get_attribute('deterministic', $vertex) == 1) { 

      print MASTER "$DT_count  % DT number\n";
      $DT_count++;
      print MASTER "${vertex}_DT\n";
      printf MASTER "%d %% Number of parents\n", (scalar $G->predecessors($vertex));

      @parents = $G->predecessors($vertex);
      if ((scalar @parents) == 0) {
        printf MASTER "   -1 %d\n", int(rand($G->get_attribute('cardinality', $vertex)));
      }
      else
      {   
        print MASTER "   -1 (mod((1+";
        for($i=0; $i<(scalar @parents); $i++)
        {
          print MASTER "p$i";
          for($j=0; $j<$i; $j++) {
            print MASTER "*c$j";
          }
          if ($i!=(scalar @parents)-1) {
            print MASTER "+";
          }
        }

        printf MASTER "),%d))\n", $G->get_attribute('cardinality', $vertex);
      } 
      print MASTER "\n";  
    } 
  }

  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  print MASTER "% Deterministic CPTs\n";
  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
  print MASTER "DETERMINISTIC_CPT_IN_FILE inline\n\n";

  print MASTER "$nmbr_deterministic_CPT  % Number of decision trees\n\n";

  $DT_count = 0;

  if ($nmbr_frames == 1) {
    print  MASTER "$DT_count  % DT number\n";
    $DT_count++; 
    print  MASTER "dummy_DtrmnstcCPT\n";
    printf MASTER "0 %% Number of parents\n";
    printf MASTER "2\n"; 
    print  MASTER "dummy_DT\n\n";
  }

  foreach $vertex (@sorted_vertices)
  {
    if ($G->get_attribute('deterministic', $vertex) == 1) { 
      print MASTER "$DT_count  % count\n";
      $DT_count++;
      print MASTER "${vertex}_DeterministicCPT\n";
      @parents = $G->predecessors($vertex);
      printf MASTER "%d %% number of parents\n", (scalar @parents); 
      foreach $parent (@parents) 
      {
        printf MASTER "%d ", $G->get_attribute('cardinality', $parent);
      }
      printf MASTER "%d\n", $G->get_attribute('cardinality', $vertex);
      print MASTER "${vertex}_DT\n";
      print MASTER "\n";
    }
  }

  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  print MASTER "% Name Collections\n";
  print MASTER "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  print MASTER "NAME_COLLECTION_IN_FILE inline 1\n";
  print MASTER "0\n"; 
  print MASTER "myglobal % name\n";
  print MASTER "91 % collection length\n"; 
  for ($i=0; $i<91; $i++)
  {
    print MASTER "gm$i\n";
  }

  close MASTER; 
}

##############################################################################
# Checks if graph is connected. 
##############################################################################
sub is_connected
{
  my $G; 
  my $vertex; 
  my $i;
  my $connected; 

  $G = pop @_;
 
  foreach $vertex (@V) { 
    $G->set_attribute('marked', $vertex, 0)
  }

  recurse_tree($G, $V[0]);
  
  $connected = 1;
  for ($i=0; ($i<$nmbr_verticies) && ($connected == 1); $i++) {
    if ($G->get_attribute('marked', $V[$i]) != 1) {
      $connected = 0;
    }
  }

  return($connected);
}

sub recurse_tree
{
  my $G;
  my $vertex;
  my @neighbors;
  my $neighbor;

  $vertex = pop @_;
  $G = pop @_;

  $G->set_attribute('marked', $vertex, 1);
  
  @neighbors = $G->predecessors($vertex);
  foreach $neighbor (@neighbors) { 
    if ($G->get_attribute('marked', $neighbor) == 0) {
      recurse_tree($G, $neighbor);
    }
  }
  
  @neighbors = $G->successors($vertex);
  foreach $neighbor (@neighbors) { 
    if ($G->get_attribute('marked', $neighbor) == 0) {
      recurse_tree($G, $neighbor);
    }
  }
}

##############################################################################
# Checks if graph is acyclic 
##############################################################################
sub is_DAG
{
  my $G;
  my $i;
  my $vertex;
  my $path_found;

  $G = pop @_;

  $path_found = 0;
  for ($i=0; ($i<$nmbr_verticies) && ($path_found == 0); $i++) 
  {
    foreach $vertex (@V) {
      $G->set_attribute('marked', $vertex, 0)
    }

    $path_found = is_directed_path( $G, $V[$i], $V[$i] );
  } 

  return(!$path_found); 
}

sub is_directed_path
{
  my $G;
  my $i;
  my $start;
  my $end;
  my $path_found;
  my @children;
  my $child;

  $end   = pop @_;
  $start = pop @_;
  $G = pop @_;

  $path_found = 0;
  @children = $G->successors($start);
  for ($i=0; ($i<(scalar @children)) && ($path_found == 0); $i++) 
  {
    $child = $children[$i];
    if ($child eq $end) {
      $path_found = 1;
    }
    elsif ($G->get_attribute('marked', $child) == 0) {
      $G->set_attribute('marked', $child, 1);
      $path_found = is_directed_path( $G, $child, $end ); 
    }
  }

  return($path_found);
}


##############################################################################
# Randomly change an edge 
##############################################################################
sub change_edge 
{
  my $rndm_v_1;
  my $rndm_v_2;
  my @changed_edges;

  $rndm_v_1 = int(rand($nodes_per_frame));

  do { 
    $rndm_v_2 = int(rand($nmbr_verticies));
  } while($rndm_v_1 == $rndm_v_2);

  if ($G->has_edge($V[$rndm_v_1], $V[$rndm_v_2])) 
  {
    @changed_edges = ( @changed_edges, delete_edge($G, $rndm_v_1, $rndm_v_2) );  
  }
  elsif ($G->has_edge( $V[$rndm_v_2], $V[$rndm_v_1])) 
  {
    @changed_edges = ( @changed_edges, delete_edge($G, $rndm_v_2, $rndm_v_1) );  
    @changed_edges = ( @changed_edges, add_edge($G, $rndm_v_1, $rndm_v_2) );  
  }
  else 
  {
    @changed_edges = ( @changed_edges, add_edge($G, $rndm_v_1, $rndm_v_2) );  
  }

  return(@changed_edges);
}


##############################################################################
# Wrapper for adding an edge 
##############################################################################
sub add_edge 
{
  my $v1_index;
  my $v2_index;
  my $G;
  my @changed_edges; 
  my $frame; 

  $v2_index = pop @_;
  $v1_index = pop @_;
  $G = pop @_;

  $G->add_edge($V[$v1_index], $V[$v2_index]);  
  @changed_edges = (1, $v1_index, $v2_index);  

  if ($nmbr_frames>1)
  {
    if (($v1_index < $nodes_per_frame) && ($v2_index < $nodes_per_frame)) 
    {
      $G->add_edge($V[$v1_index+$nodes_per_frame], $V[$v2_index+$nodes_per_frame]);  
    } 
  }

  return(@changed_edges);
}


##############################################################################
# Wrapper for removing an edge 
##############################################################################
sub delete_edge 
{
  my $v1_index;
  my $v2_index;
  my $G;
  my @changed_edges; 
  my $frame; 

  $v2_index = pop @_;
  $v1_index = pop @_;
  $G = pop @_;

  $G->delete_edge($V[$v1_index], $V[$v2_index]);  
  @changed_edges = (0, $v1_index, $v2_index);  

  if ($nmbr_frames>1)
  {
    if (($v1_index < $nodes_per_frame) && ($v2_index < $nodes_per_frame)) 
    {
      $G->delete_edge($V[$v1_index+$nodes_per_frame], $V[$v2_index+$nodes_per_frame]);  
    } 
  }

  return(@changed_edges);
}

##############################################################################
# Restore a previously changed edge 
##############################################################################
sub restore_changed_edges 
{
  my $action; 
  my $v1; 
  my $v2;

  do {
    $v2 = pop @_;
    $v1 = pop @_;
    $action = pop @_;

    if ($action == 0) {
      add_edge($G, $v1, $v2);  
    }
    else { 
      delete_edge($G, $v1, $v2);  
    }
  } while( (scalar @_)>0 );

}

##############################################################################
# Sort topologically within a frame 
##############################################################################
sub toposort 
{
  my $G;
  my @sorted_all;
  my @sorted_frame;
  my @sorted_nodes; 
  my $frame; 
  my $node;

  $G = pop @_;

  @sorted_all = $G->toposort;

  foreach $node (@sorted_all) 
  {
    if ($G->get_attribute( 'frame', $node ) == 0) {
      push @sorted_frame, $node;
    }
  }

  for ($frame=0; $frame<$nmbr_frames; $frame++)
  {
    for($i=0; $i<(scalar @sorted_frame); $i++)
    {
      $node = sprintf "%d", $sorted_frame[$i]+($frame*$nodes_per_frame); 
      push @sorted_nodes, $node; 
    } 
  }

  return(@sorted_nodes);
}

##############################################################################
# Max degree 
##############################################################################
sub max_in_degree 
{
  my $G;
  my $vertex;
  my $max_degree;

  $G = pop @_;

  $max_degree = 0;
  foreach $vertex (@V) 
  {
    if ($G->in_degree($vertex) > $max_degree) {
      $max_degree = $G->in_degree($vertex);
    }
  }

  return($max_degree);
}


##############################################################################
# get_frame_sinks 
##############################################################################
sub get_frame_sinks 
{
  my $G; 
  my @all_sinks;
  my @frame_sinks;
  my $i; 
 
  $G = pop @_;
 
  @all_sinks = $G->sink_vertices;
 
  for($i=0; $i<(scalar @all_sinks); $i++)
  { 
    if ($all_sinks[$i] < $nodes_per_frame) { 
      push @frame_sinks, $all_sinks[$i];
    }
  }    

  return(@frame_sinks);
}

##############################################################################
# set_observed
##############################################################################
sub set_observed 
{
  my $vertex;
  my $max_cardinality;

  $vertex = pop @_;

  $G->set_attribute('observed', $vertex, 1);

  if ($continuous_obs) {
    $G->set_attribute('cardinality', $vertex, 42);
    if ( $G->get_attribute('deterministic', $vertex) == 0 ) {
      $nmbr_deterministic_CPT++;
      $G->set_attribute('deterministic', $vertex, 1);
    }
  }
  else {
    $G->set_attribute('cardinality', $vertex, 50);
    if ( $G->get_attribute('deterministic', $vertex) == 1 ) {
      $nmbr_deterministic_CPT--;
      $G->set_attribute('deterministic', $vertex, 0);
    }
  }

  if ($nmbr_frames>1)
  {
    $vertex = $vertex+$nodes_per_frame;
    $G->set_attribute('observed', $vertex, 1);
    if ($continuous_obs) {
      $G->set_attribute('cardinality', $vertex, 42);
      if ( $G->get_attribute('deterministic', $vertex) == 0 ) {
        $nmbr_deterministic_CPT++;
        $G->set_attribute('deterministic', $vertex, 1);
      } 
    }
    else {
      $G->set_attribute('cardinality', $vertex, 50);
      if ( $G->get_attribute('deterministic', $vertex) == 1 ) {
        $nmbr_deterministic_CPT--;
        $G->set_attribute('deterministic', $vertex, 0);
      }
    }
  }
}

