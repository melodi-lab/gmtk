
// InferenceTask subclasses handle inference over a series of sections. The 
// Sections are responsible for doing inference within themselves via whatever
// inference algorithm the Section subclass implements. The InferenceTask 
// classes just handle the messages between Sections via instances of InterfaceSeparator.

// Specific sequence inference algorithms multiply inherit the InferenceTask
// APIs they support. Not all sequence inference algorithms support every 
// inference task, e.g., Island can't do online filtering/smoothing.


// consider ditching the InferenceTask class, making all the *Task classes
// pure abstract bases for the section inference algorithms, which could
// then take the section_inference_alg as a template parameter...

class InferenceTask {
 public:
  InferenceTask(SectionInferenceAlgorithm &section_inference_alg) : 
    section_inference_alg(section_inference_alg) {}

  virtual ~InferenceTask() {}

 protected:

  SectionInferenceAlgorithm &section_inference_alg;
};


class ProbEvidenceTask : public virtual InferenceTask {
 public:
  ProbEvidenceTask(SectionInferenceAlgorithm &section_inference_alg) : 
    InferenceTask(section_inference_alg) {}

  virtual ~ProbEvidenceTask() {}

  // compute log P(X_{0:T-1}), forward pass only, O(1) memory
  virtual logpr probEvidence() {
    InterfaceSeparator msg = section_inference_alg.computeForwardInterfaceSeparator(0);
    unsigned t;
    for (t=1; t < T-1; ++t) {
      section_inference_alg.receiveForwardInterfaceSeparator(t, msg);
      msg = section_inference_alg.computeForwardInterfaceSeparator(t);
    }
    section_inference_alg.receiveForwardInterfaceSeparator(t, msg);
    section_inference_alg.computeForwardInterfaceSeparator(t); // E's gather into root
    printf("ProbE = %f\n", section_inference_alg.probEvidence(t));
  }
};


class ViterbiTask : public virtual InferenceTask {
 public:
  ViterbiTask(SectionInferenceAlgorithm &section_inference_alg) : 
    InferenceTask(section_inference_alg) {}

  virtual ~ViterbiTask() {}

  // compute argmax P(Q_{0:T-1} | X_{0:T-1})
  //   should this be a k element vector of logprs ?
  virtual vector<logpr> viterbi(unsigned kBest=1) = 0; // returns log prob(s) of k-best MPEs
};


class ForwardBackwardTask : public virtual InferenceTask {
 public:
  ForwardBackwardTask(SectionInferenceAlgorithm &section_inference_alg) :
    InferenceTask(section_inference_alg) {}

  virtual ~ForwardBackwardTask() {}

  // compute P(Q_{0:T-1} | X_{0:T-1})
  virtual logpr forwardBackward() = 0;
};


// Note - this uses StreamSource rather than FileSource
class OnlineSmoothingTask : public virtual InferenceTask {
 public:
  OnlineSmoothingTask(SectionInferenceAlgorithm &section_inference_alg) : 
    InferenceTask(section_inference_alg) {}

  virtual ~OnlineSmoothingTask() {}

  // compute argmax P(Q_t | X_{0:t+\tau})
  virtual logpr onlineSmoothing(unsigned futureSections=0, unsigned kBest=1);
};


// Now the GMTK programs can specify which tasks they require, e.g.:
// gmtkViterbi requires inference algorithms with the ViterbiTask interface,
// gmtkOnline requires inference algorithms with the OnlineSmoothingTask 
// interface, etc.

class IslandInference : public ProbEvidenceTask, 
                        public ViterbiTask, 
                        public ForwardBackwardTask
{...};


class ArchipelagosInference : public ProbEvidenceTask, 
                              public ViterbiTask, 
                              public ForwardBackwardTask
{...};


class OnlineInference :  public ProbEvidenceTask, 
                         public OnlineSmoothingTask
{...};


// cuurent "standard" inference
class ???Inference : public ProbEvidenceTask,
                     public ViterbiTask,
                     public ForwardBackwardTask
{...};


// O(Tn) memory
class OnlyKeepSeperatorInference : public ProbEvidenceTask, 
                                   public ViterbiTask, 
                                   public ForwardBackwardTask
{...};


// Assumes no edges between sections: no receive{For,Back}wardInterfaceSeparator() calls.
// The other InferenceTask implelentations should be able to handle static models
// (they'll just pass around empty InterfaceSeparators), so this isn't strictly necessary.
class StaticInference : public ProbEvidenceTask,
                        public ViterbiTask,
                        public forwardBackwardTask
{...};


// Specific section inference algorithm implementations:

//   Pedagogical 
//   Sparse join

class SectionInferenceAlgorithm {
 public:

  virtual ~SectionInferenceAlgorithm() {}

  // Prepare to do inference on segment of length T. Note that this
  // is the analog of the current JunctionTree::unroll() which sets up the
  // PartitionStructureArray for at most  P' C' C' E' (4 sections, so
  // O(1) memory), and the PartitionTableArray for between 0 and T
  // sections depending on ZeroTable, ShortTable, or LongTable. 
  // JT::unroll() also allocates O(T) memory (optionally memory mapped)
  // to store the the Viterbi values if they aren't being written to
  // a file.

  // Since we now have the ability to write Viterbi values to files, do
  // we want to drop the ability to store them in memory? This would
  // break things, as a Viterbi file would then be a required argument.
  // It would also be slower for applications that do fit in memory
  // because of the higher overhead of file I/O operations.

  // Maybe rename this to prepareForSegment() or something to avoid confusion
  // with O(T) space graph unrolling?
  virtual unsigned unroll(unsigned T) = 0; // returns number of usable frames


  // All message actions are named from the perspective of C_t.

  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  virtual InterfaceSeparator computeForwardInterfaceSeparator(unsigned t) = 0; 

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  virtual void receiveForwardInterfaceSeparator(unsigned t, InterfaceSeparator const &msg) = 0;


  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  virtual InterfaceSeparator computeBackwardsInterfaceSeparator(unsigned t) = 0;

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  virtual void receiveBackwardInterfaceSeparator(unsigned t, InterfaceSeparator const &msg) = 0;


  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass
  virtual logpr probEvidence(unsigned t) = 0;


  // Section subclasses can manage their own message ordering w/in a section.
  // Read/write the section's inference plan (JT, msg orders, etc)
  void readInferencePlan(iDataStreamFile *f);
  void writeInferencePlan(ioDataStreamFile *f);
};



class Section {

  /*

   - A section consists of a bunch of cliques and a bunch of separators
    -- currently a section has to be a tree over cliques creating a sub-tree of a junction tree.
   - The cliques are typically maxcliques, but they need not be.
       - in the current triangulation strategies they are typically maxcliques
       - sometimes we have non-max cliques to do projections (e.g., VE stuff, 
         deep model soft training, etc.)
       - in future inference models, we might have a section consist of a non-triangulated
         graph of cliques which are not maxcliques, where the idea is to do LBP over
         these cliques. I.e., we can think of this a a graph (not nec. a tree) of cliques.
      Thus, a section is really a graph of cliques in its most general form.
    - Curently, messages between cliques go indirect via separators, i.e.,
          S = V \cap U where V and U are clique and S is the index set for the separator between
         clique consisting of variables with indices V and U.
      Messages are hugin style.
    - in future, we want to allow messages to be between cliques only, without needing necesarily
       to have separators between cliques. I.e., if we want to do LBP on just the original random
       variables, then each original model factor would be a clique and we'd send messages between
       these factors (c.f. this is essentially equivalant to inference on factor graphs).
      Thus, sections will (in future) not necessarily need to have separators between cliques
      in a junction graph 

    - the connection between sections, however, will always be via an "interfaceSeparator"
      which consists of variables constituting the intersection between two neighboring sections.
    - an interfaceSeparator might be either monolithic (one big complete graph) or factored
       (i.e., an interfaceSeparator itself might be a graph of cliques, or a tree of cliques
       or a junction tree of cliques). - GMTK should support all possibilities (but the
       construction of the topology of the interfaceSeparator might be done externally, or via
       some future program, such as what Shenjie's research might be about. Default 
       instances should be 'completed' and 'all factored' the latter case being a mean-field
       approach (and is similar to the primary instance of the  BK algorithm).


Key abstract routines that we need:

1) A routine that does the necessary computation to prepare the right-side interface separator of section $t$.
   After this is done, it should be possible, once the section $t+1$'s core data structures are ready, to
   receive a message from section $t$. 

   prepareRightInterfaceSeparator()

   On the right most section, this routine might need to be special cased (since there is no need to prepare the same form of right interface if there is no section on the right, but on the other hand we still need to do everything internally so that the section, which might span many time frames, has propagated all internal information to its extreme right).

2) A routine, taken from the perspective of section $t+1$ to receive the information in section $t$'s right interface
   separator.
   A routine named from the perspective of section $t+1$.
   receiveFromLeftNeighborsRightInterfaceSeparator()

3) When section $t+1$'s left interface has been computed during a backwards inference pass, we need a routine
   that takes the info in section $t+1$'s left interface and updates section $t$'s right interface appropriately.

   A routine named from the perspective of section $t$.
   receiveFromRightNeighborsLeftInterfaceSeparator()
      
4) Once section $t$ has received info from neighbor $t+1$ on its right, it needs to do more inference internally
   in such a way that its left interface has received that information. This happens during a backwards pass.

   prepareLeftInterfaceSeparator()

   On the left most section, this routine might need to be special cased (since there is no need to prepare the same form of left interface if there is no section on the left, but on the other hand we still need to do everything internally so that the section, which might span many time frames, has propagated all internal information to its extreme left).


So forward message passing would be of the form:

prepareRightInterfaceSeparator(t);
receiveFromLeftNeighborsRightInterfaceSeparator(t,t+1)
prepareRightInterfaceSeparator(t+1);
receiveFromLeftNeighborsRightInterfaceSeparator(t+1,t+2)
...


So backwards message passing would be of the form:

prepareLeftInterfaceSeparator()
receiveFromRightNeighborsLeftInterfaceSeparator()


  */  



};



// Initially this might be just a single interface clique as in the current
// inference implementation. Eventually, it could be a junction forest or general graph.
// Section subclasses must know how to "send" and "recieve" InterfaceSeparators in each
// direction (CE/DE).

class InterfaceSeparator {
 public:
  
  // read/write the factored interface spec
  virtual static void readInterface(iDataStreamFile *f);
  virtual static void writeInterface(ioDataStreamFile *f);
  
  // methods to make structure and values available to Section inference alg
};

