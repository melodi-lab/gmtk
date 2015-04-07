
// InferenceTask subclasses handle inference over a series of sections. The 
// Sections are responsible for doing inference within themselves via whatever
// inference algorithm the Section subclass implements. The InferenceTask 
// classes just handle the messages between Sections via instances of Separator.

// Specific sequence inference algorithms multiply inherit the InferenceTask
// APIs they support. Not all sequence inference algorithms support every 
// inference task, e.g., Island can't do online filtering/smoothing.


class InferenceTask {
 public:
  InferenceTask(SectionInferenceAlgorithm &inference_alg) : 
    inference_alg(inference_alg) {}

  virtual ~InferenceTask() {}

 protected:

  SectionInferenceAlgorithm &inference_alg;
};


class ProbEvidenceTask : public virtual InferenceTask {
 public:
  ProbEvidenceTask(SectionInferenceAlgorithm &inference_alg) : 
    InferenceTask(inference_alg) {}

  virtual ~ProbEvidenceTask() {}

  // compute log P(X_{0:T-1}), forward pass only, O(1) memory
  virtual logpr probEvidence() {
    auto_ptr<Section> P(inference_alg.getSection(0)); // Section factory method

    // Hmm... not sure about the auto_ptr. If the Sections can be
    // backed by arrays, the Section ownership should stay with
    // the SectionInferenceAlgorithm. Either the Section dtor needs to
    // do the right thing, or perhaps add Section::releaseSection()
    // method that clients have to remember to call

    Separator msg = P->computeForwardSeparator();
    for (int t=1; t < T-1; ++t) {
      auto_ptr<Section> C_t(inference_alg.getSection(t));
      C_t->receiveForwardSeparator(msg);
      msg = C_t->computeForwardSeparator();
    }
    auto_ptr<Section> E(inference_alg.getSection(T-1));
    E->receiveForwardSeparator(msg);
    (void) E->computeForwardSeparator(); // E's gather into root
    printf("ProbE = %f\n", E->probEvidence());
  }
};


class ViterbiTask : public virtual InferenceTask {
 public:
  ViterbiTask(SectionInferenceAlgorithm &inference_alg) : 
    InferenceTask(inference_alg) {}

  virtual ~ViterbiTask() {}

  // compute argmax P(Q_{0:T-1} | X_{0:T-1})
  virtual logpr viterbi(unsigned nBest=1) = 0;
};


class ForwardBackwardTask : public virtual InferenceTask {
 public:
  ForwardBackwardTask(SectionInferenceAlgorithm &inference_alg) :
    InferenceTask(inference_alg) {}

  virtual ~ForwardBackwardTask() {}

  // compute P(Q_{0:T-1} | X_{0:T-1})
  virtual logpr forwardBackward() = 0;
};


class OnlineSmoothingTask : public virtual InferenceTask {
 public:
  OnlineSmoothingTask(SectionInferenceAlgorithm &inference_alg) : 
    InferenceTask(inference_alg) {}

  virtual ~OnlineSmoothingTask() {}

  // compute argmax P(Q_t | X_{0:t+\tau})
  virtual logpr onlineSmoothing(unsigned futureSections=0);
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

class OnlineInference : public OnlineSmoothingTask {...}

// O(Tn) memory
class OnlyKeepSeperatorInference : public ProbEvidenceTask, 
                                   public ViterbiTask, 
                                   public ForwardBackwardTask
{...};

// assumes no edges between sections: no receive{For,Back}wardSeparator() calls
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

  // Prepare to do inference on segment of length T
  virtual void unroll(unsigned T) = 0;

  // factory method returning a C_t Section of the appropriate subclass
  virtual Section *getSection(unsigned section_number) = 0;
};


// Section subclasses can manage their own message ordering w/in a section

class Section {

  // read/write the section's inference plan (JT, msg orders, etc)
  void readInferencePlan(iDataStreamFile *f);
  void writeInferencePlan(ioDataStreamFile *f);

  // in the following, this is C_t

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




  // compute forward message for C'_t -> C'_{t+1} (aka gather into root)
  virtual Separator computeForwardSeparator() = 0; 

  // recieve forward message for C'_{t-1} -> C'_t (sendForwardsCrossPartitions)
  virtual void receiveForwardSeparator(Separator msg) = 0;

  // compute backward message for C'_{t-1} <- C'_t (aka scatter out of root)
  virtual Separator computeBackwardsSeparator() = 0;

  // recieve backward message for C'_t <- C'_{t+1} (sendBackwardCrossPartitions)
  virtual void receiveBackwardSeparator(Separator msg) = 0;

  // return P(Q_t | X_{?}), where ? depends on the messages C_t has seen so far:
  //        P(Q_t | X_{0:t}) in the forward pass
  //        P(Q_t | X_{0:T-1}) in a (full) backward pass
  //        P(Q_t | X_{0:t+\tau}) in a smoothing backward pass
  virtual logpr probEvidence() = 0;
};



// Initially this might be just a single interface clique as in the current
// inference implementation. Eventually, it could be a general junction forest.
// Section subclasses must know how to "send" and "recieve" Separators in each
// direction (CE/DE).

class Separator {
 public:
  
  // read/write the factored interface spec
  void readInterface(iDataStreamFile *f);
  void writeInterface(ioDataStreamFile *f);
  
  // methods to make structure and values available to Section inference alg
};

