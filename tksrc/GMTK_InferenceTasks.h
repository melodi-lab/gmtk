
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

