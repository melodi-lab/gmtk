
// SequenceInference handles inference over a series of sections. The sections
// are responsible for doing inference within themselves via some
// SectionInferenceAlgorithm. The SequenceInference just handles the messages
// between sections via instances of Separator.

// Subclasses of SequenceInference:
//   Island
//   Archipelagos
//   Pedagogical 
//   Sparse join
//   O(Tn) memory
//   Online smoothing
//   Static

// Not all subclasses implement every inference task, e.g., Island can't
// do online filtering/smoothing.

class SequenceInference {
 public:

  virtual ~SequenceInference() {}

  // compute log P(X_{0:T-1})
  virtual logpr probEvidence() {
    error("not implemented"); return 0.0;
  }

  // compute argmax P(Q_{0:T-1} | X_{0:T-1})
  virtual logpr viterbi(unsigned nBest=1) {
    error("not implemented"); return 0.0;
  }

  // compute P(Q_{0:T-1} | X_{0:T-1})
  virtual logpr forwardBackward() {
    error("not implemented"); return 0.0;
  }

  // compute argmax P(Q_t | X_{0:t+\tau})
  virtual logpr onlineSmoothing(unsigned futureSections=0) {
    error("not implemented"); return 0.0;
  }

};


// Subclasses can manage their own message ordering w/in a section

class SectionInferenceAlgorithm {
 public:

  virtual ~SectionInferenceAlgorithm() {}

  // compute forward message for C'_t -> C'_{t+1}
  virtual Separator computeForwardSeparator(Section C_t) = 0;  // aka gather into root

  // recieve forward message for C'_{t-1} -> C'_t
  virtual void receiveForwardSeparator(Section C_t, Separator msg) = 0;  // aka sendForwardAcross

  // compute backward message for C'_{t-1} <- C'_t
  virtual Separator computeBackwardsSeparator(Section C_t) = 0;  // aka scatter out of root

  // recieve backward message for C'_t <- C'_{t+1}
  virtual void receiveBackwardSeparator(section C_t, Separator msg) = 0;  // aka sendBackwardAcross

};


// Initially this might be just a single interface clique as in the current inference
// implementation. Eventually, it could be a general junction forest. SectionInferenceAlgorithm
// subclasses must know how to "send" and "recieve" Separators in each direction (CE/DE).
// Probably shouldn't be subclassed unless there can be a clean interface between 
// Separator and SectionInferenceAlgorithm.

class Separator {
  // read/write these like tri files / jt files
  
  // methods to make structure and values available to SectionInferenceAlgorithm
};

