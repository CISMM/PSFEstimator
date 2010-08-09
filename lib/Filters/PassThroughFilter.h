#ifndef _PASS_THROUGH_FILTER_H_
#define _PASS_THROUGH_FILTER_H_


/**
 * This filter does nothing but serve as a chunk of code that is placed
 * in the filtering library. Add your own image filter pipelines or ITK
 * extensions here.
 */
class PassThroughFilter {

public:
  PassThroughFilter();
  virtual ~PassThroughFilter();

};

#endif // _PASS_THROUGH_FILTER_H_