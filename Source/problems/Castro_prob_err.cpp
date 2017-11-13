#include <Castro.H>

using namespace amrex;

// If you need to take any preparatory steps
// before or after applying the problem tags,
// these hooks can be used.

void
Castro::problem_pre_tagging_hook (TagBoxArray& tags,
                                  int          clearval,
                                  int          tagval,
                                  Real         time)
{
}

void
Castro::problem_post_tagging_hook (TagBoxArray& tags,
                                   int          clearval,
                                   int          tagval,
                                   Real         time)
{
}
