
#include "AMReX_LevelBld.H"
#include "Castro.H"

class CastroBld
    :
    public LevelBld
{
    virtual void variableSetUp ();
    virtual void variableCleanUp ();
    virtual AmrLevel *operator() ();
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time);
};

CastroBld Castro_bld;

LevelBld*
getLevelBld ()
{
    return &Castro_bld;
}

void
CastroBld::variableSetUp ()
{
    Castro::variableSetUp();
}

void
CastroBld::variableCleanUp ()
{
    Castro::variableCleanUp();
}

AmrLevel*
CastroBld::operator() ()
{
    return new Castro;
}

AmrLevel*
CastroBld::operator() (Amr&            papa,
                       int             lev,
                       const Geometry& level_geom,
                       const BoxArray& ba,
                       Real            time)
{
    return new Castro(papa, lev, level_geom, ba, time);
}
