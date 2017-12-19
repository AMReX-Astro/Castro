  int cnt = 0;
  Density = cnt;
  cnt += 1;

  Xmom = cnt;
  cnt += 1;

  Ymom = cnt;
  cnt += 1;

  Zmom = cnt;
  cnt += 1;

#ifdef HYBRID_MOMENTUM
  Rmom = cnt;
  cnt += 1;
#endif

#ifdef HYBRID_MOMENTUM
  Lmom = cnt;
  cnt += 1;
#endif

#ifdef HYBRID_MOMENTUM
  Pmom = cnt;
  cnt += 1;
#endif

  Eden = cnt;
  cnt += 1;

  Eint = cnt;
  cnt += 1;

  Temp = cnt;
  cnt += 1;

  if (NumAdv > 0) {
    FirstAdv = cnt;
    cnt += NumAdv;
  }

  if (NumSpec > 0) {
    FirstSpec = cnt;
    cnt += NumSpec;
  }

  if (NumAux > 0) {
    FirstAux = cnt;
    cnt += NumAux;
  }

#ifdef SHOCK_VAR
  Shock = cnt;
  cnt += 1;
#endif

