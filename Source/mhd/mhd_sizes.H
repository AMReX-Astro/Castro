#ifndef CASTRO_MHD_SIZES_H
#define CASTRO_MHD_SIZES_H

constexpr int UMAGX = NUM_STATE;
constexpr int UMAGY = NUM_STATE+1;
constexpr int UMAGZ = NUM_STATE+2;

constexpr int NEIGN = 7;

constexpr int IEIGN_RHO = 0;
constexpr int IEIGN_U = 1;
constexpr int IEIGN_V = 2;
constexpr int IEIGN_W = 3;
constexpr int IEIGN_P = 4;

// perpendicular magnetic field components
constexpr int IEIGN_BT = 5;
constexpr int IEIGN_BTT = 6;

#endif
