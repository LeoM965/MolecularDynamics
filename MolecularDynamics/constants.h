#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants {
    const double ANGSTROM_TO_METERS = 1e-10;
    const double KCAL_PER_MOL_TO_JOULES = 4184.0;
    const double AVOGADRO = 6.022e23;
    const double ATOMIC_MASS_UNIT = 1.661e-27;
    const double BOLTZMANN = 1.381e-23;
    const double PI = 3.141592653589793;

    const double BOND_FORCE_CONST = 350.0 * KCAL_PER_MOL_TO_JOULES / (ANGSTROM_TO_METERS * ANGSTROM_TO_METERS);
    const double ANGLE_FORCE_CONST = 40.0 * KCAL_PER_MOL_TO_JOULES;
    const double MAX_FORCE = 1.0e-9;

    const double TETRAHEDRAL_BOND_LENGTH = 1.09;
    const double TETRAHEDRAL_S = TETRAHEDRAL_BOND_LENGTH / sqrt(3.0);
    const double TETRAHEDRAL_ANGLE = 109.5;

    const double NH3_BOND_LENGTH = 1.01;
    const double NH3_ANGLE = 107.0;
    const double NH3_ANGLE_RAD = NH3_ANGLE * PI / 180.0;
    const double NH3_COS_ANGLE = cos(NH3_ANGLE_RAD);
    const double NH3_HEIGHT_FACTOR = sqrt((1.0 + 2.0 * NH3_COS_ANGLE) / 3.0);
    const double NH3_RADIUS_FACTOR = sqrt((2.0 * (1.0 - NH3_COS_ANGLE)) / 3.0);

    const double WATER_BOND_LENGTH = 0.96;
    const double WATER_ANGLE = 104.5;
    const double WATER_HALF_ANGLE_RAD = WATER_ANGLE * PI / 360.0;
    const double WATER_DX = WATER_BOND_LENGTH * sin(WATER_HALF_ANGLE_RAD);
    const double WATER_DY = WATER_BOND_LENGTH * cos(WATER_HALF_ANGLE_RAD);

    const double CO2_BOND_LENGTH = 1.16;
    const double ACETYLENE_BOND_LENGTHS[3] = { 1.06, 1.20, 1.06 };
    const double HCL_BOND_LENGTH = 1.27;

    const double DEG_TO_RAD = PI / 180.0;

    const double LJ_EPSILON[6] = { 0.07, 0.03, 0.15, 0.17, 0.05, 0.23 }; // C, H, O, N, Be, Cl
    const double LJ_SIGMA[6] = { 3.47, 2.65, 3.07, 3.25, 2.80, 3.60 };   // C, H, O, N, Be, Cl
    const double LJ_CUTOFF = 10.0; // Distanta de cutoff in angstromi
}

#endif