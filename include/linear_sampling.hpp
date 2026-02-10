#ifndef LINEAR_SAMPLING_HPP
#define LINEAR_SAMPLING_HPP

#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "fem.hpp"

class LinearSampling {
public:

    // -------------------------------------------------------------------------
    // Calcul de la composante du champs diffracté sur un bord 
    // -------------------------------------------------------------------------
    static void compute_u_s(const MeshP2& mesh, int n_mode, int boundary_tag, double k0, double h, vector<complexe>& u_s, vector<complexe>& u_n) {
        

    }
    // -------------------------------------------------------------------------
    // Projection du champ diffracté n sur le mode m
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Calcul de F
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Calcul de F++ (back scattering)
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Calcul de G 
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Generation d'image matlab de log(1/||h||)
    // -------------------------------------------------------------------------
}
#endif // LINEAR_SAMPLING_HPP
