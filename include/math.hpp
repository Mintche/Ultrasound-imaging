#ifndef MATH_HPP
#define MATH_HPP
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <cmath>
#include <complex>
using namespace std;
typedef complex<double> complexe;
//---------------------------------------------------------------------------
//  Opérations sur vector<T>
//---------------------------------------------------------------------------
template<typename T> vector<T> operator+(const vector<T>& u, const vector<T>& v)
{
    vector<T> w(u);
    auto itv = v.begin();
    for(auto itw=w.begin(); itw!=w.end(); ++itw, ++itv) *itw+=*itv;
    return w;
}
template<typename T> vector<T> operator-(const vector<T>& u, const vector<T>& v)
{
    vector<T> w(u);
    auto itv = v.begin();
    for(auto itw=w.begin(); itw!=w.end(); ++itw, ++itv) *itw-=*itv;
    return w;
}
template<typename T> vector<T> operator*(const vector<T>& u, const T& s)
{
    vector<T> w(u);
    for(auto& wi : w) wi*=s;
    return w;
}
template<typename T> vector<T> operator*(const T& s, const vector<T>& u)
{
    vector<T> w(u);
    for(auto& wi : w) wi*=s;
    return w;
}
template<typename T> vector<T> operator/(const vector<T>& u, const T& s)
{
    vector<T> w(u);
    for(auto& wi : w) wi/=s;
    return w;
}

template<typename T> T operator|(const vector<T>& u, const vector<T>& v)
{
    T s= T(0);
    auto itv = v.begin();
    for(auto itu=u.begin(); itu!=u.end(); ++itu, ++itv) s+=*itu * *itv;
    return s;
}

template<>
complexe operator|(const vector<complexe>& u, const vector<complexe>& v)
{
    complexe s= complexe(0);
    auto itv = v.begin();
    for(auto itu=u.begin(); itu!=u.end(); ++itu, ++itv) s+=conj(*itu) * (*itv);
    return s;
}

template<typename T> T norm(const vector<T>&u)
{
    return sqrt(abs(u|u));
}
template<typename T> ostream& operator<<(ostream& os,const vector<T>& v)
{
  os<<"(";
  auto itv=v.begin();
  for(;itv!=v.end()-1;++itv) os<<(*itv)<<",";
  os<<(*itv)<<")";
  return os;
}

//---------------------------------------------------------------------------
//     classe FullMatrix
//---------------------------------------------------------------------------

template <typename T>
class FullMatrix {

protected:

    vector<T> coefs;
    int n_rows;
    int n_cols;
    
public:

    FullMatrix(int n, int m) : n_rows(n), n_cols(m), coefs(n * m, T(0)) {}

    // Accès

    T& operator()(int i, int j) { return coefs[i * n_cols + j]; }
    const T& operator()(int i, int j) const { return coefs[i * n_cols + j]; }

    // Opérateurs
    
    vector<T> operator*(const vector<T>& x) const {
        // Ajout d'une vérification de la taille du vecteur
        if (x.size() != static_cast<size_t>(n_cols)){
            throw std::invalid_argument("Erreur: Le nombre de colonnes de la matrice doit être égal à la taille du vecteur.");
        }
        vector<T> res(this->n_rows, T(0));
        for(int i = 0; i < n_rows; ++i) {
            for(int j = 0; j < n_cols; ++j) {
                res[i] += (*this)(i, j) * x[j];
            }
        }
        return res;
    }

    void operator+=(const FullMatrix<T>& M){
        // vérification de la taille
        if (n_rows != M.n_rows || n_cols != M.n_cols){
            throw invalid_argument("Erreur: Addition de matrices de tailles différentes.");
        }
        for (size_t i = 0; i < coefs.size(); i++){
            coefs[i] += M.coefs[i];
        }
    }

    void operator-=(const FullMatrix<T>& M){
        // vérification de la taille
        if (n_rows != M.n_rows || n_cols != M.n_cols){
            throw invalid_argument("Erreur: Soustraction de matrices de tailles différentes.");
        }
        for (size_t i = 0; i < coefs.size(); i++){
            coefs[i] -= M.coefs[i];
        }
    }

    // Résolution via pivot de Gauss

    void solve(vector<T>& x, const vector<T>& b){
        // La résolution par pivot de Gauss ne s'applique qu'aux matrices carrées
        if (n_rows != n_cols){
            throw logic_error("Erreur: solve() ne peut être appelée que pour des matrices carrées.");
        }
        if (b.size() != static_cast<size_t>(n_rows)) {
            throw invalid_argument("Erreur: La taille du vecteur b doit correspondre au nombre de lignes de la matrice.");
        }

        int n = n_rows;
        FullMatrix<T> A = *this; 
        x = b; 

        for (int k = 0; k < n; ++k) {
            
            // --- AJOUT PIVOT ---
            int pivot = k;
            auto max_val = std::abs(A(k,k));
            
            // Recherche du meilleur pivot dans la colonne k
            for(int i = k + 1; i < n; ++i) {
                if(std::abs(A(i,k)) > max_val) {
                    max_val = std::abs(A(i,k));
                    pivot = i;
                }
            }
            
            // Si le pivot est nul, la matrice est singulière
            if (max_val < 1e-12) { 
                // Lancer une exception pour une matrice singulière
                throw runtime_error("Erreur: Matrice singulière ou presque singulière.");
            }

            // Échange des lignes dans A
            if (pivot != k) {
                for (int j = k; j < n; ++j) { // On peut commencer à j=k (avant c'est 0)
                    swap(A(k,j), A(pivot,j));
                }
                // Échange dans le vecteur second membre x
                swap(x[k], x[pivot]);
            }

            // Suite normale de l'élimination
            for (int i = k + 1; i < n; ++i) {
                T factor = A(i, k) / A(k, k);
                for (int j = k; j < n; ++j) { // Optimisation: commencer à j=k
                    A(i, j) -= factor * A(k, j);
                }
                x[i] -= factor * x[k];
            }
        }

        // Substitution arrière
        for (int i = n - 1; i >= 0; --i) {
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A(i, j) * x[j];
            }
            x[i] /= A(i, i);
        }
    }

    void print(ostream& os) const {
        os << "(";
        for (int i = 0; i < n_rows; i++){
            for (int j = 0; j < n_cols; j++){
                os << (*this)(i,j) << " ";
            }
            os <<endl;
        }
        os << ")";
    }
};

template<typename T>
FullMatrix<T> operator+(const FullMatrix<T>& A, const FullMatrix<T>& B){
    FullMatrix<T> R = A;
    R += B;
    return R;
}

template<typename T>
FullMatrix<T> operator-(const FullMatrix<T>& A, const FullMatrix<T>& B){
    FullMatrix<T> R = A;
    R -= B;
    return R;
}

//---------------------------------------------------------------------------
//     Classe ProfileMatrix
//---------------------------------------------------------------------------

template <typename T>
class ProfileMatrix {
private:

    vector<T> coefs;
    vector<size_t> p;
    vector<size_t> q;
    int n;

public:

    ProfileMatrix(const vector<size_t>& p_in) : n(p_in.size()), p(p_in) {
        q.resize(n);
        size_t total_size = 0;
        for (size_t i = 0; i < n; i++) {
            q[i] = i - p[i] + 1;
            total_size += q[i];
        }
        coefs.resize(total_size, T(0));
    }

    // Accès

    T operator()(int i, int j) const {
        if (i >= j) {
            if (j >= p[i]){
                size_t index = 0;
                for (size_t k = 0; k < i; k++) {
                    index += q[k];
                }
                index += j - p[i];
                return coefs[index];
            }
            else{
                return T(0); 
            }
        } else {
            // Exploite la symétrie de la matrice
            if (i >= p[j]){
                size_t index = 0;
                for (size_t k = 0; k < j; k++) {
                    index += q[k];
                }
                index += i - p[j];
                return coefs[index];
            }
            else{
                return T(0); 
            }
        }
    }
    
    // Accès écriture (Nécessaire pour remplir la matrice)
    T& operator()(int i, int j) {
        if (i < j) std::swap(i, j); // Symétrie : on se ramène au triangle inférieur
        
        if (j < static_cast<int>(p[i])) {
            throw runtime_error("Erreur: Tentative d'ecriture hors du profil.");
        }

        size_t index = 0;
        for (size_t k = 0; k < i; k++) {
            index += q[k];
        }
        index += j - p[i];
        return coefs[index];
    }

    // Opérateurs

    vector<T> operator*(const vector<T>& x) const {
        if (x.size() != static_cast<size_t>(n)) {
            // dimensions
            throw invalid_argument("Matrix and vector dimensions do not match for multiplication.");
        }
        
        vector<T> res(n, T(0));
        size_t index_offset = 0;

        for (size_t i = 0; i < n; i++) {
            for (size_t j = p[i]; j < i; j++) {
                const T& val = coefs[index_offset + j - p[i]];
                res[i] += val * x[j];
                res[j] += val * x[i];
            }
            const T& diag_val = coefs[index_offset + i - p[i]];
            res[i] += diag_val * x[i];

            index_offset += q[i];
        }
        return res;
    }

    void operator+=(const ProfileMatrix<T>& M){
        for (size_t i = 0;i<coefs.size();i++){
            coefs[i] += M.coefs[i];
        }
    }

    void operator-=(const ProfileMatrix<T>& M){
        for (size_t i = 0;i<coefs.size();i++){
            coefs[i] -= M.coefs[i];
        }
    }

    void solve(vector<T>& x, const vector<T>& b) {
        if (b.size() != static_cast<size_t>(n)) throw invalid_argument("Dimension mismatch");
        
        // Copie pour LDL^T (on ne modifie pas la matrice d'origine)
        ProfileMatrix<T> A_ldl(*this);
        
        // Pré-calcul des offsets pour accès rapide
        vector<size_t> offsets(n);
        size_t current = 0;
        for(int i=0; i<n; ++i) {
            offsets[i] = current;
            current += A_ldl.q[i];
        }

        // Factorisation LDL^T
        for (int i = 0; i < n; ++i) {
            T d_val = T(0);
            int pi = static_cast<int>(A_ldl.p[i]);
            
            for (int j = pi; j < i; ++j) {
                T sum = T(0);
                int pj = static_cast<int>(A_ldl.p[j]);
                int start_k = (pi > pj) ? pi : pj;

                for (int k = start_k; k < j; ++k) {
                    int pk = static_cast<int>(A_ldl.p[k]);
                    T L_ik = A_ldl.coefs[offsets[i] + k - pi];
                    T L_jk = A_ldl.coefs[offsets[j] + k - pj];
                    T D_kk = A_ldl.coefs[offsets[k] + k - pk];
                    sum += L_ik * L_jk * D_kk;
                }
                
                size_t idx_ij = offsets[i] + j - pi;
                size_t idx_jj = offsets[j] + j - pj;
                
                T L_ij = (A_ldl.coefs[idx_ij] - sum) / A_ldl.coefs[idx_jj];
                A_ldl.coefs[idx_ij] = L_ij;
                d_val += L_ij * L_ij * A_ldl.coefs[idx_jj];
            }
            
            size_t idx_ii = offsets[i] + i - pi;
            A_ldl.coefs[idx_ii] -= d_val;
            if (abs(A_ldl.coefs[idx_ii]) < 1e-14) throw runtime_error("Pivot nul");
        }

        // Résolution
        x = b;
        // Descente L z = b
        for (int i = 0; i < n; ++i) {
            T sum = T(0);
            int pi = static_cast<int>(A_ldl.p[i]);
            for (int j = pi; j < i; ++j) 
                sum += A_ldl.coefs[offsets[i] + j - pi] * x[j];
            x[i] -= sum;
        }
        // Diagonale D y = z
        for (int i = 0; i < n; ++i) 
            x[i] /= A_ldl.coefs[offsets[i] + i - static_cast<int>(A_ldl.p[i])];
        // Remontée L^T x = y
        for (int i = n - 1; i >= 0; --i) {
            int pi = static_cast<int>(A_ldl.p[i]);
            for (int j = pi; j < i; ++j) 
                x[j] -= A_ldl.coefs[offsets[i] + j - pi] * x[i];
        }
    }
};

template<typename T>

ProfileMatrix<T> operator+(const ProfileMatrix<T>& A,const ProfileMatrix<T>& B){
    ProfileMatrix<T> R(A);
    R += B;
    return R;
}

template<typename T>

ProfileMatrix<T> operator-(const ProfileMatrix<T>& A,const ProfileMatrix<T>& B){
    ProfileMatrix<T> R(A);
    R -= B;
    return R;
}

#endif //MATH_HPP
