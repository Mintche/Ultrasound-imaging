#ifndef MATH_HPP
#define MATH_HPP
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <cmath>
using namespace std;
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
    T s=0.;
    auto itv = v.begin();
    for(auto itu=u.begin(); itu!=u.end(); ++itu, ++itv) s+=*itu * *itv;
    return s;
}
template<typename T> T norme(const vector<T>&u)
{
    return sqrt(u|u);
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

    T& operator()(int i, int j) { return coefs[i * this->n_cols + j]; }
    const T& operator()(int i, int j) const { return coefs[i * this->n_cols + j]; }

    // Produit Matrice-Vecteur
    
    vector<T> operator*(const vector<T>& x) const {
        vector<T> res(this->n_rows, T(0));
        for(int i = 0; i < this->n_rows; ++i) {
            for(int j = 0; j < this->n_cols; ++j) {
                res[i] += (*this)(i, j) * x[j];
            }
        }
        return res;
    }

    // Résolution via pivot de Gauss (LU)

    void solve(vector<T>& x, const vector<T>& b){

        int n = this->n_rows;
        FullMatrix<T> A = *this; // Copie de la matrice
        x = b; // Initialisation du vecteur solution

        // Elimination de Gauss
        for (int k = 0; k < n; ++k) {
            for (int i = k + 1; i < n; ++i) {
                T factor = A(i, k) / A(k, k);
                for (int j = k; j < n; ++j) {
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
        for (int i = 0; i < this->n_rows; i++){
            for (int j = 0; j < this->n_cols; j++){
                os << (*this)(i,j) << " ";
            }
            os <<endl;
        }
        os << ")";
    }
};

//---------------------------------------------------------------------------
//     Classe ProfileMatrix
//---------------------------------------------------------------------------

template <typename T>
class ProfileMatrix {
private:

    vector<T> coefs;
    vector<T> p;
    vector<T> q;

public:

    ProfileMatrix(vector<T>& p_in, int n_rows) : p(p_in), q(p_in.size(), size_t(0)) {

        // Initialisation de q
        q[0] = 0;

        for (size_t i = 1; i < p.size(); i++) {
            q[i] = n_rows - p[i] + 1;
        }

        // Initialisation de coefs
        size_t total_size = 0;
        for (size_t i = 0; i < p.size(); i++) {
            total_size += q[i];
        }
        coefs.resize(total_size, T(0));
    }

    // Accès

    T& operator()(int i, int j) {
        if (i >= j) {
            size_t index = 0;
            for (size_t k = 0; k < i; k++) {
                index += q[k];
            }
            index += (j - (i - p[i]));
            return coefs[index];
        } else {
            size_t index = 0;
            for (size_t k = 0; k < j; k++) {
                index += q[k];
            }
            index += (i - (j - p[j]));
            return coefs[index];
        }
    }

    vector<T> operator*(const vector<T>& x) const override {
        vector<T> res(this->n_rows, T(0));
        for (int i = 0; i < this->n_rows; ++i) {
            for (int j = 0; j < this->n_cols; ++j) {
                res[i] += (*this)(i, j) * x[j];
            }
        }
        return res;
    }

    void solve(vector<T>& x, const vector<T>& b) override {
        // ATTENTION: Helmholtz (k0^2) n'est pas définie positive. Gauss est plus sûr.
    }
    
    void print(ostream& os) const override { /* ... */ }
};

//---------------------------------------------------------------------------
//     Solver
//---------------------------------------------------------------------------


#endif //MATH_HPP
