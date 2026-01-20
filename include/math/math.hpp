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
//     classe abstraite Matrix
//---------------------------------------------------------------------------

template <typename T>

class Matrix {

protected:

    vector<T> mat;
    int n_rows;
    int n_cols;

public:

    Matrix(int rows, int cols) : n_rows(rows), n_cols(cols), mat(rows * cols, T(0)) {}
    
    // Destructeur

    virtual ~Matrix() {} 

    int rows() const { return n_rows; }
    int cols() const { return n_cols; }

    // Méthodes Virtuelles Pures (pas de méthodes d'accés pour la performance)

    virtual vector<T> operator*(const vector<T>& x) const = 0;

    virtual void solve(vector<T>& x, const vector<T>& b) = 0;
    
    virtual void print(ostream& os) const = 0;
};


//---------------------------------------------------------------------------
//     classe FullMatrix
//---------------------------------------------------------------------------

template <typename T>

class FullMatrix : public Matrix<T> {

public:

    FullMatrix(int n, int m) : Matrix<T>(n, m) {}

    // Accès

    T& operator()(int i, int j) { return data[i * this->n_cols + j]; }
    const T& operator()(int i, int j) const { return data[i * this->n_cols + j]; }

    // Produit Matrice-Vecteur
    
    vector<T> operator*(const vector<T>& x) const override {
        vector<T> res(this->n_rows, T(0));
        for(int i = 0; i < this->n_rows; ++i) {
            for(int j = 0; j < this->n_cols; ++j) {
                res[i] += (*this)(i, j) * x[j];
            }
        }
        return res;
    }

    // Résolution via pivot de Gauss (LU)

    void solve(vector<T>& x, const vector<T>& b) override {

        // ICI : Implémenter l'élimination de Gauss
        // C'est requis par le sujet pour les matrices pleines

    }

    void print(ostream& os) const override {
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
class ProfileMatrix : public Matrix<T> {
private:

    vector<T> p;
    vector<T> q;

public:

    ProfileMatrix(int n, int m) : Matrix<T>(n, m), p(n,int(0)), q(n,int(0)) {}

    //remplissage de p et q (si on a la matrice ou si on sait qu'on fait de la FEM et qu'on a les triangles)

    // Accès

    T& operator()(int i, int j) {/* ... */}

    vector<T> operator*(const vector<T>& x) const override {
        //
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
