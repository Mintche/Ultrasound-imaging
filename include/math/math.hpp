#ifndef MATH_HPP
#define MATH_HPP
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <cmath>
using namespace std;
//---------------------------------------------------------------------------
//  vector<T> algebra
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
//     classe Pint ( pair<int,int>)
//---------------------------------------------------------------------------

typedef pair<int,int> Pint;

inline bool operator<(const Pint& ij1, const Pint& ij2)
{
  if ( ij1.first < ij2.first) return true;
  if ( ij1.first == ij2.first and ij1.second < ij2.second ) return true;
  return false;
}
inline ostream& operator<<(ostream& os,const Pint& ij)
{
    os << "(" << ij.first << "," << ij.second << ")";
    return os;
}

//---------------------------------------------------------------------------
//     classe Sparse<T> héritée de map<Pint,T>
//---------------------------------------------------------------------------

template <typename T>
class Sparse : public map<Pint,T>
{public:
  int m,n; //dimensions de la matrice
  Sparse(int mi=0,int ni=0);
  T& operator()(int i, int j);
  T operator()(int i, int j) const;
  void supprime(int i,int j);
  Sparse<T>& operator+=(const Sparse<T>& M);  // +=M
  Sparse<T>& operator-=(const Sparse<T>& M);  // -=M
  Sparse<T>& operator*=(const T& s);          // *=s
  Sparse<T>& operator/=(const T& s);          // /=s
  void remplissage() const; // affiche le remplissage
};

// implémentations des fonctions associées à la classe Sparse

template <typename T>

Sparse<T>::Sparse(int mi, int ni) //il herite du constructeur de map donc il initialise la map aussi
{
    m = mi;
    n = ni;
}

template<typename T>
T Sparse<T>::operator()(int i, int j) const
{
    if (i <= m && j <= n && i > 0 && j > 0) {
        auto it = this->find(Pint(i, j)); // Utiliser find
        if (it != this->end()) return it->second;
    }
    return T();
}

template<typename T>

T& Sparse<T>::operator()(int i, int j)
{
    if (i <= m && j <= n && i > 0 && j > 0) return (*this)[Pint(i,j)];
    if ( i > m && j > n){
        m = i;
        n = j;
        return (*this)[Pint(i,j)];
    }
    if (j > n){
        n = j;
        return (*this)[Pint(i,j)];
    }
    if (i > m){
        m = i;
        return (*this)[Pint(i,j)];
    }
    return (*this)[Pint(i,j)];
}

template<typename T>
ostream& operator<<(ostream& os,const Sparse<T>& S)
{
    os << "Matrcie Sparse de " << S.size() << " éléments." << endl;

    for (auto& it : S){
        os << it.first << " : " << it.second << endl;
    }
    return os;
}

template<typename T>

void Sparse<T>::remplissage() const
{
    cout << "  ";
    for(int i =1 ; i<=m; i++) cout << i << " ";
    cout << endl;

    for(int i =1 ; i<=m; i++){
        cout << i << " ";
        for (int j=1; j<=n; j++){
            auto it = (*this).find(Pint(i,j));
            if ( it != (*this).end()) cout << "x ";
            else cout << "  ";
        }
        cout << endl;
    }
}

template<typename T>

Sparse<T>& Sparse<T>::operator+=(const Sparse<T>& M)
{
    for (auto& it : M){
        (*this)(it.first.first,it.first.second) += it.second;
    }
    return *this;
}

template<typename T>

Sparse<T>& Sparse<T>::operator-=(const Sparse<T>& M)
{
    for (auto& it : M){
        (*this)(it.first.first,it.first.second) -= it.second;
    }
    return *this;
}

template<typename T>

Sparse<T>& Sparse<T>::operator*=(const T& s)
{
    for (auto& it : *this){
        it.second *= s;
    }
    return *this;
}

template<typename T>

Sparse<T>& Sparse<T>::operator/=(const T& s)
{
    for (auto& it : *this){
        it.second /= s;
    }
    return *this;
}

template<typename T>

Sparse<T> operator+(const Sparse<T>& A,const Sparse<T>& B)
{
    Sparse<T> R(A);
    return R+=B;
}

template<typename T>

Sparse<T> operator-(const Sparse<T>& A,const Sparse<T>& B)
{
    Sparse<T> R(A);
    return R-=B;
}

template<typename T>

Sparse<T> operator*(const Sparse<T>& A,const T& s)
{
    Sparse<T> R(A);
    return R*=s;
}

template<typename T>

Sparse<T> operator*(const T& s,const Sparse<T>& A)
{
    return A*s;
}

template<typename T>

Sparse<T> operator/(const Sparse<T>& A,const T& s)
{
    Sparse<T> R(A);
    return R/=s;
}

template<typename T>
vector<T> operator*(const Sparse<T>& A, const vector<T>& b)
{
    if (A.n != b.size()) {
        cerr << "Erreur dimensions: A(" << A.m << "x" << A.n 
             << ") * b(" << b.size() << ")" << endl;
        return vector<T>(); 
    }
    vector<T> x(A.m, T(0)); 
    for (auto const& element : A) {
        int i = element.first.first;
        int j = element.first.second;
        T val = element.second;

        x[i-1] += val * b[j-1];
    }

    return x;
}

template<typename T>

vector<T> gradConj(const Sparse<T>& A, const vector <T>& b , const vector <T>& x0 ,int kmax , const T& eps=0.0001)
{
    vector<T> x = x0;
    vector<T> G = A*x0 - b;
    vector<T> W = G;
    for (int k=0;k<kmax;k++){
        T theta = (G|W)/(A*W|W);
        x = x - theta*W;
        vector<T> temp = G;
        G = G - theta*A*W;
        if (norme(G)<=eps*norme(b)) cout << "Convegence à " << eps <<" près"; return x;
        W= G + (norme(G)/norme(temp))*W;
    }
    return x; 
}
#endif //MATH_HPP
