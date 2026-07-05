#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace std;
using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_dec_float_100;
using Float = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>>;

cpp_int parse_int(const string &s){
    cpp_int x=0; int i=0; bool neg=false; if(s.size() && s[0]=='-'){neg=true;i=1;}
    for(; i<(int)s.size(); ++i){ char c=s[i]; if(c>='0'&&c<='9'){ x *= 10; x += c-'0'; }}
    return neg? -x:x;
}

Float toF(const cpp_int &x){ return x.convert_to<Float>(); }

cpp_int roundF(Float x){
    bool neg=false; if(x<0){ neg=true; x=-x; }
    Float y = floor(x + Float("0.5"));
    // convert_to<cpp_int> should truncate toward zero
    string s = y.str(0, std::ios_base::fixed);
    cpp_int r = parse_int(s);
    return neg ? -r : r;
}

Float dotBF(const vector<cpp_int>& a, const vector<Float>& b){
    Float s=0;
    int n=a.size();
    for(int i=0;i<n;i++) if(a[i]!=0) s += toF(a[i]) * b[i];
    return s;
}
Float dotFF(const vector<Float>& a, const vector<Float>& b){
    Float s=0; int n=a.size(); for(int i=0;i<n;i++) s += a[i]*b[i]; return s;
}

void compute_gso(const vector<vector<cpp_int>>& B, vector<vector<Float>>& mu, vector<vector<Float>>& bs, vector<Float>& norm){
    int m=B.size(), n=B[0].size();
    mu.assign(m, vector<Float>(m));
    bs.assign(m, vector<Float>(n));
    norm.assign(m, Float(0));
    for(int i=0;i<m;i++){
        for(int c=0;c<n;c++) bs[i][c]=toF(B[i][c]);
        for(int j=0;j<i;j++){
            if(norm[j]==0){ mu[i][j]=0; continue; }
            Float d = dotBF(B[i], bs[j]);
            mu[i][j] = d / norm[j];
            Float mij = mu[i][j];
            if(mij != 0){
                for(int c=0;c<n;c++) bs[i][c] -= mij * bs[j][c];
            }
        }
        norm[i]=dotFF(bs[i],bs[i]);
    }
}

int bitlen_abs(const cpp_int& x){
    if(x==0) return 0; cpp_int y=x<0?-x:x; return boost::multiprecision::msb(y)+1;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    Float delta = Float("0.99");
    if(argc>=2) delta = Float(argv[1]);
    int m,n; if(!(cin>>m>>n)) return 1;
    vector<vector<cpp_int>> B(m, vector<cpp_int>(n));
    string s;
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){ cin>>s; B[i][j]=parse_int(s); }
    vector<vector<Float>> mu, bs; vector<Float> norm;
    compute_gso(B,mu,bs,norm);
    int k=1; long long iter=0, swaps=0, reds=0;
    while(k<m){
        iter++;
        // size reduction, repeatedly with recompute after each row change
        for(int j=k-1;j>=0;j--){
            cpp_int r = roundF(mu[k][j]);
            if(r != 0){
                for(int c=0;c<n;c++) B[k][c] -= r*B[j][c];
                reds++;
                compute_gso(B,mu,bs,norm);
            }
        }
        Float lhs = norm[k];
        Float rhs = (delta - mu[k][k-1]*mu[k][k-1]) * norm[k-1];
        if(lhs >= rhs){
            k++;
        }else{
            swap(B[k],B[k-1]); swaps++;
            compute_gso(B,mu,bs,norm);
            k = max(k-1,1);
        }
        if(iter % 10000 == 0){ cerr << "iter "<<iter<<" k "<<k<<" swaps "<<swaps<<" reds "<<reds<<"\n"; }
        if(iter > 2000000){ cerr << "too many iterations\n"; break; }
    }
    cerr << "done iter "<<iter<<" swaps "<<swaps<<" reds "<<reds<<"\n";
    cout << m << " " << n << "\n";
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            if(j) cout << ' ';
            cout << B[i][j];
        }
        cout << "\n";
    }
    return 0;
}
